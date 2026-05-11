#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>

const int    PDG_PHI       = 333;
const int    PDG_KP        = 321;
const double NHCAL_ETA_MIN = -4.05;
const double NHCAL_ETA_MAX = -1.20;


// pseudorapidity from 3-momentum
double pseudorapidity(double px, double py, double pz) {
    double momentum_magnitude = std::sqrt(px*px + py*py + pz*pz);
    if (momentum_magnitude == 0 || std::fabs(pz) >= momentum_magnitude) return std::nan("");
    return 0.5 * std::log((momentum_magnitude + pz) / (momentum_magnitude - pz));
}

// find the K+ and K- MC indices for a phi -> K+ K- decay
bool find_phi_to_kaons(const TTreeReaderArray<int>&          mc_pdg,
                       const TTreeReaderArray<int>&          mc_gen_status,
                       const TTreeReaderArray<unsigned int>& mc_daughters_begin,
                       const TTreeReaderArray<unsigned int>& mc_daughters_end,
                       const TTreeReaderArray<int>&          mc_daughter_index,
                       int& kaon_plus_mc_index,
                       int& kaon_minus_mc_index)
{
    kaon_plus_mc_index  = -1;
    kaon_minus_mc_index = -1;
    for (size_t mc_particle_index = 0; mc_particle_index < mc_pdg.GetSize(); ++mc_particle_index) {
        if (mc_pdg[mc_particle_index] != PDG_PHI) continue;
        int gen_status_value = mc_gen_status[mc_particle_index];
        if (gen_status_value != 1 && gen_status_value != 2) continue;
        for (unsigned int daughter_slot = mc_daughters_begin[mc_particle_index];
                          daughter_slot < mc_daughters_end[mc_particle_index];
                          ++daughter_slot) {
            int child_mc_index = mc_daughter_index[daughter_slot];
            if (mc_pdg[child_mc_index] ==  PDG_KP) kaon_plus_mc_index  = child_mc_index;
            if (mc_pdg[child_mc_index] == -PDG_KP) kaon_minus_mc_index = child_mc_index;
        }
        if (kaon_plus_mc_index >= 0 && kaon_minus_mc_index >= 0) return true;
    }
    return false;
}

// MC index -> reco-particle index, -1 if unmatched
int reco_index_of(int mc_index,
                  const TTreeReaderArray<int>&          assoc_sim_id,
                  const TTreeReaderArray<int>&          assoc_rec_id)
{
    for (size_t assoc_slot = 0; assoc_slot < assoc_sim_id.GetSize(); ++assoc_slot)
        if (assoc_sim_id[assoc_slot] == mc_index) return assoc_rec_id[assoc_slot];
    return -1;
}

// does the given track index appear in match_proj (via the projection table)?
bool track_in_match_list(int track_index,
                         const TTreeReaderArray<int>& match_proj,
                         const TTreeReaderArray<int>& proj_track_idx)
{
    for (size_t match_slot = 0; match_slot < match_proj.GetSize(); ++match_slot) {
        int projection_index = match_proj[match_slot];
        if (projection_index < 0 || (size_t) projection_index >= proj_track_idx.GetSize()) continue;
        if (proj_track_idx[projection_index] == track_index) return true;
    }
    return false;
}

// percent of `numerator` against `denominator`, safe against division by zero
double pct(long numerator, long denominator) {
    return denominator > 0 ? 100.0 * numerator / denominator : 0.0;
}


// ---------------------------------------------------------------------------
// MAIN MACRO
// ---------------------------------------------------------------------------
void per_kaon_match_overlap(const char* file_path = "") {

    if (file_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'per_kaon_match_overlap.C(\"FILE.root\")'\n";
        return;
    }

    // ---------- 1. open the file -----------------------------------------
    TFile* root_file = TFile::Open(file_path);
    if (!root_file || root_file->IsZombie()) { std::cerr << "cannot open\n"; return; }

    // ---------- 2. bind every needed branch ------------------------------
    TTreeReader reader("events", root_file);
    TTreeReaderArray<int>           mc_pdg            (reader, "MCParticles.PDG");
    TTreeReaderArray<int>           mc_gen_status     (reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double>        mc_px             (reader, "MCParticles.momentum.x");
    TTreeReaderArray<double>        mc_py             (reader, "MCParticles.momentum.y");
    TTreeReaderArray<double>        mc_pz             (reader, "MCParticles.momentum.z");
    TTreeReaderArray<unsigned int>  mc_daughters_begin(reader, "MCParticles.daughters_begin");
    TTreeReaderArray<unsigned int>  mc_daughters_end  (reader, "MCParticles.daughters_end");
    TTreeReaderArray<int>           mc_daughter_index (reader, "_MCParticles_daughters.index");
    TTreeReaderArray<int>           assoc_sim_id      (reader, "_ReconstructedChargedParticleAssociations_sim.index");
    TTreeReaderArray<int>           assoc_rec_id      (reader, "_ReconstructedChargedParticleAssociations_rec.index");
    TTreeReaderArray<unsigned int>  reco_tracks_begin (reader, "ReconstructedChargedParticles.tracks_begin");
    TTreeReaderArray<unsigned int>  reco_tracks_end   (reader, "ReconstructedChargedParticles.tracks_end");
    TTreeReaderArray<int>           reco_tracks_idx   (reader, "_ReconstructedChargedParticles_tracks.index");
    TTreeReaderArray<int>           proj_track_idx    (reader, "_CalorimeterTrackProjections_track.index");
    TTreeReaderArray<int>           nhcal_match_proj  (reader, "_HcalEndcapNTrackClusterMatches_track.index");
    TTreeReaderArray<int>           ecaln_match_proj  (reader, "_EcalEndcapNTrackClusterMatches_track.index");

    std::cout << "opened " << file_path << " (" << reader.GetEntries(true) << " events)\n";

    // ---------- 3. four mutually-exclusive Venn counters -----------------
    long n_kaons_in_eta_and_reco = 0;
    long n_nhcal_only            = 0;
    long n_ecaln_only            = 0;
    long n_in_both_calos         = 0;
    long n_in_neither_calo       = 0;

    // ---------- 4. main event loop ---------------------------------------
    while (reader.Next()) {

        // (a) find the truth K+ K- pair
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;
        int kaon_mc_indices[2] = {kaon_plus_mc_index, kaon_minus_mc_index};

        // (b) for each kaon, classify by (nHcal? ECalN?)
        for (int slot = 0; slot < 2; ++slot) {
            int mc_kaon_index = kaon_mc_indices[slot];

            // (b.i) truth eta cut
            double kaon_truth_eta = pseudorapidity(mc_px[mc_kaon_index],
                                                   mc_py[mc_kaon_index],
                                                   mc_pz[mc_kaon_index]);
            if (!((NHCAL_ETA_MIN <= kaon_truth_eta) && (kaon_truth_eta < NHCAL_ETA_MAX))) continue;

            // (b.ii) reco cut
            int kaon_reco_index = reco_index_of(mc_kaon_index, assoc_sim_id, assoc_rec_id);
            if (kaon_reco_index < 0) continue;
            ++n_kaons_in_eta_and_reco;

            // (b.iii) get the underlying track index
            unsigned int tracks_begin_for_reco = reco_tracks_begin[kaon_reco_index];
            unsigned int tracks_end_for_reco   = reco_tracks_end[kaon_reco_index];
            if (tracks_end_for_reco <= tracks_begin_for_reco) continue;
            int track_index = reco_tracks_idx[tracks_begin_for_reco];

            // (b.iv) which calorimeters claim this track?
            bool kaon_in_nhcal = track_in_match_list(track_index, nhcal_match_proj, proj_track_idx);
            bool kaon_in_ecaln = track_in_match_list(track_index, ecaln_match_proj, proj_track_idx);
            if (kaon_in_nhcal && kaon_in_ecaln)             ++n_in_both_calos;
            else if (kaon_in_nhcal)                         ++n_nhcal_only;
            else if (kaon_in_ecaln)                         ++n_ecaln_only;
            else                                            ++n_in_neither_calo;
        }
    }

    // ---------- 5. report + plot -----------------------------------------
    long denom = n_kaons_in_eta_and_reco;
    std::cout << "kaons in eta AND reco : " << n_kaons_in_eta_and_reco << "\n";
    std::cout << "  nHcal only          : " << n_nhcal_only      << "  (" << pct(n_nhcal_only,      denom) << "%)\n";
    std::cout << "  ECalN only          : " << n_ecaln_only      << "  (" << pct(n_ecaln_only,      denom) << "%)\n";
    std::cout << "  both calos          : " << n_in_both_calos   << "  (" << pct(n_in_both_calos,   denom) << "%)\n";
    std::cout << "  neither             : " << n_in_neither_calo << "  (" << pct(n_in_neither_calo, denom) << "%)\n";

    TH1F* h = new TH1F("h", Form("Per-kaon nHcal-vs-ECalN overlap  (N=%ld);;%% of kaons in eta AND reco",
                                 n_kaons_in_eta_and_reco), 4, 0, 4);
    h->SetBinContent(1, pct(n_nhcal_only,      denom));
    h->SetBinContent(2, pct(n_ecaln_only,      denom));
    h->SetBinContent(3, pct(n_in_both_calos,   denom));
    h->SetBinContent(4, pct(n_in_neither_calo, denom));
    h->GetXaxis()->SetBinLabel(1, "nHcal only");
    h->GetXaxis()->SetBinLabel(2, "ECalN only");
    h->GetXaxis()->SetBinLabel(3, "both calos");
    h->GetXaxis()->SetBinLabel(4, "neither");
    h->SetFillColor(kAzure - 4);
    h->SetMinimum(0);

    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    h->Draw("BAR");
    TLatex label; label.SetTextAlign(22); label.SetTextSize(0.04);
    for (int i = 1; i <= 4; ++i)
        label.DrawLatex(h->GetBinCenter(i), h->GetBinContent(i) + 1.5,
                        Form("%.1f%%", h->GetBinContent(i)));
    canvas->SaveAs("per_kaon_match_overlap.pdf");
    std::cout << "wrote per_kaon_match_overlap.pdf\n";

    root_file->Close();
}
