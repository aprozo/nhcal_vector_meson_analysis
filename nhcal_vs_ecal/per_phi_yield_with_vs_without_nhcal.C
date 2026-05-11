#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>

const int PDG_PHI = 333;
const int PDG_KP  = 321;


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

// for one reco-particle index, fill flags telling which backward calo(s) saw it
void compute_calo_match_flags(int reco_index,
                              const TTreeReaderArray<unsigned int>& reco_tracks_begin,
                              const TTreeReaderArray<unsigned int>& reco_tracks_end,
                              const TTreeReaderArray<int>&          reco_tracks_idx,
                              const TTreeReaderArray<int>&          proj_track_idx,
                              const TTreeReaderArray<int>&          nhcal_match_proj,
                              const TTreeReaderArray<int>&          ecaln_match_proj,
                              bool& kaon_in_nhcal,
                              bool& kaon_in_ecaln)
{
    kaon_in_nhcal = false;
    kaon_in_ecaln = false;
    unsigned int tracks_begin_for_reco = reco_tracks_begin[reco_index];
    unsigned int tracks_end_for_reco   = reco_tracks_end[reco_index];
    if (tracks_end_for_reco <= tracks_begin_for_reco) return;
    int track_index = reco_tracks_idx[tracks_begin_for_reco];
    kaon_in_nhcal = track_in_match_list(track_index, nhcal_match_proj, proj_track_idx);
    kaon_in_ecaln = track_in_match_list(track_index, ecaln_match_proj, proj_track_idx);
}


// ---------------------------------------------------------------------------
// MAIN MACRO
// ---------------------------------------------------------------------------
void per_phi_yield_with_vs_without_nhcal(const char* file_path = "") {

    if (file_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'per_phi_yield_with_vs_without_nhcal.C(\"FILE.root\")'\n";
        return;
    }

    // ---------- 1. open the file -----------------------------------------
    TFile* root_file = TFile::Open(file_path);
    if (!root_file || root_file->IsZombie()) { std::cerr << "cannot open\n"; return; }

    // ---------- 2. bind every needed branch ------------------------------
    TTreeReader reader("events", root_file);
    TTreeReaderArray<int>           mc_pdg            (reader, "MCParticles.PDG");
    TTreeReaderArray<int>           mc_gen_status     (reader, "MCParticles.generatorStatus");
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

    // ---------- 3. counters for the two cut definitions ------------------
    long n_phi_with_both_reco  = 0;
    long n_phi_ecaln_only_cut  = 0;    // both Ks ECalN-matched (no-nHcal baseline)
    long n_phi_either_calo_cut = 0;    // each K in either calo  (with nHcal)

    // ---------- 4. main event loop ---------------------------------------
    while (reader.Next()) {

        // (a) find the truth K+ K- pair
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;

        // (b) walk MC -> reco for both kaons; require both to be reco
        int kaon_plus_reco_index  = reco_index_of(kaon_plus_mc_index,  assoc_sim_id, assoc_rec_id);
        int kaon_minus_reco_index = reco_index_of(kaon_minus_mc_index, assoc_sim_id, assoc_rec_id);
        if (kaon_plus_reco_index < 0 || kaon_minus_reco_index < 0) continue;
        ++n_phi_with_both_reco;

        // (c) compute the four per-kaon calo-match flags
        bool kaon_plus_in_nhcal,  kaon_plus_in_ecaln;
        bool kaon_minus_in_nhcal, kaon_minus_in_ecaln;
        compute_calo_match_flags(kaon_plus_reco_index,
            reco_tracks_begin, reco_tracks_end, reco_tracks_idx,
            proj_track_idx, nhcal_match_proj, ecaln_match_proj,
            kaon_plus_in_nhcal, kaon_plus_in_ecaln);
        compute_calo_match_flags(kaon_minus_reco_index,
            reco_tracks_begin, reco_tracks_end, reco_tracks_idx,
            proj_track_idx, nhcal_match_proj, ecaln_match_proj,
            kaon_minus_in_nhcal, kaon_minus_in_ecaln);

        // (d) apply the two cuts and bump the appropriate counters
        if (kaon_plus_in_ecaln && kaon_minus_in_ecaln)
            ++n_phi_ecaln_only_cut;
        if ((kaon_plus_in_nhcal  || kaon_plus_in_ecaln)
         && (kaon_minus_in_nhcal || kaon_minus_in_ecaln))
            ++n_phi_either_calo_cut;
    }

    // ---------- 5. report + plot -----------------------------------------
    long gain_abs = n_phi_either_calo_cut - n_phi_ecaln_only_cut;
    double gain_rel = n_phi_ecaln_only_cut > 0 ? (double) gain_abs / n_phi_ecaln_only_cut : 0.0;

    std::cout << "phi -> KK with both kaons reco     : " << n_phi_with_both_reco << "\n";
    std::cout << "  without nHcal (both Ks ECalN)    : " << n_phi_ecaln_only_cut << "\n";
    std::cout << "  with    nHcal (each K in either) : " << n_phi_either_calo_cut << "\n";
    std::cout << "  extra kept thanks to nHcal       : " << gain_abs
              << "  (+" << 100*gain_rel << "%)\n";

    TH1F* h = new TH1F("h",
        Form("Adding nHcal to the ECalN baseline keeps +%.0f%% more phi;;# of phi -> KK kept",
             100 * gain_rel), 2, 0, 2);
    h->SetBinContent(1, n_phi_ecaln_only_cut);
    h->SetBinContent(2, n_phi_either_calo_cut);
    h->GetXaxis()->SetBinLabel(1, "ECalN only (no nHcal)");
    h->GetXaxis()->SetBinLabel(2, "ECalN U nHcal (with nHcal)");
    h->SetFillColor(kAzure - 4);
    h->SetMinimum(0);

    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    h->Draw("BAR");
    TLatex label; label.SetTextAlign(22); label.SetTextSize(0.04);
    label.DrawLatex(h->GetBinCenter(1), n_phi_ecaln_only_cut  + 0.03 * n_phi_either_calo_cut,
                    Form("%ld", n_phi_ecaln_only_cut));
    label.DrawLatex(h->GetBinCenter(2), n_phi_either_calo_cut + 0.03 * n_phi_either_calo_cut,
                    Form("%ld", n_phi_either_calo_cut));
    canvas->SaveAs("per_phi_yield_with_vs_without_nhcal.pdf");
    std::cout << "wrote per_phi_yield_with_vs_without_nhcal.pdf\n";

    root_file->Close();
}
