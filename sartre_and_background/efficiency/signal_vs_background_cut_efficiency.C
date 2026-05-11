
#include <iostream>
#include <vector>
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

// reco index -> parent track has an entry in `match_proj`?
bool is_track_in_match_list(int reco_index,
                            const TTreeReaderArray<unsigned int>& reco_tracks_begin,
                            const TTreeReaderArray<unsigned int>& reco_tracks_end,
                            const TTreeReaderArray<int>&          reco_tracks_idx,
                            const TTreeReaderArray<int>&          proj_track_idx,
                            const TTreeReaderArray<int>&          match_proj)
{
    if (reco_index < 0) return false;
    unsigned int tracks_begin_for_reco = reco_tracks_begin[reco_index];
    unsigned int tracks_end_for_reco   = reco_tracks_end[reco_index];
    if (tracks_end_for_reco <= tracks_begin_for_reco) return false;
    int track_index = reco_tracks_idx[tracks_begin_for_reco];
    for (size_t match_slot = 0; match_slot < match_proj.GetSize(); ++match_slot) {
        int projection_index = match_proj[match_slot];
        if (projection_index < 0 || (size_t) projection_index >= proj_track_idx.GetSize()) continue;
        if (proj_track_idx[projection_index] == track_index) return true;
    }
    return false;
}


// ---------------------------------------------------------------------------
// SIGNAL processing: for every truth phi -> KK, check the cut on both reco kaons
// ---------------------------------------------------------------------------
void process_signal(const char* signal_path, long& n_signal_total, long& n_signal_kept) {

    std::cout << "signal: " << signal_path << "\n";
    TFile* signal_file = TFile::Open(signal_path);
    TTreeReader reader("events", signal_file);

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

    while (reader.Next()) {

        // (a) find the K+ K- pair
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;
        ++n_signal_total;

        // (b) walk MC -> reco
        int kaon_plus_reco_index  = reco_index_of(kaon_plus_mc_index,  assoc_sim_id, assoc_rec_id);
        int kaon_minus_reco_index = reco_index_of(kaon_minus_mc_index, assoc_sim_id, assoc_rec_id);
        if (kaon_plus_reco_index < 0 || kaon_minus_reco_index < 0) continue;

        // (c) apply the 'both kaons nHcal-matched' cut
        if (is_track_in_match_list(kaon_plus_reco_index,  reco_tracks_begin, reco_tracks_end,
                                   reco_tracks_idx, proj_track_idx, nhcal_match_proj)
         && is_track_in_match_list(kaon_minus_reco_index, reco_tracks_begin, reco_tracks_end,
                                   reco_tracks_idx, proj_track_idx, nhcal_match_proj))
            ++n_signal_kept;
    }
    std::cout << "signal: kept " << n_signal_kept << "/" << n_signal_total << "\n";
    signal_file->Close();
}


// ---------------------------------------------------------------------------
// BACKGROUND processing: every OS reco pair, apply same cut
// ---------------------------------------------------------------------------
void process_background(const char* bg_path, long& n_bg_total, long& n_bg_kept) {

    std::cout << "bg:     " << bg_path << "\n";
    TFile* bg_file = TFile::Open(bg_path);
    TTreeReader reader("events", bg_file);

    TTreeReaderArray<float>         reco_charge      (reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<unsigned int>  reco_tracks_begin(reader, "ReconstructedChargedParticles.tracks_begin");
    TTreeReaderArray<unsigned int>  reco_tracks_end  (reader, "ReconstructedChargedParticles.tracks_end");
    TTreeReaderArray<int>           reco_tracks_idx  (reader, "_ReconstructedChargedParticles_tracks.index");
    TTreeReaderArray<int>           proj_track_idx   (reader, "_CalorimeterTrackProjections_track.index");
    TTreeReaderArray<int>           nhcal_match_proj (reader, "_HcalEndcapNTrackClusterMatches_track.index");

    while (reader.Next()) {

        // (a) split reco-particle indices by charge sign
        std::vector<int> positive_reco_indices, negative_reco_indices;
        for (size_t reco_index = 0; reco_index < reco_charge.GetSize(); ++reco_index) {
            if (reco_charge[reco_index] > 0) positive_reco_indices.push_back((int) reco_index);
            if (reco_charge[reco_index] < 0) negative_reco_indices.push_back((int) reco_index);
        }

        // (b) try every (+,-) pair and apply the same cut
        for (int positive_reco_index : positive_reco_indices) {
            bool positive_track_in_nhcal = is_track_in_match_list(
                positive_reco_index, reco_tracks_begin, reco_tracks_end,
                reco_tracks_idx, proj_track_idx, nhcal_match_proj);
            for (int negative_reco_index : negative_reco_indices) {
                ++n_bg_total;
                if (positive_track_in_nhcal &&
                    is_track_in_match_list(negative_reco_index, reco_tracks_begin, reco_tracks_end,
                                           reco_tracks_idx, proj_track_idx, nhcal_match_proj))
                    ++n_bg_kept;
            }
        }
    }
    std::cout << "bg:     kept " << n_bg_kept << "/" << n_bg_total << "\n";
    bg_file->Close();
}


// ---------------------------------------------------------------------------
// MAIN MACRO
// ---------------------------------------------------------------------------
void signal_vs_background_cut_efficiency(const char* signal_path = "",
                                         const char* bg_path     = "") {

    if (signal_path[0] == '\0' || bg_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'signal_vs_background_cut_efficiency.C(\"SIGNAL.root\",\"BG.root\")'\n";
        return;
    }

    // ---------- 1. count signal kept + total -----------------------------
    long n_signal_total = 0, n_signal_kept = 0;
    process_signal(signal_path, n_signal_total, n_signal_kept);

    // ---------- 2. count bg kept + total ---------------------------------
    long n_bg_total = 0, n_bg_kept = 0;
    process_background(bg_path, n_bg_total, n_bg_kept);

    // ---------- 3. report efficiencies and ratio -------------------------
    double eff_signal = n_signal_total > 0 ? (double) n_signal_kept / n_signal_total : 0.0;
    double eff_bg     = n_bg_total     > 0 ? (double) n_bg_kept     / n_bg_total     : 0.0;
    double ratio      = eff_bg > 0 ? eff_signal / eff_bg : 0.0;

    std::cout << "signal kept     : " << 100*eff_signal << "%\n";
    std::cout << "bg kept         : " << 100*eff_bg     << "%\n";
    std::cout << "eff(S)/eff(BG)  : " << ratio << "  (>1 means the cut helps)\n";

    // ---------- 4. draw + save 2-bar comparison --------------------------
    TH1F* h = new TH1F("h", Form("Cut efficiency:  S/BG ratio = %.2f;;%% kept", ratio),
                       2, 0, 2);
    h->SetBinContent(1, 100*eff_signal);
    h->SetBinContent(2, 100*eff_bg);
    h->GetXaxis()->SetBinLabel(1, "signal (both Ks)");
    h->GetXaxis()->SetBinLabel(2, "bg (both tracks)");
    h->SetFillColor(kAzure - 4);
    h->SetMinimum(0);

    TCanvas* canvas = new TCanvas("canvas", "", 700, 500);
    h->Draw("BAR");
    TLatex label; label.SetTextAlign(22); label.SetTextSize(0.04);
    label.DrawLatex(h->GetBinCenter(1), 100*eff_signal + 2, Form("%.1f%%", 100*eff_signal));
    label.DrawLatex(h->GetBinCenter(2), 100*eff_bg     + 2, Form("%.1f%%", 100*eff_bg));
    canvas->SaveAs("signal_vs_background_cut_efficiency.pdf");
    std::cout << "wrote signal_vs_background_cut_efficiency.pdf\n";
}
