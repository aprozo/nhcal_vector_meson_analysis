// reco_mass_with_both_kaons_in_nhcal.C
// ====================================
// Headline plot: m(K+K-) for events where both reconstructed kaons
// were matched to a cluster in the backward HCal.
//
// Run with:
//   root -l -b -q 'reco_mass_with_both_kaons_in_nhcal.C("FILE.root")'

#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>

const int    PDG_PHI       = 333;
const int    PDG_KP        = 321;
const double KAON_MASS_GEV = 0.493677;


// ---------------------------------------------------------------------------
// find_phi_to_kaons:
//   return the K+ and K- MC indices from the first phi -> KK decay.
// ---------------------------------------------------------------------------
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
        int kaon_plus_local = -1, kaon_minus_local = -1;
        for (unsigned int daughter_slot = mc_daughters_begin[mc_particle_index];
                          daughter_slot < mc_daughters_end[mc_particle_index];
                          ++daughter_slot) {
            int child_mc_index = mc_daughter_index[daughter_slot];
            if (mc_pdg[child_mc_index] ==  PDG_KP) kaon_plus_local  = child_mc_index;
            if (mc_pdg[child_mc_index] == -PDG_KP) kaon_minus_local = child_mc_index;
        }
        if (kaon_plus_local >= 0 && kaon_minus_local >= 0) {
            kaon_plus_mc_index  = kaon_plus_local;
            kaon_minus_mc_index = kaon_minus_local;
            return true;
        }
    }
    return false;
}


// ---------------------------------------------------------------------------
// reco_index_of:
//   look up an MC particle in the MC<->reco association table and
//   return the index of the matched reconstructed charged particle
//   (or -1 if there is no match).
// ---------------------------------------------------------------------------
int reco_index_of(int mc_index,
                  const TTreeReaderArray<int>&          assoc_sim_id,
                  const TTreeReaderArray<int>&          assoc_rec_id)
{
    for (size_t assoc_slot = 0; assoc_slot < assoc_sim_id.GetSize(); ++assoc_slot)
        if (assoc_sim_id[assoc_slot] == mc_index) return assoc_rec_id[assoc_slot];
    return -1;
}


// ---------------------------------------------------------------------------
// is_kaon_matched_in_nhcal:
//   given a reco charged particle index, follow it to its parent track
//   and ask whether that track appears in HcalEndcapNTrackClusterMatches.
// ---------------------------------------------------------------------------
bool is_kaon_matched_in_nhcal(int reco_index,
                              const TTreeReaderArray<unsigned int>& reco_tracks_begin,
                              const TTreeReaderArray<unsigned int>& reco_tracks_end,
                              const TTreeReaderArray<int>&          reco_tracks_idx,
                              const TTreeReaderArray<int>&          proj_track_idx,
                              const TTreeReaderArray<int>&          nhcal_match_proj)
{
    if (reco_index < 0) return false;
    unsigned int tracks_begin_for_reco = reco_tracks_begin[reco_index];
    unsigned int tracks_end_for_reco   = reco_tracks_end[reco_index];
    if (tracks_end_for_reco <= tracks_begin_for_reco) return false;
    int track_index = reco_tracks_idx[tracks_begin_for_reco];
    for (size_t match_slot = 0; match_slot < nhcal_match_proj.GetSize(); ++match_slot) {
        int projection_index = nhcal_match_proj[match_slot];
        if (projection_index < 0 || (size_t) projection_index >= proj_track_idx.GetSize()) continue;
        if (proj_track_idx[projection_index] == track_index) return true;
    }
    return false;
}


// ---------------------------------------------------------------------------
// MAIN MACRO
// ---------------------------------------------------------------------------
void reco_mass_with_both_kaons_in_nhcal(const char* file_path = "") {

    if (file_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'reco_mass_with_both_kaons_in_nhcal.C(\"FILE.root\")'\n";
        return;
    }

    // ---------- 1. open the file -----------------------------------------
    TFile* root_file = TFile::Open(file_path);
    if (!root_file || root_file->IsZombie()) { std::cerr << "cannot open\n"; return; }

    // ---------- 2. bind every needed branch ------------------------------
    TTreeReader reader("events", root_file);
    TTreeReaderArray<int>           mc_pdg             (reader, "MCParticles.PDG");
    TTreeReaderArray<int>           mc_gen_status      (reader, "MCParticles.generatorStatus");
    TTreeReaderArray<unsigned int>  mc_daughters_begin (reader, "MCParticles.daughters_begin");
    TTreeReaderArray<unsigned int>  mc_daughters_end   (reader, "MCParticles.daughters_end");
    TTreeReaderArray<int>           mc_daughter_index  (reader, "_MCParticles_daughters.index");
    TTreeReaderArray<int>           assoc_sim_id       (reader, "_ReconstructedChargedParticleAssociations_sim.index");
    TTreeReaderArray<int>           assoc_rec_id       (reader, "_ReconstructedChargedParticleAssociations_rec.index");
    TTreeReaderArray<float>         reco_px            (reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float>         reco_py            (reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float>         reco_pz            (reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<unsigned int>  reco_tracks_begin  (reader, "ReconstructedChargedParticles.tracks_begin");
    TTreeReaderArray<unsigned int>  reco_tracks_end    (reader, "ReconstructedChargedParticles.tracks_end");
    TTreeReaderArray<int>           reco_tracks_idx    (reader, "_ReconstructedChargedParticles_tracks.index");
    TTreeReaderArray<int>           proj_track_idx     (reader, "_CalorimeterTrackProjections_track.index");
    TTreeReaderArray<int>           nhcal_match_proj   (reader, "_HcalEndcapNTrackClusterMatches_track.index");

    std::cout << "opened " << file_path << " (" << reader.GetEntries(true) << " events)\n";

    // ---------- 3. output histogram --------------------------------------
    TH1F* h_mKK = new TH1F("h_mKK",
        "Reco m(K+K-), both kaons nHcal-matched;m(K+K-) [GeV];entries / bin",
        120, 0.98, 1.10);

    // ---------- 4. main event loop ---------------------------------------
    int n_filled = 0;
    while (reader.Next()) {

        // (a) find the truth K+ and K- from a phi decay
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;

        // (b) walk MC -> reco association to find the reconstructed kaons
        int kaon_plus_reco_index  = reco_index_of(kaon_plus_mc_index,  assoc_sim_id, assoc_rec_id);
        int kaon_minus_reco_index = reco_index_of(kaon_minus_mc_index, assoc_sim_id, assoc_rec_id);
        if (kaon_plus_reco_index < 0 || kaon_minus_reco_index < 0) continue;

        // (c) require BOTH reco kaons to be matched to a nHcal cluster
        bool kaon_plus_in_nhcal  = is_kaon_matched_in_nhcal(kaon_plus_reco_index,
            reco_tracks_begin, reco_tracks_end, reco_tracks_idx,
            proj_track_idx, nhcal_match_proj);
        bool kaon_minus_in_nhcal = is_kaon_matched_in_nhcal(kaon_minus_reco_index,
            reco_tracks_begin, reco_tracks_end, reco_tracks_idx,
            proj_track_idx, nhcal_match_proj);
        if (!(kaon_plus_in_nhcal && kaon_minus_in_nhcal)) continue;

        // (d) build reco four-vectors with the kaon mass hypothesis
        double kaon_plus_px = reco_px[kaon_plus_reco_index];
        double kaon_plus_py = reco_py[kaon_plus_reco_index];
        double kaon_plus_pz = reco_pz[kaon_plus_reco_index];
        double kaon_plus_energy = std::sqrt(kaon_plus_px*kaon_plus_px
                                          + kaon_plus_py*kaon_plus_py
                                          + kaon_plus_pz*kaon_plus_pz
                                          + KAON_MASS_GEV*KAON_MASS_GEV);

        double kaon_minus_px = reco_px[kaon_minus_reco_index];
        double kaon_minus_py = reco_py[kaon_minus_reco_index];
        double kaon_minus_pz = reco_pz[kaon_minus_reco_index];
        double kaon_minus_energy = std::sqrt(kaon_minus_px*kaon_minus_px
                                           + kaon_minus_py*kaon_minus_py
                                           + kaon_minus_pz*kaon_minus_pz
                                           + KAON_MASS_GEV*KAON_MASS_GEV);

        // (e) m(K+K-) and fill the histogram
        double total_energy = kaon_plus_energy + kaon_minus_energy;
        double total_px     = kaon_plus_px     + kaon_minus_px;
        double total_py     = kaon_plus_py     + kaon_minus_py;
        double total_pz     = kaon_plus_pz     + kaon_minus_pz;
        double m_squared = total_energy*total_energy
                         - (total_px*total_px + total_py*total_py + total_pz*total_pz);
        double m_kaon_pair = std::sqrt(std::max(m_squared, 0.0));
        h_mKK->Fill(m_kaon_pair);
        ++n_filled;
    }
    std::cout << "phi -> KK with BOTH kaons nHcal-matched: " << n_filled << "\n";

    // ---------- 5. draw + save the plot ----------------------------------
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    h_mKK->SetLineWidth(2);
    h_mKK->Draw("HIST");
    TLine pdg_line(1.020, 0, 1.020, h_mKK->GetMaximum());
    pdg_line.SetLineStyle(2); pdg_line.SetLineColor(kGray + 1);
    pdg_line.Draw();
    canvas->SaveAs("reco_mass_with_both_kaons_in_nhcal.pdf");
    std::cout << "wrote reco_mass_with_both_kaons_in_nhcal.pdf\n";

    root_file->Close();
}
