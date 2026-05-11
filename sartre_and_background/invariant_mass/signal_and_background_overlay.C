// signal_and_background_overlay.C
// ===============================
// Uses TTreeReader + TTreeReaderArray; named helper functions.
//
// Run with:
//   root -l -b -q 'signal_and_background_overlay.C("SIGNAL.root","BG.root")'

#include <iostream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

const int    PDG_PHI       = 333;
const int    PDG_KP        = 321;
const double KAON_MASS_GEV = 0.493677;


// invariant mass of a 2-body system from individual four-momenta
double invariant_mass_two_body(double e1, double px1, double py1, double pz1,
                               double e2, double px2, double py2, double pz2)
{
    double total_e  = e1  + e2;
    double total_px = px1 + px2;
    double total_py = py1 + py2;
    double total_pz = pz1 + pz2;
    double m_squared = total_e*total_e
                     - (total_px*total_px + total_py*total_py + total_pz*total_pz);
    return std::sqrt(std::max(m_squared, 0.0));
}

// energy of a particle under the kaon mass hypothesis
double energy_with_kaon_mass(double px, double py, double pz) {
    return std::sqrt(px*px + py*py + pz*pz + KAON_MASS_GEV*KAON_MASS_GEV);
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


// ---------------------------------------------------------------------------
// SIGNAL processing: fill h_signal with reco m(K+K-) for every truth phi -> KK
// ---------------------------------------------------------------------------
void process_signal(const char* signal_path, TH1F* h_signal) {

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
    TTreeReaderArray<float>         reco_px           (reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float>         reco_py           (reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float>         reco_pz           (reader, "ReconstructedChargedParticles.momentum.z");

    int n_filled = 0;
    while (reader.Next()) {

        // (a) find the truth K+ and K- from a phi decay
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;

        // (b) MC -> reco lookup for both kaons
        int kaon_plus_reco_index  = reco_index_of(kaon_plus_mc_index,  assoc_sim_id, assoc_rec_id);
        int kaon_minus_reco_index = reco_index_of(kaon_minus_mc_index, assoc_sim_id, assoc_rec_id);
        if (kaon_plus_reco_index < 0 || kaon_minus_reco_index < 0) continue;

        // (c) fill m(K+K-) from reco four-vectors (kaon mass hypothesis)
        double kaon_plus_energy  = energy_with_kaon_mass(reco_px[kaon_plus_reco_index],
                                                        reco_py[kaon_plus_reco_index],
                                                        reco_pz[kaon_plus_reco_index]);
        double kaon_minus_energy = energy_with_kaon_mass(reco_px[kaon_minus_reco_index],
                                                        reco_py[kaon_minus_reco_index],
                                                        reco_pz[kaon_minus_reco_index]);
        double m = invariant_mass_two_body(
            kaon_plus_energy,
            reco_px[kaon_plus_reco_index],  reco_py[kaon_plus_reco_index],  reco_pz[kaon_plus_reco_index],
            kaon_minus_energy,
            reco_px[kaon_minus_reco_index], reco_py[kaon_minus_reco_index], reco_pz[kaon_minus_reco_index]);
        h_signal->Fill(m);
        ++n_filled;
    }
    std::cout << "signal phi -> KK reco: " << n_filled << "\n";
    signal_file->Close();
}


// ---------------------------------------------------------------------------
// BACKGROUND processing: fill h_bg with m(K+K-) for every OS reco pair
// ---------------------------------------------------------------------------
void process_background(const char* bg_path, TH1F* h_bg) {

    std::cout << "bg: " << bg_path << "\n";
    TFile* bg_file = TFile::Open(bg_path);
    TTreeReader reader("events", bg_file);

    TTreeReaderArray<float>         reco_charge(reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<float>         reco_px    (reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float>         reco_py    (reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float>         reco_pz    (reader, "ReconstructedChargedParticles.momentum.z");

    long n_pairs = 0;
    while (reader.Next()) {

        // (a) split reco-particle indices into positive and negative charges
        std::vector<int> positive_reco_indices, negative_reco_indices;
        for (size_t reco_index = 0; reco_index < reco_charge.GetSize(); ++reco_index) {
            if (reco_charge[reco_index] > 0) positive_reco_indices.push_back((int) reco_index);
            if (reco_charge[reco_index] < 0) negative_reco_indices.push_back((int) reco_index);
        }

        // (b) form every (+,-) pair with the kaon mass hypothesis
        for (int positive_reco_index : positive_reco_indices) {
            double positive_energy = energy_with_kaon_mass(reco_px[positive_reco_index],
                                                           reco_py[positive_reco_index],
                                                           reco_pz[positive_reco_index]);
            for (int negative_reco_index : negative_reco_indices) {
                double negative_energy = energy_with_kaon_mass(reco_px[negative_reco_index],
                                                               reco_py[negative_reco_index],
                                                               reco_pz[negative_reco_index]);
                double m = invariant_mass_two_body(
                    positive_energy,
                    reco_px[positive_reco_index], reco_py[positive_reco_index], reco_pz[positive_reco_index],
                    negative_energy,
                    reco_px[negative_reco_index], reco_py[negative_reco_index], reco_pz[negative_reco_index]);
                h_bg->Fill(m);
                ++n_pairs;
            }
        }
    }
    std::cout << "bg opposite-sign pairs: " << n_pairs << "\n";
    bg_file->Close();
}


// ---------------------------------------------------------------------------
// MAIN MACRO
// ---------------------------------------------------------------------------
void signal_and_background_overlay(const char* signal_path = "",
                                   const char* bg_path     = "") {

    if (signal_path[0] == '\0' || bg_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'signal_and_background_overlay.C(\"SIGNAL.root\",\"BG.root\")'\n";
        return;
    }

    // ---------- 1. histograms --------------------------------------------
    TH1F* h_signal = new TH1F("h_signal", "signal: phi -> KK reco;m(K+K-) [GeV];entries / bin",
                              120, 0.98, 1.10);
    TH1F* h_bg     = new TH1F("h_bg",     "bg: OS K-hyp pairs", 120, 0.98, 1.10);

    // ---------- 2. fill both --------------------------------------------
    process_signal    (signal_path, h_signal);
    process_background(bg_path,     h_bg);

    // ---------- 3. draw + save the overlay -------------------------------
    TCanvas* canvas = new TCanvas("canvas", "", 900, 600);
    h_bg->SetFillColorAlpha(kRed, 0.30);
    h_bg->SetLineColor(kRed);
    h_bg->Draw("HIST");
    h_signal->SetLineColor(kBlue + 1);
    h_signal->SetLineWidth(2);
    h_signal->Draw("HIST SAME");
    h_bg->SetTitle("Signal (Sartre) vs combinatorial background (Pythia8 NC DIS);m(K+K-) [GeV];entries / bin");

    TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.88);
    legend->AddEntry(h_signal, "signal: phi -> KK reco", "l");
    legend->AddEntry(h_bg,     "bg: DIS OS K-hyp pairs", "f");
    legend->Draw();
    canvas->SaveAs("signal_and_background_overlay.pdf");
    std::cout << "wrote signal_and_background_overlay.pdf\n";
}
