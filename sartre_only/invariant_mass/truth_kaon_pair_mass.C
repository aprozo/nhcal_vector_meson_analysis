// truth_kaon_pair_mass.C
// ======================
// Uses TTreeReader + TTreeReaderArray, the modern clean ROOT pattern.
//
// Run with:
//   root -l -b -q 'truth_kaon_pair_mass.C("FILE.root")'

#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>

// physics constants used below
const int PDG_PHI = 333;
const int PDG_KP  = 321;


// ---------------------------------------------------------------------------
// find_phi_to_kaons:
//   walk the MC particle list, return through reference parameters the
//   indices of K+ and K- coming from the FIRST phi(1020) -> K+K- decay.
//   Returns true if both kaons were found.
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

        // is this an MC phi(1020) that decayed (status 1 = stable, 2 = decayed)?
        if (mc_pdg[mc_particle_index] != PDG_PHI) continue;
        int gen_status_value = mc_gen_status[mc_particle_index];
        if (gen_status_value != 1 && gen_status_value != 2) continue;

        // scan the phi's daughter list for a K+ and a K-
        int kaon_plus_local  = -1;
        int kaon_minus_local = -1;
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
// invariant_mass_two_body:
//   sqrt( (E1+E2)^2 - |p1+p2|^2 ) for a two-particle system.
// ---------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------
// MAIN MACRO
// ---------------------------------------------------------------------------
void truth_kaon_pair_mass(const char* file_path = "") {

    if (file_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'truth_kaon_pair_mass.C(\"FILE.root\")'\n";
        return;
    }

    // ---------- 1. open the file -----------------------------------------
    TFile* root_file = TFile::Open(file_path);
    if (!root_file || root_file->IsZombie()) {
        std::cerr << "could not open " << file_path << "\n";
        return;
    }

    // ---------- 2. bind every needed branch to a typed array -------------
    TTreeReader reader("events", root_file);
    TTreeReaderArray<int>           mc_pdg            (reader, "MCParticles.PDG");
    TTreeReaderArray<int>           mc_gen_status     (reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double>        mc_px             (reader, "MCParticles.momentum.x");
    TTreeReaderArray<double>        mc_py             (reader, "MCParticles.momentum.y");
    TTreeReaderArray<double>        mc_pz             (reader, "MCParticles.momentum.z");
    TTreeReaderArray<double>        mc_mass           (reader, "MCParticles.mass");
    TTreeReaderArray<unsigned int>  mc_daughters_begin(reader, "MCParticles.daughters_begin");
    TTreeReaderArray<unsigned int>  mc_daughters_end  (reader, "MCParticles.daughters_end");
    TTreeReaderArray<int>           mc_daughter_index (reader, "_MCParticles_daughters.index");

    std::cout << "opened " << file_path << " (" << reader.GetEntries(true) << " events)\n";

    // ---------- 3. output histogram --------------------------------------
    TH1F* h_mKK = new TH1F("h_mKK",
        "Truth invariant mass, phi -> K+K-;m(K+K-) [GeV];entries / bin",
        120, 0.98, 1.10);

    // ---------- 4. main event loop ---------------------------------------
    int n_found = 0;
    while (reader.Next()) {

        // (a) find the K+ K- pair from a phi(1020) -> KK decay
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;

        // (b) build the K+ truth four-vector (energy from mass + p)
        double kaon_plus_px = mc_px[kaon_plus_mc_index];
        double kaon_plus_py = mc_py[kaon_plus_mc_index];
        double kaon_plus_pz = mc_pz[kaon_plus_mc_index];
        double kaon_plus_mass = mc_mass[kaon_plus_mc_index];
        double kaon_plus_energy = std::sqrt(kaon_plus_px*kaon_plus_px
                                          + kaon_plus_py*kaon_plus_py
                                          + kaon_plus_pz*kaon_plus_pz
                                          + kaon_plus_mass*kaon_plus_mass);

        // (c) build the K- truth four-vector
        double kaon_minus_px = mc_px[kaon_minus_mc_index];
        double kaon_minus_py = mc_py[kaon_minus_mc_index];
        double kaon_minus_pz = mc_pz[kaon_minus_mc_index];
        double kaon_minus_mass = mc_mass[kaon_minus_mc_index];
        double kaon_minus_energy = std::sqrt(kaon_minus_px*kaon_minus_px
                                           + kaon_minus_py*kaon_minus_py
                                           + kaon_minus_pz*kaon_minus_pz
                                           + kaon_minus_mass*kaon_minus_mass);

        // (d) m(K+K-) and fill the histogram
        double m_kaon_pair = invariant_mass_two_body(
            kaon_plus_energy,  kaon_plus_px,  kaon_plus_py,  kaon_plus_pz,
            kaon_minus_energy, kaon_minus_px, kaon_minus_py, kaon_minus_pz);
        h_mKK->Fill(m_kaon_pair);
        ++n_found;
    }
    std::cout << "found " << n_found << " phi -> K+ K- decays\n";

    // ---------- 5. draw + save the plot ----------------------------------
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    h_mKK->SetLineWidth(2);
    h_mKK->Draw("HIST");
    TLine pdg_line(1.020, 0, 1.020, h_mKK->GetMaximum());
    pdg_line.SetLineStyle(2); pdg_line.SetLineColor(kGray + 1);
    pdg_line.Draw();
    canvas->SaveAs("truth_kaon_pair_mass.pdf");
    std::cout << "wrote truth_kaon_pair_mass.pdf\n";

    root_file->Close();
}
