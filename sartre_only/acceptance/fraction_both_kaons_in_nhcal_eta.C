// fraction_both_kaons_in_nhcal_eta.C
// ==================================
// Uses TTreeReader + TTreeReaderArray.
//
// Run with:
//   root -l -b -q 'fraction_both_kaons_in_nhcal_eta.C("FILE.root")'

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


// ---------------------------------------------------------------------------
// MAIN MACRO
// ---------------------------------------------------------------------------
void fraction_both_kaons_in_nhcal_eta(const char* file_path = "") {

    if (file_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'fraction_both_kaons_in_nhcal_eta.C(\"FILE.root\")'\n";
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

    std::cout << "opened " << file_path << " (" << reader.GetEntries(true) << " events)\n";

    // ---------- 3. counters ----------------------------------------------
    long n_phi_total       = 0;
    long n_phi_any_in_eta  = 0;
    long n_phi_both_in_eta = 0;

    // ---------- 4. main event loop ---------------------------------------
    while (reader.Next()) {

        // (a) get the K+ and K- MC indices for this event's phi -> KK
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;
        ++n_phi_total;

        // (b) ask whether each kaon points at the nHcal acceptance
        double kaon_plus_eta  = pseudorapidity(mc_px[kaon_plus_mc_index],
                                               mc_py[kaon_plus_mc_index],
                                               mc_pz[kaon_plus_mc_index]);
        double kaon_minus_eta = pseudorapidity(mc_px[kaon_minus_mc_index],
                                               mc_py[kaon_minus_mc_index],
                                               mc_pz[kaon_minus_mc_index]);
        bool kaon_plus_in_eta  = (NHCAL_ETA_MIN <= kaon_plus_eta)  && (kaon_plus_eta  < NHCAL_ETA_MAX);
        bool kaon_minus_in_eta = (NHCAL_ETA_MIN <= kaon_minus_eta) && (kaon_minus_eta < NHCAL_ETA_MAX);

        // (c) bump the right counters
        if (kaon_plus_in_eta || kaon_minus_in_eta) ++n_phi_any_in_eta;
        if (kaon_plus_in_eta && kaon_minus_in_eta) ++n_phi_both_in_eta;
    }

    // ---------- 5. print + plot ------------------------------------------
    double p_any  = n_phi_total > 0 ? (double) n_phi_any_in_eta  / n_phi_total : 0.0;
    double p_both = n_phi_total > 0 ? (double) n_phi_both_in_eta / n_phi_total : 0.0;
    std::cout << "phi -> KK decays total          : " << n_phi_total << "\n";
    std::cout << "  >= 1 kaon in nHcal eta range  : " << n_phi_any_in_eta
              << "   (" << 100*p_any  << "%)\n";
    std::cout << "  both kaons in nHcal eta range : " << n_phi_both_in_eta
              << "   (" << 100*p_both << "%)\n";

    TH1F* h = new TH1F("h",
        Form("Truth-eta acceptance in nHcal  (N_phi=%ld);;%% of phi -> KK", n_phi_total),
        2, 0, 2);
    h->SetBinContent(1, 100*p_any);
    h->SetBinContent(2, 100*p_both);
    h->GetXaxis()->SetBinLabel(1, ">= 1 K in eta");
    h->GetXaxis()->SetBinLabel(2, "both Ks in eta");
    h->SetFillColor(kAzure - 4);
    h->SetMinimum(0);

    TCanvas* canvas = new TCanvas("canvas", "", 700, 500);
    h->Draw("BAR");
    TLatex label; label.SetTextAlign(22); label.SetTextSize(0.04);
    label.DrawLatex(h->GetBinCenter(1), 100*p_any  + 2, Form("%.1f%%", 100*p_any));
    label.DrawLatex(h->GetBinCenter(2), 100*p_both + 2, Form("%.1f%%", 100*p_both));
    canvas->SaveAs("fraction_both_kaons_in_nhcal_eta.pdf");
    std::cout << "wrote fraction_both_kaons_in_nhcal_eta.pdf\n";

    root_file->Close();
}
