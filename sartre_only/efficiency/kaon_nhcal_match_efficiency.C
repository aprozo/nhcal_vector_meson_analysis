// kaon_nhcal_match_efficiency.C
// =============================
//
// For every truth kaon from a phi -> K+ K- decay, walk three sequential
// conditions:
//
//     A.  truth-eta in nHcal acceptance              "in eta"
//     B.  reconstructed as a charged particle         "reco"
//     C.  parent track matched to a cluster in
//         HcalEndcapN  (nHcal)  OR  EcalEndcapN       "match"
//
// A minimum-pT cut MIN_PT_GEV is applied on the kaon AFTER the eta cut,
// so the denominator is "kaons in nHcal eta AND pT > MIN_PT_GEV".
//
// Produces two PDFs:
//   kaon_nhcal_match_efficiency.pdf         - integrated chain (4 bars)
//   kaon_nhcal_match_efficiency_vs_pt.pdf   - three efficiencies vs kaon p_T:
//        tracking eff    = N(in-eta AND reco)            / N(in-eta)
//        nHcal match eff = N(in-eta AND reco AND nHcal)  / N(in-eta)
//        ECalN match eff = N(in-eta AND reco AND ECalN)  / N(in-eta)


#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

const int    PDG_PHI       = 333;
const int    PDG_KP        = 321;
const double NHCAL_ETA_MIN = -4.05;
const double NHCAL_ETA_MAX = -1.20;
const double MIN_PT_GEV    = 0.2;   // skip kaons softer than this (curl in B-field)


// pseudorapidity from 3-momentum (uses |p|, not E - this is eta, not y)
double pseudorapidity(double px, double py, double pz) {
    double momentum_magnitude = std::sqrt(px*px + py*py + pz*pz);
    if (momentum_magnitude == 0 || std::fabs(pz) >= momentum_magnitude) return std::nan("");
    return 0.5 * std::log((momentum_magnitude + pz) / (momentum_magnitude - pz));
}

// transverse momentum p_T = sqrt(px^2 + py^2)
double transverse_momentum(double px, double py) {
    return std::sqrt(px*px + py*py);
}

// safe division: 0 if denominator is 0
double safe_div(double numerator, double denominator) {
    return denominator > 0 ? numerator / denominator : 0.0;
}

// find the two kaon MC indices for a phi -> K+ K- decay
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

// reco-particle index -> parent track has an entry in `match_proj`?
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
// MAIN MACRO
// ---------------------------------------------------------------------------
void kaon_nhcal_match_efficiency(const char* file_path = "") {

    if (file_path[0] == '\0') {
        std::cerr << "usage: root -l -b -q 'kaon_nhcal_match_efficiency.C(\"FILE.root\")'\n";
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
    std::cout << "min p_T cut on kaons : " << MIN_PT_GEV << " GeV\n";

    // ---------- 3. counters for the integrated chain ---------------------
    long n_kaons_total                     = 0;
    long n_kaons_in_eta_passing_pt         = 0;   // denominator after eta + pT cut
    long n_kaons_reco                      = 0;   //   + reconstructed
    long n_kaons_reco_and_nhcal_matched    = 0;   //   + nHcal match
    long n_kaons_reco_and_ecaln_matched    = 0;   //   + ECalN match

    // ---------- 4. pT histograms for the differential plot ---------------
    // Sartre coherent phi at 18x110 produces soft kaons: the bulk lives
    // in 0 - 2 GeV. Bin finely there.
    const int    n_pt_bins = 15;
    const double pt_min    = 0.0;
    const double pt_max    = 3.0;
    TH1F* h_pt_denom       = new TH1F("h_pt_denom",
        ";kaon p_{T} [GeV];kaons", n_pt_bins, pt_min, pt_max);
    TH1F* h_pt_reco        = new TH1F("h_pt_reco",
        ";kaon p_{T} [GeV];kaons", n_pt_bins, pt_min, pt_max);
    TH1F* h_pt_nhcal_match = new TH1F("h_pt_nhcal_match",
        ";kaon p_{T} [GeV];kaons", n_pt_bins, pt_min, pt_max);
    TH1F* h_pt_ecaln_match = new TH1F("h_pt_ecaln_match",
        ";kaon p_{T} [GeV];kaons", n_pt_bins, pt_min, pt_max);

    // ---------- 5. main event loop ---------------------------------------
    while (reader.Next()) {

        // (a) find the two truth kaons from the phi -> KK decay
        int kaon_plus_mc_index, kaon_minus_mc_index;
        if (!find_phi_to_kaons(mc_pdg, mc_gen_status,
                               mc_daughters_begin, mc_daughters_end, mc_daughter_index,
                               kaon_plus_mc_index, kaon_minus_mc_index))
            continue;
        int kaon_mc_indices[2] = {kaon_plus_mc_index, kaon_minus_mc_index};

        // (b) walk each kaon through eta + pT + reco + match
        for (int slot = 0; slot < 2; ++slot) {
            int mc_kaon_index = kaon_mc_indices[slot];
            ++n_kaons_total;

            // (b.i) truth eta in nHcal acceptance
            double kaon_truth_eta = pseudorapidity(mc_px[mc_kaon_index],
                                                   mc_py[mc_kaon_index],
                                                   mc_pz[mc_kaon_index]);
            bool kaon_in_eta = (NHCAL_ETA_MIN <= kaon_truth_eta) && (kaon_truth_eta < NHCAL_ETA_MAX);
            if (!kaon_in_eta) continue;

            // (b.ii) minimum-pT cut
            double kaon_truth_pt = transverse_momentum(mc_px[mc_kaon_index],
                                                      mc_py[mc_kaon_index]);
            if (kaon_truth_pt < MIN_PT_GEV) continue;
            ++n_kaons_in_eta_passing_pt;
            h_pt_denom->Fill(kaon_truth_pt);

            // (b.iii) reconstructed as a charged particle
            int kaon_reco_index = reco_index_of(mc_kaon_index, assoc_sim_id, assoc_rec_id);
            if (kaon_reco_index < 0) continue;
            ++n_kaons_reco;
            h_pt_reco->Fill(kaon_truth_pt);

            // (b.iv) reco track has a nHcal cluster match
            if (is_track_in_match_list(kaon_reco_index,
                                       reco_tracks_begin, reco_tracks_end, reco_tracks_idx,
                                       proj_track_idx, nhcal_match_proj)) {
                ++n_kaons_reco_and_nhcal_matched;
                h_pt_nhcal_match->Fill(kaon_truth_pt);
            }
            // (b.v) reco track has an ECalN cluster match (parallel to b.iv)
            if (is_track_in_match_list(kaon_reco_index,
                                       reco_tracks_begin, reco_tracks_end, reco_tracks_idx,
                                       proj_track_idx, ecaln_match_proj)) {
                ++n_kaons_reco_and_ecaln_matched;
                h_pt_ecaln_match->Fill(kaon_truth_pt);
            }
        }
    }

    // ---------- 6. print integrated efficiencies -------------------------
    double p_reco_given_denom        = safe_div(n_kaons_reco,                   n_kaons_in_eta_passing_pt);
    double p_nhcal_match_given_denom = safe_div(n_kaons_reco_and_nhcal_matched, n_kaons_in_eta_passing_pt);
    double p_ecaln_match_given_denom = safe_div(n_kaons_reco_and_ecaln_matched, n_kaons_in_eta_passing_pt);

    std::cout << "kaons total                                : " << n_kaons_total << "\n";
    std::cout << "  in nHcal eta AND p_T > " << MIN_PT_GEV << " GeV "
              << " : " << n_kaons_in_eta_passing_pt << "  (denominator)\n";
    std::cout << "    AND reconstructed                      : " << n_kaons_reco
              << "  (" << 100*p_reco_given_denom << "%)\n";
    std::cout << "      AND nHcal cluster-matched            : " << n_kaons_reco_and_nhcal_matched
              << "  (" << 100*p_nhcal_match_given_denom << "%)\n";
    std::cout << "      AND ECalN cluster-matched            : " << n_kaons_reco_and_ecaln_matched
              << "  (" << 100*p_ecaln_match_given_denom << "%)\n";

    // ---------- 7. plot A: integrated 4-bar chain ------------------------
    TH1F* h_integrated = new TH1F("h_integrated",
        Form("Per-kaon nHcal vs ECalN match chain (p_{T} > %.2f GeV);;"
             "%% of kaons in nHcal #eta with p_{T} > %.2f GeV",
             MIN_PT_GEV, MIN_PT_GEV),
        4, 0, 4);
    h_integrated->SetBinContent(1, 100.0);                                // denominator = 100% by definition
    h_integrated->SetBinContent(2, 100 * p_reco_given_denom);
    h_integrated->SetBinContent(3, 100 * p_nhcal_match_given_denom);
    h_integrated->SetBinContent(4, 100 * p_ecaln_match_given_denom);
    h_integrated->GetXaxis()->SetBinLabel(1, "in #eta AND p_{T} cut");
    h_integrated->GetXaxis()->SetBinLabel(2, "+ reco");
    h_integrated->GetXaxis()->SetBinLabel(3, "+ nHcal match");
    h_integrated->GetXaxis()->SetBinLabel(4, "+ ECalN match");
    h_integrated->SetFillColor(kAzure - 4);
    h_integrated->SetMinimum(0);
    h_integrated->SetStats(0);

    TCanvas* canvas_integrated = new TCanvas("canvas_integrated", "", 800, 600);
    h_integrated->Draw("BAR");
    TLatex label_integrated; label_integrated.SetTextAlign(22); label_integrated.SetTextSize(0.04);
    for (int i = 1; i <= 4; ++i)
        label_integrated.DrawLatex(h_integrated->GetBinCenter(i), h_integrated->GetBinContent(i) + 2,
                                   Form("%.1f%%", h_integrated->GetBinContent(i)));
    canvas_integrated->SaveAs("kaon_nhcal_match_efficiency.pdf");
    std::cout << "wrote kaon_nhcal_match_efficiency.pdf\n";

    // ---------- 8. plot B: efficiencies vs kaon p_T ----------------------
    // TH1::Divide(num, denom, c1, c2, "b") gives binomial errors.
    TH1F* h_eff_reco        = (TH1F*) h_pt_reco       ->Clone("h_eff_reco");
    h_eff_reco       ->Divide(h_pt_reco,        h_pt_denom, 1.0, 1.0, "b");
    h_eff_reco       ->SetTitle(Form("nHcal vs ECalN match chain vs kaon p_{T} "
                                     "(denom = kaons in nHcal #eta, p_{T} > %.2f GeV);"
                                     "kaon p_{T} [GeV];efficiency", MIN_PT_GEV));

    TH1F* h_eff_nhcal_match = (TH1F*) h_pt_nhcal_match->Clone("h_eff_nhcal_match");
    h_eff_nhcal_match->Divide(h_pt_nhcal_match, h_pt_denom, 1.0, 1.0, "b");

    TH1F* h_eff_ecaln_match = (TH1F*) h_pt_ecaln_match->Clone("h_eff_ecaln_match");
    h_eff_ecaln_match->Divide(h_pt_ecaln_match, h_pt_denom, 1.0, 1.0, "b");

    h_eff_reco       ->SetMarkerStyle(20); h_eff_reco       ->SetMarkerColor(kBlack);
    h_eff_reco       ->SetLineColor(kBlack);
    h_eff_nhcal_match->SetMarkerStyle(21); h_eff_nhcal_match->SetMarkerColor(kRed + 1);
    h_eff_nhcal_match->SetLineColor(kRed + 1);
    h_eff_ecaln_match->SetMarkerStyle(22); h_eff_ecaln_match->SetMarkerColor(kGreen + 2);
    h_eff_ecaln_match->SetLineColor(kGreen + 2);
    h_eff_reco       ->GetYaxis()->SetRangeUser(0.0, 1.05);
    h_eff_reco       ->SetStats(0);
    h_eff_nhcal_match->SetStats(0);
    h_eff_ecaln_match->SetStats(0);

    TCanvas* canvas_pt = new TCanvas("canvas_pt", "", 800, 600);
    canvas_pt->SetGridy(1);
    h_eff_reco       ->Draw("E1");
    h_eff_nhcal_match->Draw("E1 SAME");
    h_eff_ecaln_match->Draw("E1 SAME");

    TLegend* legend_pt = new TLegend(0.45, 0.70, 0.88, 0.88);
    legend_pt->AddEntry(h_eff_reco,        "tracking eff:  reco | in #eta",         "lep");
    legend_pt->AddEntry(h_eff_nhcal_match, "nHcal match eff:  nHcal | in #eta",     "lep");
    legend_pt->AddEntry(h_eff_ecaln_match, "ECalN match eff:  ECalN | in #eta",     "lep");
    legend_pt->Draw();

    canvas_pt->SaveAs("kaon_nhcal_match_efficiency_vs_pt.pdf");
    std::cout << "wrote kaon_nhcal_match_efficiency_vs_pt.pdf\n";

    root_file->Close();
}
