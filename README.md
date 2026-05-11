#  φ → K⁺K⁻ in nHcal

- **`sartre_only/`** — signal-only studies on Sartre coherent φ → K⁺K⁻.
- **`sartre_and_background/`** — adds the Pythia8 DIS combinatorial bg.
- **`nhcal_vs_ecal/`** — compares the backward HCal vs backward ECal
  cluster matching as a tag for kaons.

Inside `sartre_only/` and `sartre_and_background/` there are three
subfolders by study type:

- `invariant_mass/` — m(K⁺K⁻) histograms.
- `efficiency/` — fraction of kaons (or pairs) surviving a cut.
- `acceptance/` — geometric distributions (truth η, etc.).


Every analysis here walks the same chain of indices in the ROOT file:

```text
MC kaon (truth)
  │  via ReconstructedChargedParticleAssociations
  ▼
reconstructed charged particle  (a track + momentum)
  │  via _ReconstructedChargedParticles_tracks.index
  ▼
parent track  (an entry in CentralCKFTracks)
  │  via _CalorimeterTrackProjections_track.index
  ▼
track projection at the calorimeter face
  │  via _HcalEndcapNTrackClusterMatches_track.index
  ▼
nHcal cluster  (the kaon was "seen in nHcal")
```


## Data
Stream directly with XRootD (no download needed) or copy:

```text
root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.10.2/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_phi_ab_eAu_1.0000.eicrecon.edm4eic.root
```

Background:

```text
root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/26.03.0/epic_craterlake/DIS/NC/18x275/minQ2=1/pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0000.eicrecon.edm4eic.root
```

## How to run

```bash
cd sartre_only/invariant_mass/
root -l -b -q 'truth_kaon_pair_mass.C("root://dtn-eic.jlab.org//volatile/eic/.../sartre_bnonsat_Au_phi_ab_eAu_1.0000.eicrecon.edm4eic.root")'
```

