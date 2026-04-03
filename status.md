# HH→bbττ Scouting Analysis — Current Status

*Last updated: 2026-03-20*

## Overview

Search for non-resonant Higgs pair production in the bbττ final state using
CMS Run 3 (2024) scouting data at √s = 13.6 TeV. The analysis uses
ScoutingPFJetRecluster jets with ParticleNet b-tagging in a fully hadronic
selection (≥4 jets, dijet mass windows for H→bb and H→ττ candidates).

**Luminosity**: 103.965 fb⁻¹ (from brilcalc, Run2024 golden JSON)


## Branch Structure

| Branch | Purpose | Status |
|--------|---------|--------|
| `scouting` | **Active** — main development branch | Current |
| `quickanddirtyselection` | Previous working branch (pre-restructure) | Archived |
| `lite` | Pre-scouting era (ABCD, gen matching, 2025) | Archived |


## MC Samples

### Backgrounds

| Process | Samples | Binning | Generator | Status |
|---------|---------|---------|-----------|--------|
| DY (Z/γ*→ℓℓ) | 5 | pT(ℓℓ) bins (40-100, 100, 200, 400, 600+) | amcatnloFXFX | Active |
| TT (tt̄) | 3 | By W decay (2L2Nu, 4Q, LNu2Q) | powheg | Active |
| QCD multijet | 7 | HT bins (40-70 to 2000+) | madgraphMLM | **New** |
| **Missing** | — | QCD HT 400-600, 600-800 | — | Not produced |

### Signal

| Process | Coupling | xsec [pb] | Status |
|---------|----------|-----------|--------|
| ggF HH→bbττ (SM) | kl=1, c2=0 | 0.002506 | Active |
| 9 BSM benchmarks | Various kl, c2 | Placeholder | Commented out in config |


## Decay Mode Classification (12 MC groups)

All MC events get a single `decayMode` integer from `Ana::GlobalDecayMode()`:

| Mode | Process | Channel | Color |
|------|---------|---------|-------|
| 1-5 | DY | ee, μμ, τhτh, τμτh, τeτh | gold → tomato |
| 10-12 | TT | hadronic, semi-lep, di-lep | green → teal |
| 20-22 | Signal | bbτhτh, bbτμτh, bbτeτh | red, pink, purple |
| 30 | QCD | multijet (no gen matching) | cyan |


## Selection (Cutflow)

```
Trigger → ≥4 jets → 4 jets pT > 20 → 2 b-tags (BvsAll > 0.1)
→ H→bb mass [100,150] GeV → H→ττ mass [40,150] GeV
```

Three trigger paths analyzed in parallel:
- **NoTrigger** — preselection only (no trigger requirement)
- **DST_JetHT** — `DST_PFScouting_JetHT` (scouting hadronic)
- **PARKING_HH** — OR of 11 HLT quad-jet/b-tag paths (parking)


## Plot Types

| Type | Description | Output |
|------|-------------|--------|
| Stacked | MC stacked histogram (12 groups) | `stacked/{trig}_{var}.png` |
| Shape | Normalized MC overlay | `shape/{trig}_{var}.png` |
| Trigger overlay | All-MC summed, 3 triggers overlaid | `overlay/{var}.png` |
| Per-channel overlay | Signal-only per decay channel | `overlay/{ch}_{var}.png` |
| Eff stacked | Stacked + trigger efficiency panel | `eff/{trig}/{var}.png` |
| Sig stacked | Stacked + bin-by-bin S/√B panel | `sig/{trig}/{var}.png` |
| Cum sig | Stacked + cumulative S/√B (right-to-left) | `cuml_sig/{trig}/{var}.png` |


## Code Architecture

```
Ana/
├── cutflow_TrigEff.py      Orchestrator (~590 LOC): CLI, phase sequencing
├── make_pdf.py             PDF comparison tool (side-by-side triggers)
├── genAna.py               Gen-level analysis (signal MC only)
├── analysis/               Core analysis modules (extracted from cutflow_TrigEff.py)
│   ├── config.py           YAML config loaders (variables, triggers, decay modes)
│   ├── definitions.py      RDataFrame .Define() chains (kinematics, gen matching)
│   ├── histograms.py       Histogram booking & materialisation
│   ├── cutflow.py          Cutflow booking, extraction, markdown formatting
│   ├── cache.py            Histogram cache save/load/invalidation
│   └── plots.py            Plot orchestration (Phase 3, 3.5, 5)
├── elements/               C++ engine (compiled via ACLiC)
│   ├── common.h            Shared types, PDG constants, utilities
│   ├── GenMatching.C       Gen-level matching (DY, TT, Signal, GlobalDecayMode)
│   └── RecoObjects.C       Reco selection (leptons, b-jets, dijet pairing)
├── config/                 YAML configuration (no magic numbers in C++)
│   ├── samples.yaml        MC samples, file paths, cross sections
│   ├── objects.yaml        Object quality cuts (muon, electron, btag)
│   ├── acceptance.yaml     Fiducial/preselection cuts
│   ├── regions.yaml        Signal region mass windows
│   ├── cuts.yaml           User-defined selection cuts
│   ├── variables.yaml      Plot variable definitions (bins, ranges, labels)
│   └── triggers.yaml       Trigger definitions + brilcalc paths
├── utils/                  Python utilities
│   ├── data.py             Sample loading, XCache, file discovery
│   ├── plotting.py         Matplotlib/mplhep low-level plot functions
│   ├── triggers.py         Trigger bit definitions (DST, Parking)
│   └── skim.py             Slim ROOT file writer (branch auto-detection)
├── scripts/                Shell scripts
│   ├── run_xsec.sh         GenXSecAnalyzer wrapper
│   ├── run_brilcalc.sh     Luminosity calculation
│   └── check_prescale.sh   Trigger prescale checks
├── legacy/                 Old code (kept for reference, not used)
├── des/                    Design docs (separate git repo)
├── notes.md                Working notes + activity log
└── status.md               This file
```


## Performance

- **Event loop**: Single `RunGraphs()` call batches all Phase 1+2 actions
- **Histogram cache**: `.hist_cache.root` auto-invalidates on script/config/C++ changes
- **Threading**: ROOT ImplicitMT (4 threads default, 32 available)
- **Typical run time**: ~270s for 10 files/sample, 1 variable, skip-cutflow


## What's Next

- [ ] Run with data overlay (drop `--no-data`)
- [ ] Generate all plots for all variables (drop `--plot-vars`)
- [ ] Investigate missing QCD HT bins (400-600, 600-800) — request production?
- [ ] Add leptonic triggers (DST_Muon, DST_Electron, PARKING_Muon, PARKING_EG)
- [ ] Enable BSM signal benchmarks (9 coupling points in config, commented out)
- [x] Performance optimization: merge 7 RunGraphs into 1 (done)
- [x] Generalize `_build_groups()` to be fully config-driven (done — analysis/histograms.py)
- [x] Organize plots into subfolders (stacked/, shape/, overlay/) (done)
- [x] Modularise cutflow_TrigEff.py into analysis/ package (done — 2,280→590 LOC)
