# HH→bbττ Scouting Analysis — Current Status

*Last updated: 2026-04-09*

## Overview

Search for non-resonant Higgs pair production in the bbττ final state using
CMS Run 3 (2024) scouting data at √s = 13.6 TeV. Three decay channels:
τhτh (fully hadronic), τμτh, and τeτh. Uses UParT AK4 tagger for b/tau
identification and COI (Candidate-of-Interest) matching for Higgs candidate
reconstruction.

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

Per-channel cuts defined in `config/cuts.yaml` with named cards:

**Common**: nJets ≥ 2, jet pT > 20, b-tag scores > 0.8, dphi(bb,ττ) > 1.5, dphi(MET,τ₀) < 1.4, MT(τ₀,MET) < 100, D_ζ > -50

**τhτh**: nJets ≥ 4, lepton veto, tau scores > 0.3, mbb ∈ [70,150], mττ ∈ [50,150]

**τμτh**: nJets ≥ 3, require good muon, ΔR(μ,τ) > 0.5, tau score > 0.3, mass windows

**τeτh**: nJets ≥ 3, require good electron, ΔR(e,τ) > 0.5, tau score > 0.3, mass windows

**QCD rejection variables**: D_ζ, MT(τ₀,MET), dphi(MET,τ₀), score_product, tau_OS, MET_significance

**Lepton selection**: `Ana::SelectMuon()` / `Ana::SelectElectron()` in C++, cuts from `config/objects.yaml`

Four trigger paths:
- **NoTrigger** — preselection only
- **DST_JetHT** — scouting hadronic
- **DST_Muon** / **DST_Electron** — scouting leptonic
- **PARKING_HH** — parking HH (11 HLT paths)


## Plot Types

| Type | Description | Output |
|------|-------------|--------|
| Stacked | MC stacked histogram (12 groups) | `stacked/{trig}_{var}.png` |
| Shape | Normalized MC overlay | `shape/{trig}_{var}.png` |
| Trigger overlay | All-MC summed, triggers overlaid | `overlay/{var}.png` |
| Exclusive overlay | Events firing ONLY one trigger | `overlay_exclusive/{var}.png` |
| Per-channel overlay | Signal-only per decay channel | `overlay/{ch}_{var}.png` |
| Eff stacked | Stacked + trigger efficiency panel | `eff/{trig}/{var}.png` |
| Sig stacked | Stacked + bin-by-bin S/√B panel | `sig/{trig}/{var}.png` |
| Cum sig | Stacked + cumulative S/√B (right-to-left) | `cuml_sig/{trig}/{var}.png` |
| Trig efficiency | CMS-style efficiency vs gen_mHH | `trigger_efficiency_mHH.png` |

## BDT Classification

| Mode | Description | Objective |
|------|-------------|-----------|
| Binary | Signal vs background | `binary:logistic` |
| Multi-class | 12 decay modes as separate classes | `multi:softprob` |

**BDT plots**: confusion matrix, per-class ROC, feature importance, score distribution, significance scan (with stat uncertainties), overtraining check, loss curves.

**Best results** (binary, tauhtauh preselection): AUC = 0.93, Z_A = 0.087. Nominal analysis achieves Z_A ≈ 0.245 — gap is inherent to scouting reconstruction.

## Cut Optimization

RGS (Random Grid Search) + Grid Search with GPU/JIT/numpy backends. Signal events as candidate thresholds, Asimov significance as figure of merit.


## Code Architecture

14 user scripts, 11 analysis modules, 8 utility modules:

```
scouting_HHbbtautau/
├── User scripts (analysis)
│   ├── classify.py             BDT classifier (binary + multi-class)
│   ├── optimize.py             Cut optimization (RGS + Grid Search)
│   ├── trigger_study.py        Trigger overlays + cutflow per trigger
│   ├── trig_eff.py             Trigger efficiency vs gen_mHH
│   ├── coi_study.py            COI tagging study (with cuts)
│   ├── HHbbtt.py               Gen-matched analysis
│   ├── kinematics.py           AK4 jet kinematics
│   ├── ak8.py                  AK8 fat-jet studies
│   └── ...                     (+ 6 more: allMC_wocuts, coi_study_wocuts, etc.)
├── analysis/                   Framework core (11 modules)
│   ├── runner.py               setup(), load_and_run(), PlotResult, Plotter
│   ├── event_loop.py           book_and_run(), RunGraphs, EventLoopResult
│   ├── data_loading.py         MC/data loading, SKIM flag, auto-skim
│   ├── definitions.py          RDataFrame .Define() chains (~42 KB)
│   ├── cutflow.py              Per-trigger cutflow with per-channel S/√B
│   ├── histograms.py           Histogram booking (inclusive + exclusive triggers)
│   ├── plots.py                Plot orchestration + parallel rendering
│   ├── config.py               YAML config loaders
│   ├── constants.py            Branch names, physics constants
│   ├── variables.py            PlotVar dataclass
│   └── setup.py                ROOT init, macro loading
├── utils/                      Library modules (8 modules)
│   ├── classify.py             BDT: train, evaluate, significance scan, plots
│   ├── optimize.py             Cut optimizer: RGS, Grid, GPU/JIT/numpy
│   ├── plotting.py             Stacked/shape/overlay/efficiency plots, CMS style
│   ├── skim.py                 Slim files with file locking + atomic writes
│   ├── data.py                 Sample loading, brilcalc, file discovery
│   ├── triggers.py             Trigger bit definitions
│   ├── logging.py              Timestamped run logging
│   └── pdf.py                  PDF generation
├── elements/                   C++ engine
│   ├── common.h                Shared types, PDG constants
│   ├── GenMatching.C           Gen matching, GlobalDecayMode()
│   └── RecoObjects.C           SelectMuon(), SelectElectron(), dijet pairing
├── config/                     YAML configuration
│   ├── samples.yaml            MC samples, cross sections, file paths
│   ├── analysis.yaml           Decay modes, sig/bkg mode lists
│   ├── cuts.yaml               Named cut cards (common, tauhtauh, taumutauh, tauetauh)
│   ├── objects.yaml            Lepton quality cuts (pT, η, iso, dxy, dz)
│   └── triggers.yaml           Trigger definitions + brilcalc paths
├── legacy/                     Old code (kept for reference)
├── docs/                       Wiki + design docs
├── notes.md                    Working notes + activity log
└── status.md                   This file
```


## Performance

- **Event loop**: Single `RunGraphs()` call batches all Phase 1+2 actions
- **Histogram cache**: `.hist_cache.root` auto-invalidates on script/config/C++ changes
- **Threading**: ROOT ImplicitMT (4 threads default, 32 available)
- **Typical run time**: ~270s for 10 files/sample, 1 variable, skip-cutflow


## What's Next

- [ ] N-1 plots (show each cut's impact individually)
- [ ] Stat uncertainties in cutflow table (√Σw²)
- [ ] Per-channel BDTs (separate for τhτh, τμτh, τeτh)
- [ ] mHH reconstruction (4-body invariant mass from COI candidates)
- [ ] DNN classifier (PyTorch, compare with BDT)
- [ ] Bayesian hyperparameter optimization (Optuna)
- [ ] k-fold cross-validation for more reliable Z_A
- [ ] Enable BSM signal benchmarks (9 coupling points in config)
- [x] Multi-class BDT (12 decay modes, XGBoost multi:softprob)
- [x] Binary BDT classifier (XGBoost GPU, Z_A = 0.087)
- [x] Cut optimizer (RGS + Grid Search, GPU/JIT/numpy)
- [x] Trigger efficiency study (CMS-style mHH plot)
- [x] Trigger overlays (inclusive + exclusive)
- [x] Lepton selection (C++ SelectMuon/SelectElectron)
- [x] QCD rejection (D_ζ, MT, dphi, score_product → QCD = 0)
- [x] Skim safety (file locking, SKIM=False flag)
- [x] Per-trigger cuts with per-channel cutflow
- [x] Cumulative S/√B plots
- [x] Modularise into analysis/ package
- [x] Merge RunGraphs into single event loop
