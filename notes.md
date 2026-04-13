# Working Notes — HH→bbττ Scouting Analysis

---

## Activity Log

### 2026-04-09

- **Trigger efficiency study** (`trig_eff.py`): CMS-style trigger efficiency vs gen_mHH plot
  - Filled cyan histogram: signal mHH distribution before trigger
  - Efficiency curves with Clopper-Pearson errors for DST_JetHT and PARKING_HH
  - DST_JetHT: ~45% at mHH=300 GeV, rising to ~100% above 700 GeV
  - PARKING_HH: ~15% at mHH=300 GeV, rising to ~75% at 1000 GeV
  - New function `plot_trigger_efficiency_overlay()` in `utils/plotting.py`
- **Documentation audit and update**: all docs synced with current codebase state
- **Gitignore cleanup**: added `.claude/`, `.vscode/`, `docs/`, `*.root`, `*.npz`, editor/OS files
- Removed Claude co-author tags from all 4 past commits via `git filter-branch`

### 2026-04-07

- **Multi-class BDT** (`classify.py` with `MULTICLASS = True`):
  - 12 decay modes as separate classes (5 DY + 3 TT + 1 QCD + 3 signal)
  - XGBoost `multi:softprob`, per-class weight normalization
  - Stratified train/val/test split preserving class proportions
  - Confusion matrix, per-class ROC (12 one-vs-rest curves), signal score distribution
  - Significance scan on signal score = sum of signal class probabilities
  - Z_A = 0.070 ± 0.014 (with stat uncertainties on S, B, Z_A)
  - QCD class dropped automatically (0 events after preselection)
- **Removed 3 preselection cuts** per advisor: `dphi_MET_tau0 < 1.4`, `MT_tau0_MET < 100`, `D_zeta > -50`; kept as BDT features to check importance
- **Trigger study script** (`trigger_study.py`): per-trigger cutflow + inclusive/exclusive trigger overlays
- **Cutflow name fix**: shows `name: expr` instead of just `name` (e.g. `min_jets: nJets >= 4`)

### 2026-04-05

- **BDT classifier runs** (XGBoost, GPU on A100):
  - Run 1 (with QCD, 5 files): `scale_pos_weight = sum_bkg/sum_sig ≈ 1 billion` → BDT scored everything near 1.0, Z_A = 0.003 despite AUC = 0.999
  - Fixed: normalize training weights so signal total = background total. Physics weights used only for significance scan.
  - Run 2 (no QCD, lr=0.05, 1000 trees, depth=4): AUC = 0.923, Z_A = 0.055
  - Run 3 (no QCD, lr=0.01, 3000 trees, depth=4): AUC = 0.920, Z_A = 0.073
  - Run 4 (depth=10, lr=0.01, 10000 trees): AUC = 0.926, Z_A = 0.073 (early stopped at 3832)
  - Run 5 (all 22 files per sample, depth=10, lr=0.01, 10000 trees): AUC = 0.885, Z_A = 0.063 — **no early stopping fired** (9995/10000), best threshold at 0.995

- **Diagnosis**: `max_depth=10` overfits. Threshold at 0.995 = knife-edge discrimination. Val loss barely moves. Shallower trees (depth=4) generalize better.

- **70/15/15 train/val/test split** implemented with `train_test_split()` from `utils/optimize.py`

- **Training features**: loss curves, tqdm progress bar, GPU acceleration all working

- **Comparison to nominal**: nominal analysis (AN-2025/103) achieves Z_A ≈ 0.245. Scouting BDT is 3-4× worse, expected due to limited reconstruction (no full tracking, coarser calo, UParT proxy for tau ID). Scouting ceiling estimated at Z_A ~ 0.15-0.18.

- **Data leak identified**: `b_coi0_TauVsAll` and `tau_coi0_BvsAll` are defined as `-1.f` for background samples in `define_coi_matching()`. Excluded from BDT features.

- **Skim safety**: identified race condition in auto-skim (no file locking). Added `fcntl.flock()` + atomic temp-file-then-rename to `ensure_slim()`. Added `SKIM = False` flag so scripts can skip auto-skimming entirely (use raw EOS + Define chains).

- **Next steps**: add more features (AK8 H-taggers, subleading pTs, tau_OS, centrality), fix hyperparams (depth=4), k-fold CV, optuna tuning, per-channel BDTs

### 2026-04-04

- With the following cuts, **QCD multijet is reduced to exactly zero**:

```yaml
common:
  min_jets: "nJets >= 2"
  jet_pt: "ak4_pt0 > 20 && ak4_pt1 > 20"
  b_coi0_score: "b_coi0_score > 0.8"
  b_coi1_score: "b_coi1_score > 0.8"
  dphi_bb_tautau: "dphi_bb_tautau > 1.5708"
  dphi_MET_tau0: "dphi_MET_tau0 < 1.5708"
  MT_tau0_MET: "MT_tau0_MET < 100"

tauhtauh:
  tau_coi0_score: "tau_coi0_score > 0.3"
  tau_coi1_score: "tau_coi1_score > 0.3"
  mbb_window: "mbb_coi > 70 && mbb_coi < 150"
  mtautau_window: "mtautau_coi > 50 && mtautau_coi < 150"
```

### 2026-04-01

- Added cumulative S/√B plots (`cuml_sig/` directory, `ratio="cuml_significance"`)
  - Right-to-left integral: for each bin edge x, integrate S and B from x to ∞
  - Both standalone function and parallel `render_task` path
- Fixed parallel plotting: `ProcessPoolExecutor` with `spawn` context
  - Added `if __name__ == '__main__':` guards to all 9 user scripts
  - Moved heavy imports (`analysis.runner`) inside guard so spawn workers stay lightweight
  - Removed top-level `import ROOT` from `utils/plotting.py` (workers don't need ROOT)
  - Emptied `utils/__init__.py` to prevent eager import chain
- Fixed CMS label overlapping between panels in two-panel plots (sig, cuml_sig, eff)
- Fixed `tight_layout` warning for gridspec figures
- Condensed negative-yield warnings to single summary line per trigger (was 20+ lines)
- Extracted `_setup_significance_figure()` helper to deduplicate sig/cuml_sig code
- Cleaned up `skim.py` to reuse `collect_all_expressions()` from `utils/skim.py`
- Created Foam wiki (`wiki/`) with 16 interconnected notes
- Added `.vscode/` config for Foam extension

### 2026-03-27

- Major modularisation refactoring of `cutflow_TrigEff.py` (2,280 → 590 LOC)
- Created `analysis/` package with 7 modules:
  - `config.py` — YAML config loaders (variables, triggers, decay modes)
  - `definitions.py` — all RDataFrame .Define() chains (kinematics, gen matching)
  - `histograms.py` — histogram booking & materialisation
  - `cutflow.py` — cutflow booking, extraction, markdown formatting
  - `cache.py` — histogram cache save/load/invalidation
  - `plots.py` — plot orchestration (Phase 3, 3.5, 5)
- Created `config/variables.yaml` — all PLOT_VARS, GEN_PLOT_VARS, PLOT_VARS_2D,
  DECAY_MODES, CUTFLOW_STEPS moved from inline Python to YAML
- Created `config/triggers.yaml` — TRIG_LIST, exclusive triggers, brilcalc paths
- Updated `des/architecture.md` and `status.md` to reflect new structure

### 2026-03-20

- Ran GenXSecAnalyzer for all 7 QCD-4Jets HT-binned samples — cross sections
  updated in `config/samples.yaml` (were off by 10-100x from placeholder values)
- Added QCD multijet (mode 30) to DECAY_MODES LUT — now 12 MC groups total
  (5 DY + 3 TT + 3 Signal + 1 QCD). QCD skips decayMode filter (no gen matching)
- Cutflow table now shows DY, TT, QCD, Signal columns with S/√B including QCD in B
- Renamed `--max-files` → `--max-mc-files` for clarity (separate from `--max-data-files`)
- Added AK8 fat jet variables (17 new): pT, eta, mass, soft-drop mass,
  ScoutGlobalParT tagger scores (Xbb, Xτhτh, Xτμτh, Xτeτh, QCD)
  for leading + subleading jets. AK8 pT > 150 GeV cut applied.
- Created `des/cpp_engine.md` — full C++ engine reference doc
- Created `notes.md` (activity log) and `status.md` (project summary)
- XCache data loading tested — currently hanging (infrastructure issue)
- Committed and pushed all changes to `scouting` branch

### 2026-03-19

- Implemented `GlobalDecayMode()` C++ function in `elements/GenMatching.C`
  - Unified decay mode numbering: 1-5 (DY), 10-12 (TT), 20-22 (Signal), 30 (QCD)
  - Added NanoAOD type overloads (RVec<short>, RVec<unsigned short>) to fix
    segfaults from implicit type conversion + ImplicitMT
- Fixed `common.h` symbol conflicts with legacy `HHbbtautauAnaElements.C`
  - Moved `deltaPhi`, `deltaR`, `isLastCopy`, `isHardProcess`, `fromHardProcess`
    to anonymous namespaces in GenMatching.C and RecoObjects.C
- Fixed segfault after histogram cache save — added `gc.collect()` +
  `ROOT.gROOT.GetListOfFiles().Clear()` barrier before Phase 3 plotting
- Verified pixi env works end-to-end (compilation + event loop + plotting)

### 2026-03-18 – 2026-03-19

- Restructured `Ana/` for scouting analysis
  - Created `elements/` folder with `common.h`, `GenMatching.C`, `RecoObjects.C`
  - Moved legacy code (`HHbbtautauAnaElements.C`, `Particle.h`) to `legacy/`
  - Removed legacy Jupyter notebooks from `legacy/`
- Added `output/` and `logs/` to `.gitignore`
- Set up pixi environment at `/work/users/das214/pixi/ana/`
  - Pixi can't build under `/home/` (Purdue AF policy)
  - Workflow: `pixi shell` from `/work/`, then `cd` to `/home/.../Ana/`
- Added trigger overlay plots (Phase 3.5) and per-channel signal overlays
- Added `--plot-vars`, `--recache`, `--skip-cutflow` CLI flags
- Implemented histogram cache (`.hist_cache.root`) with auto-invalidation
- Added NoTrigger (preselection only) alongside DST_JetHT and PARKING_HH
- Created `make_pdf.py` comparison tool (side-by-side triggers, `--compare` mode)
- YAML configs: `config/objects.yaml`, `config/acceptance.yaml`, `config/regions.yaml`

### 2026-03-11

- Updated luminosity to brilcalc value (103.965 fb⁻¹)
- Fixed b-tagger: BvsAll excludes `prob_bb` (AK8-only variable)
- Added `--max-mc-files` for data file limiting (replaces Range() which breaks ImplicitMT)
- Created PDF maker tool

### 2026-03-08 – 2026-03-09

- Added kinematic plots (jet pT, eta, phi, HT, m4j, dijet masses)
- Implemented b-jet selection with ParticleNet BvsAll discriminant
- Added dark theme support
- Enabled ROOT ImplicitMT (32 threads)

### 2026-02-26 – 2026-02-27

- Refactored codebase, added GenXSecAnalyzer scripts (`scripts/run_xsec.sh`)
- Added Run2024 scouting data
- Fixed yield plot errors

### 2026-02-06 – 2026-02-13

- Initial scouting analysis: trigger efficiency, cutflow
- AK4 and AK8 jet implementations
- Added acceptance conditions and `.gitignore`
- Cleaned up `__pycache__` and `.ipynb_checkpoints`

### 2026-01-19

- Gen-level matching work (on `lite` branch, pre-scouting)
  - bb and ττ component matching, dR requirements
  - FatJet matching, tagger variables
  - Electron reco function, ID thresholds

### 2025 (lite branch — pre-scouting era)

- ABCD method implementation, signal region cuts
- C++ anaConfig for cut definitions
- Fit distribution exports
- Repository pruning (Nov 2025)

### 2022 (original repo — BDT/ML era)

- BDT training with variable importance
- ROC curves, FOM computation
- Input variable correlations

---

## Environment Setup

- **Pixi env**: `/work/users/das214/pixi/ana/` (activate: `cd /work/users/das214/pixi/ana && pixi shell`)
- **CMSSW env**: `source /cvmfs/cms.cern.ch/cmsset_default.sh && cd CMSSW_15_0_15/src && eval \`scramv1 runtime -sh\``
- **Never mix** CMSSW and pixi in the same shell — CMSSW pollutes LD_LIBRARY_PATH and causes segfaults
- `export PIXI_CACHE_DIR=/work/users/das214/.pixi-cache` (add to `~/.bashrc`)

## QCD Cross Sections (GenXSecAnalyzer, 2026-03-20)

Ran via `scripts/run_xsec.sh` on MINIAODSIM, 1M events per sample.

| HT bin [GeV] | xsec [pb] | uncertainty [pb] |
|---------------|-----------|-------------------|
| 40-70 | 3.113e+08 | ± 1.042e+06 |
| 100-200 | 2.513e+07 | ± 3.436e+05 |
| 200-400 | 1.975e+06 | ± 8.945e+03 |
| 800-1000 | 3.027e+03 | ± 1.426e+01 |
| 1000-1200 | 8.951e+02 | ± 4.858e+00 |
| 1500-2000 | 1.271e+02 | ± 8.289e-01 |
| 2000+ | 2.672e+01 | ± 4.135e-01 |

Missing HT bins: 400-600 and 600-800 (not produced on EOS).

## Global Decay Mode Numbering

```
 0  = unknown / not classified
 1  = DY  Z→ee
 2  = DY  Z→μμ
 3  = DY  Z→τhτh
 4  = DY  Z→τμτh
 5  = DY  Z→τeτh
10  = TT  fully hadronic
11  = TT  semi-leptonic
12  = TT  fully leptonic
20  = HH  →bb τhτh
21  = HH  →bb τμτh
22  = HH  →bb τeτh
30  = QCD multijet (no gen-level decay matching)
```

## CLI Quick Reference

```bash
# Full run (all vars, all plot types)
python cutflow_TrigEff.py --theme light --max-mc-files 10 --no-data --nthreads 4

# Quick test (1 var, shapes only, skip cutflow)
python cutflow_TrigEff.py --theme light --max-mc-files 10 --no-data --nthreads 4 \
    --plot-vars "ak4_pt0" --plot-type shape --skip-cutflow

# With data overlay (drop --no-data, uses XCache)
python cutflow_TrigEff.py --theme light --max-mc-files 10 --nthreads 4 \
    --recache --overwrite --plot-vars "ak4_pt0"

# Force rebuild everything
python cutflow_TrigEff.py --theme light --max-mc-files 10 --no-data --nthreads 4 \
    --recache --overwrite

# GenXSecAnalyzer (needs cmsenv in separate shell, NOT pixi)
cd scripts/ && source run_xsec.sh <index>
```

## Compilation Notes

- `root -l -q -e 'gROOT->LoadMacro("elements/GenMatching.C+")'` — standalone test (from Ana/ dir)
- `root++ ` (double plus) forces recompile; single `+` reuses cached .so
- ACLiC `cd`s into the file's directory — use `#include "common.h"` not `"elements/common.h"`
- Never run standalone LoadMacro when cutflow_TrigEff.py uses `SetBuildDir()` — creates conflicting .so files

## Known Issues

- Scouting ParticleNet b-tagging is fundamentally poor (limited tracking info, wp=0.1)
- `Range()` incompatible with `EnableImplicitMT()` — use `--max-mc-files` instead
- Segfault after cache save with 10 files — fixed with `gc.collect()` + `ROOT.gROOT.GetListOfFiles().Clear()`
- DY PTLL-600 has 3080 files (large stats) — may dominate memory usage
