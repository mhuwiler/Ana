# Working Notes — HH→bbττ Scouting Analysis

---

## Activity Log

### 2026-03-20

- Ran GenXSecAnalyzer for all 7 QCD-4Jets HT-binned samples — cross sections
  updated in `config/samples.yaml` (were off by 10-100x from placeholder values)
- Added QCD multijet (mode 30) to DECAY_MODES LUT — now 12 MC groups total
  (5 DY + 3 TT + 3 Signal + 1 QCD). QCD skips decayMode filter (no gen matching)
- Cutflow table now shows DY, TT, QCD, Signal columns with S/√B including QCD in B
- Created `des/cpp_engine.md` — full C++ engine reference doc
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
