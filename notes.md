# Working Notes

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
python cutflow_TrigEff.py --theme light --max-files 10 --no-data --nthreads 4

# Quick test (1 var, shapes only, skip cutflow)
python cutflow_TrigEff.py --theme light --max-files 10 --no-data --nthreads 4 \
    --plot-vars "ak4_pt0" --plot-type shape --skip-cutflow

# Force rebuild everything
python cutflow_TrigEff.py --theme light --max-files 10 --no-data --nthreads 4 \
    --recache --overwrite

# GenXSecAnalyzer (needs cmsenv, NOT pixi)
cd scripts/ && source run_xsec.sh <index>
```

## Compilation Notes

- `root -l -q -e 'gROOT->LoadMacro("elements/GenMatching.C+")'` — standalone test (from Ana/ dir)
- `root++ ` (double plus) forces recompile; single `+` reuses cached .so
- ACLiC `cd`s into the file's directory — use `#include "common.h"` not `"elements/common.h"`
- Never run standalone LoadMacro when cutflow_TrigEff.py uses `SetBuildDir()` — creates conflicting .so files

## Known Issues

- Scouting ParticleNet b-tagging is fundamentally poor (limited tracking info, wp=0.1)
- `Range()` incompatible with `EnableImplicitMT()` — use `--max-files` instead
- Segfault after cache save with 10 files — fixed with `gc.collect()` + `ROOT.gROOT.GetListOfFiles().Clear()`
- DY PTLL-600 has 3080 files (large stats) — may dominate memory usage
