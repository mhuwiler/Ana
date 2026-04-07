#!/usr/bin/env python3
"""Cut optimization for tauhtauh channel using RGS and Grid Search.

Strategy (following standard CMS practice):
  - Mass windows fixed from physics (mH = 125 GeV) → applied as preselection
  - Tagger working points optimized to maximize Asimov significance Z_A

Methods:
  1. RGS            — signal events as candidate thresholds (sequential)
  2. RGS Parallel   — same, split across CPU cores
  3. Grid Search    — evenly spaced 4D grid (sequential)
  4. Grid Parallel  — same, split across CPU cores
"""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT", "QCD"]
TRIGGERS     = ["NoTrigger"]
DATA         = False
SKIM         = True
NTHREADS     = 4
MAX_MC_FILES = 22
MAX_EVENTS   = 500000

# Preselection: all cuts from cuts.yaml EXCEPT the tagger scores (which we optimize)
# CUTS = [
#     "nJets >= 4",
#     "ak4_pt0 > 20 && ak4_pt1 > 20",
#     "dphi_bb_tautau > 1.5708",
#     "dphi_MET_tau0 < 1.5708",
#     "MT_tau0_MET < 100",
#     "D_zeta > -50",
#     "!has_good_muon && !has_good_electron",   # lepton veto (tauhtauh)
#     "mbb_coi > 70 && mbb_coi < 150",
#     "mtautau_coi > 50 && mtautau_coi < 150",
# ]


# # Preselection (tauhtauh channel, minus tagger scores — BDT handles those)
CUTS = [
    "nJets >= 4",
    "ak4_pt0 > 20 && ak4_pt1 > 20",
    "dphi_bb_tautau > 1.5",
    "dphi_MET_tau0 < 1.4",
    # TODO: "dphi_MET_tau1 < 1.4", (form the second tau) (this can go for the leptonic tau)
    "MT_tau0_MET < 100", # This also does not go fo rthe hadronic tau (check)
    "D_zeta > -50",
    # "!has_good_muon && !has_good_electron",
    "mbb_coi > 70 && mbb_coi < 150",
    "mtautau_coi > 50 && mtautau_coi < 150",
]




# Tagger scores to optimize: (column_name, direction, (scan_lo, scan_hi))
# All one-sided ">" cuts — RGS uses signal events as candidate thresholds
OPT_VARS = [
    ("b_coi0_score",    ">", (0.0, 1.0)),
    ("b_coi1_score",    ">", (0.0, 1.0)),
    ("tau_coi0_score",  ">", (0.0, 1.0)),
    ("tau_coi1_score",  ">", (0.0, 1.0)),
]

# Optimization settings
N_RGS_SAMPLES = 50000     # number of signal events to sample for RGS
N_GRID_STEPS  = 20         # steps per variable for grid search (20^4 = 160k points)
N_WORKERS     = 64         # parallel workers (numpy fallback only; JIT/GPU handle parallelism internally)
LUMI          = 103.965    # fb^-1
SEED          = 42
BACKEND       = "auto"     # "auto" (GPU→JIT→numpy), "gpu", "jit", "numpy"

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import os
    import numpy as np
    from analysis.runner import setup, load_and_run
    from utils.optimize import (
        extract_arrays, train_test_split, evaluate_on_test,
        run_rgs, run_rgs_parallel,
        run_grid, run_grid_parallel,
        print_results, plot_fom_ranking, save_results,
    )

    ctx = setup(__file__, samples=SAMPLES, triggers=TRIGGERS, cuts=CUTS, data=DATA,
                skim=SKIM, nthreads=NTHREADS, max_mc_files=MAX_MC_FILES, max_events=MAX_EVENTS)
    result = load_and_run(ctx, [])

    # Extract numpy arrays from RDataFrames (post mass-window preselection)
    print("\nExtracting per-event arrays via AsNumpy (this may take ~30s)...", flush=True)
    mc_sel = result.trig_selections["NoTrigger"]["mc_sel"]
    sig, sig_w, bkg, bkg_w = extract_arrays(mc_sel, OPT_VARS)
    print(f"\nSignal: {len(sig)} events, Background: {len(bkg)} events")
    print(f"Total signal weight (× lumi): {np.sum(sig_w) * LUMI:.2f}, "
          f"Total bkg weight (× lumi): {np.sum(bkg_w) * LUMI:.2g}")

    # Train/test split (50/50)
    print("\nSplitting into train/test sets...")
    (sig_train, sig_w_train, bkg_train, bkg_w_train,
     sig_test, sig_w_test, bkg_test, bkg_w_test) = train_test_split(
        sig, sig_w, bkg, bkg_w, test_fraction=0.5, seed=SEED)

    # ── RGS (parallel / GPU / JIT) — optimize on TRAIN ──
    best_rgs = run_rgs_parallel(sig_train, sig_w_train, bkg_train, bkg_w_train,
                                OPT_VARS, LUMI, n_samples=N_RGS_SAMPLES, seed=SEED,
                                n_workers=N_WORKERS, backend=BACKEND)
    evaluate_on_test(best_rgs, sig_test, sig_w_test, bkg_test, bkg_w_test, LUMI)
    print_results(best_rgs, OPT_VARS)

    # ── Grid Search (parallel / GPU / JIT) — optimize on TRAIN ──
    best_grid = run_grid_parallel(sig_train, sig_w_train, bkg_train, bkg_w_train,
                                  OPT_VARS, LUMI, n_steps=N_GRID_STEPS,
                                  n_workers=N_WORKERS, backend=BACKEND)
    evaluate_on_test(best_grid, sig_test, sig_w_test, bkg_test, bkg_w_test, LUMI)
    print_results(best_grid, OPT_VARS)

    # ── Save results + plots ──
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "optimize")
    plots_dir = os.path.join(out_dir, "plots")
    save_results({"RGS": best_rgs, "Grid": best_grid}, OPT_VARS, out_dir)
    plot_fom_ranking(best_rgs, os.path.join(plots_dir, "rgs_fom.png"))
    plot_fom_ranking(best_grid, os.path.join(plots_dir, "grid_fom.png"))
