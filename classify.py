#!/usr/bin/env python3
"""BDT classifier for HH→bbττ signal vs background discrimination.

Trains an XGBoost BDT on ~23 features after preselection cuts.
70/15/15 train/val/test split. GPU-accelerated if available.
"""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT", "QCD"]
TRIGGERS     = ["NoTrigger"]
DATA         = False
SKIM         = True
NTHREADS     = 32
MAX_MC_FILES = 22
MAX_EVENTS   = None

# Preselection (tauhtauh channel, minus tagger scores — BDT handles those)
CUTS = [
    "nJets >= 4",
    # "ak4_pt0 > 20 && ak4_pt1 > 20",
    # "dphi_bb_tautau > 1.5708",
    # "dphi_MET_tau0 < 1.5708",
    # "MT_tau0_MET < 100",
    # "D_zeta > -50",
    # "!has_good_muon && !has_good_electron",
    # "mbb_coi > 70 && mbb_coi < 150",
    # "mtautau_coi > 50 && mtautau_coi < 150",
]

# Features for BDT (23 variables)
FEATURES = [
    "b_coi0_score", "b_coi1_score", "tau_coi0_score", "tau_coi1_score",
    "mbb_coi", "mtautau_coi",
    "dR_bb_coi", "dR_tautau_coi", "dR_bb_tautau", "dphi_bb_tautau",
    "HT", "MET", "MET_significance", "nJets",
    "score_product", "D_zeta", "MT_tau0_MET",
    "b_coi0_pt", "tau_coi0_pt",
]

LUMI     = 103.965
SEED     = 42
USE_GPU  = True

# BDT hyperparameters (override defaults from utils/classify.py)
BDT_PARAMS = {
    "max_depth": 10,
    "learning_rate": 0.01,
    "n_estimators": 10_000,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 10,
    "gamma": 1.0,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
    "early_stopping_rounds": 50,
}

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import os
    import numpy as np
    from analysis.runner import setup, load_and_run
    from utils.classify import (
        extract_features, train_bdt, significance_scan,
        plot_roc, plot_feature_importance, plot_score_distribution,
        plot_significance_scan, plot_overtraining, plot_loss_curves, save_model,
    )
    from utils.optimize import train_test_split
    from sklearn.metrics import roc_curve, roc_auc_score

    # 1. Setup + load
    ctx = setup(__file__, samples=SAMPLES, triggers=TRIGGERS, cuts=CUTS,
                data=DATA, skim=SKIM, nthreads=NTHREADS, max_mc_files=MAX_MC_FILES,
                max_events=MAX_EVENTS)
    result = load_and_run(ctx, [])

    # 2. Extract features
    print("\nExtracting features via AsNumpy...")
    mc_sel = result.trig_selections["NoTrigger"]["mc_sel"]
    X_sig, w_sig, X_bkg, w_bkg = extract_features(mc_sel, FEATURES)
    print(f"\nSignal: {len(X_sig)} events, Background: {len(X_bkg)} events")
    print(f"Signal weighted (× lumi): {np.sum(w_sig) * LUMI:.2f}")
    print(f"Background weighted (× lumi): {np.sum(w_bkg) * LUMI:.2g}")

    # 3. Train/Val/Test split (70/15/15)
    print("\nSplitting 70/15/15 (train/val/test)...")
    # First: 85% train+val, 15% test
    (sig_tv, w_sig_tv, bkg_tv, w_bkg_tv,
     sig_te, w_sig_te, bkg_te, w_bkg_te) = train_test_split(
        X_sig, w_sig, X_bkg, w_bkg, test_fraction=0.15, seed=SEED)
    # Second: from the 85%, take ~17.6% as val → 70% train, 15% val
    (sig_tr, w_sig_tr, bkg_tr, w_bkg_tr,
     sig_val, w_sig_val, bkg_val, w_bkg_val) = train_test_split(
        sig_tv, w_sig_tv, bkg_tv, w_bkg_tv, test_fraction=0.176, seed=SEED+1)

    # Combine into X, y, w
    X_train = np.vstack([sig_tr, bkg_tr])
    y_train = np.concatenate([np.ones(len(sig_tr)), np.zeros(len(bkg_tr))])
    w_train = np.concatenate([w_sig_tr, w_bkg_tr])

    X_val = np.vstack([sig_val, bkg_val])
    y_val = np.concatenate([np.ones(len(sig_val)), np.zeros(len(bkg_val))])
    w_val = np.concatenate([w_sig_val, w_bkg_val])

    X_test = np.vstack([sig_te, bkg_te])
    y_test = np.concatenate([np.ones(len(sig_te)), np.zeros(len(bkg_te))])
    w_test = np.concatenate([w_sig_te, w_bkg_te])

    print(f"  Train: {len(X_train)} ({len(sig_tr)} sig + {len(bkg_tr)} bkg)")
    print(f"  Val:   {len(X_val)} ({len(sig_val)} sig + {len(bkg_val)} bkg)")
    print(f"  Test:  {len(X_test)} ({len(sig_te)} sig + {len(bkg_te)} bkg)")

    # 4. Train BDT
    print("\nTraining BDT...")
    model, eval_results = train_bdt(X_train, y_train, w_train,
                                     X_val, y_val, w_val,
                                     FEATURES, use_gpu=USE_GPU, params=BDT_PARAMS)

    # 5. Score test set
    scores_test = model.predict_proba(X_test)[:, 1]
    scores_sig = scores_test[y_test == 1]
    scores_bkg = scores_test[y_test == 0]

    # 6. ROC + AUC on test set
    fpr, tpr, _ = roc_curve(y_test, scores_test, sample_weight=np.abs(w_test))
    auc_val = roc_auc_score(y_test, scores_test, sample_weight=np.abs(w_test))
    print(f"\n  Test AUC = {auc_val:.4f}")

    # 7. Significance scan on TEST set (using PHYSICS weights, not normalized)
    print("\nSignificance scan on test set...")
    scan = significance_scan(scores_sig, w_sig_te, scores_bkg, w_bkg_te, LUMI)
    print(f"  Best threshold: score > {scan['best_threshold']:.4f}")
    print(f"  S = {scan['best_S']:.4g},  B = {scan['best_B']:.4g},  Z_A = {scan['best_Z']:.4f}")

    # 8. Plots + save
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "classify")
    plots_dir = os.path.join(out_dir, "plots")

    plot_roc(fpr, tpr, auc_val, os.path.join(plots_dir, "roc.png"))
    plot_feature_importance(model, FEATURES, os.path.join(plots_dir, "importance.png"))
    plot_score_distribution(scores_sig, w_sig_te, scores_bkg, w_bkg_te,
                            os.path.join(plots_dir, "score_dist.png"))
    plot_significance_scan(scan, os.path.join(plots_dir, "significance_scan.png"))
    plot_overtraining(model, X_train, y_train, X_test, y_test,
                      os.path.join(plots_dir, "overtraining.png"))
    plot_loss_curves(eval_results, getattr(model, 'best_iteration', None),
                     os.path.join(plots_dir, "loss_curves.png"))
    save_model(model, FEATURES, scan, auc_val, out_dir)

    print(f"\n{'=' * 60}")
    print(f"  BDT Classification Complete")
    print(f"{'=' * 60}")
    print(f"  AUC = {auc_val:.4f}")
    print(f"  Best cut: BDT score > {scan['best_threshold']:.4f}")
    print(f"  S = {scan['best_S']:.4g},  B = {scan['best_B']:.4g},  Z_A = {scan['best_Z']:.4f}")
    print(f"  Results saved to: {out_dir}")
