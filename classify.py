#!/usr/bin/env python3
"""BDT classifier for HH→bbττ signal vs background discrimination.

Trains an XGBoost BDT on ~23 features after preselection cuts.
70/15/15 train/val/test split. GPU-accelerated if available.
"""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT"] #, "QCD"]
TRIGGERS     = ["NoTrigger"]
DATA         = False
SKIM         = True
NTHREADS     = 32
MAX_MC_FILES = 22
MAX_EVENTS   = None

# Preselection
CUTS = [
    "nJets >= 4",
    "ak4_pt0 > 20 && ak4_pt1 > 20",
    "dphi_bb_tautau > 1.5",
    "mbb_coi > 70 && mbb_coi < 150",
    "mtautau_coi > 50 && mtautau_coi < 150",
]

# Features for BDT (20 variables)
FEATURES = [
    "b_coi0_score", "b_coi1_score", "tau_coi0_score", "tau_coi1_score",
    "mbb_coi", "mtautau_coi",
    "dR_bb_coi", "dR_tautau_coi", "dR_bb_tautau", "dphi_bb_tautau",
    "HT", "MET", "MET_significance", "nJets",
    "score_product",
    "D_zeta", "MT_tau0_MET", "dphi_MET_tau0",
    "b_coi0_pt", "tau_coi0_pt",
]

LUMI     = 103.965
SEED     = 42
USE_GPU  = True

# Multi-class classification: all 12 decay modes as separate classes
MULTICLASS = False

CLASS_MAP = {
    1: 0, 2: 1, 3: 2, 4: 3, 5: 4,       # DY (5 classes)
    10: 5, 11: 6, 12: 7,                  # TT (3 classes)
    30: 8,                                 # QCD (1 class)
    20: 9, 21: 10, 22: 11,                # Signal (3 classes)
}
CLASS_NAMES = [
    "DY ee", "DY \u03bc\u03bc", "DY \u03c4h\u03c4h", "DY \u03c4\u03bc\u03c4h", "DY \u03c4e\u03c4h",
    "TT had", "TT semi", "TT dilep",
    "QCD",
    "HH \u03c4h\u03c4h", "HH \u03c4\u03bc\u03c4h", "HH \u03c4e\u03c4h",
]
SIG_CLASSES = [9, 10, 11]

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
        extract_features_multiclass, train_bdt_multiclass, significance_scan_multiclass,
        plot_roc, plot_feature_importance, plot_score_distribution,
        plot_significance_scan, plot_overtraining, plot_loss_curves, save_model,
        plot_confusion_matrix_multiclass, plot_roc_multiclass,
        plot_score_distribution_multiclass, plot_overtraining_multiclass,
    )
    from utils.optimize import train_test_split, train_test_split_multiclass
    from sklearn.metrics import roc_curve, roc_auc_score

    # 1. Setup + load
    ctx = setup(__file__, samples=SAMPLES, triggers=TRIGGERS, cuts=CUTS,
                data=DATA, skim=SKIM, nthreads=NTHREADS, max_mc_files=MAX_MC_FILES,
                max_events=MAX_EVENTS)
    result = load_and_run(ctx, [])

    mc_sel = result.trig_selections["NoTrigger"]["mc_sel"]
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "classify")
    plots_dir = os.path.join(out_dir, "plots")

    if MULTICLASS:
        # ── Multi-class flow ──────────────────────────────────────────────
        n_classes = len(CLASS_NAMES)

        # 2. Extract features with class labels
        print("\nExtracting features via AsNumpy (multi-class)...")
        X, y, w, present_classes = extract_features_multiclass(mc_sel, FEATURES, CLASS_MAP)

        # Remap class names and signal indices to match contiguous labels
        active_names = [CLASS_NAMES[c] for c in present_classes]
        active_sig = [i for i, c in enumerate(present_classes) if c in SIG_CLASSES]
        n_classes = len(active_names)
        print(f"  Active classes ({n_classes}): {active_names}")
        print(f"  Signal class indices: {active_sig}")

        # 3. Stratified 70/15/15 split
        print("\nSplitting 70/15/15 (train/val/test, stratified)...")
        X_tv, y_tv, w_tv, X_te, y_te, w_te = train_test_split_multiclass(
            X, y, w, test_fraction=0.15, seed=SEED)
        X_tr, y_tr, w_tr, X_val, y_val, w_val = train_test_split_multiclass(
            X_tv, y_tv, w_tv, test_fraction=0.176, seed=SEED+1)
        print(f"  Train: {len(X_tr)}, Val: {len(X_val)}, Test: {len(X_te)}")

        # 4. Train multi-class BDT
        print("\nTraining multi-class BDT...")
        model, eval_results = train_bdt_multiclass(
            X_tr, y_tr, w_tr, X_val, y_val, w_val,
            FEATURES, n_classes=n_classes, class_names=active_names,
            use_gpu=USE_GPU, params=BDT_PARAMS)

        # 5. Score test set
        proba_test = model.predict_proba(X_te)
        sig_score = sum(proba_test[:, c] for c in active_sig)
        y_pred = np.argmax(proba_test, axis=1)

        # 6. Significance scan (using PHYSICS weights)
        print("\nSignificance scan on test set...")
        scan = significance_scan_multiclass(proba_test, w_te, y_te, LUMI, active_sig)
        print(f"  Best threshold: sig_score > {scan['best_threshold']:.4f}")
        print(f"  S = {scan['best_S']:.4g} +/- {scan['best_dS']:.4g}")
        print(f"  B = {scan['best_B']:.4g} +/- {scan['best_dB']:.4g}")
        print(f"  Z_A = {scan['best_Z']:.4f} +/- {scan['best_dZ']:.4f}")

        # 7. Plots
        plot_confusion_matrix_multiclass(y_te, y_pred, active_names,
                                          os.path.join(plots_dir, "confusion_matrix.png"))
        plot_roc_multiclass(y_te, proba_test, active_names,
                            os.path.join(plots_dir, "roc.png"))
        plot_feature_importance(model, FEATURES, os.path.join(plots_dir, "importance.png"))
        plot_score_distribution_multiclass(sig_score, w_te, y_te, active_names, active_sig,
                                           os.path.join(plots_dir, "score_dist.png"))
        plot_significance_scan(scan, os.path.join(plots_dir, "significance_scan.png"))
        plot_overtraining_multiclass(model, X_tr, y_tr, X_te, y_te, active_sig,
                                      os.path.join(plots_dir, "overtraining.png"))
        plot_loss_curves(eval_results, getattr(model, 'best_iteration', None),
                         os.path.join(plots_dir, "loss_curves.png"))
        save_model(model, FEATURES, scan, None, out_dir)

        print(f"\n{'=' * 60}")
        print(f"  Multi-class BDT Classification Complete ({n_classes} classes)")
        print(f"{'=' * 60}")
        print(f"  Best cut: sig_score > {scan['best_threshold']:.4f}")
        print(f"  S = {scan['best_S']:.4g} +/- {scan['best_dS']:.4g}")
        print(f"  B = {scan['best_B']:.4g} +/- {scan['best_dB']:.4g}")
        print(f"  Z_A = {scan['best_Z']:.4f} +/- {scan['best_dZ']:.4f}")
        print(f"  Results saved to: {out_dir}")

    else:
        # ── Binary flow (unchanged) ──────────────────────────────────────
        # 2. Extract features
        print("\nExtracting features via AsNumpy...")
        X_sig, w_sig, X_bkg, w_bkg = extract_features(mc_sel, FEATURES)
        print(f"\nSignal: {len(X_sig)} events, Background: {len(X_bkg)} events")
        print(f"Signal weighted (x lumi): {np.sum(w_sig) * LUMI:.2f}")
        print(f"Background weighted (x lumi): {np.sum(w_bkg) * LUMI:.2g}")

        # 3. Train/Val/Test split (70/15/15)
        print("\nSplitting 70/15/15 (train/val/test)...")
        (sig_tv, w_sig_tv, bkg_tv, w_bkg_tv,
         sig_te, w_sig_te, bkg_te, w_bkg_te) = train_test_split(
            X_sig, w_sig, X_bkg, w_bkg, test_fraction=0.15, seed=SEED)
        (sig_tr, w_sig_tr, bkg_tr, w_bkg_tr,
         sig_val, w_sig_val, bkg_val, w_bkg_val) = train_test_split(
            sig_tv, w_sig_tv, bkg_tv, w_bkg_tv, test_fraction=0.176, seed=SEED+1)

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

        # 6. ROC + AUC
        fpr, tpr, _ = roc_curve(y_test, scores_test, sample_weight=np.abs(w_test))
        auc_val = roc_auc_score(y_test, scores_test, sample_weight=np.abs(w_test))
        print(f"\n  Test AUC = {auc_val:.4f}")

        # 7. Significance scan
        print("\nSignificance scan on test set...")
        scan = significance_scan(scores_sig, w_sig_te, scores_bkg, w_bkg_te, LUMI)
        print(f"  Best threshold: score > {scan['best_threshold']:.4f}")
        print(f"  S = {scan['best_S']:.4g},  B = {scan['best_B']:.4g},  Z_A = {scan['best_Z']:.4f}")

        # 8. Plots + save
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
