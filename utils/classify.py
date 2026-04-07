"""BDT classifier for signal/background discrimination.

Provides:
  extract_features()      — RDataFrame → numpy arrays
  train_bdt()             — XGBoost training with GPU + progress bar
  significance_scan()     — Z_A vs BDT score threshold
  plot_roc/importance/score_dist/significance_scan/overtraining/loss_curves
  save_model()            — model + metadata
"""

import math
import os
import time

import numpy as np
import yaml

from analysis.cutflow import asimov_significance


# ═══════════════════════════════════════════════════════════════════════════════
#  Data extraction
# ═══════════════════════════════════════════════════════════════════════════════

def extract_features(mc_sel, feature_names, sig_modes=(20, 21, 22)):
    """Extract per-event feature arrays from RDataFrames."""
    extract_cols = list(set(feature_names)) + ["w", "decayMode"]

    sig_arrays, sig_weights = [], []
    bkg_arrays, bkg_weights = [], []

    for si, (sample_name, rdf) in enumerate(mc_sel.items()):
        print(f"  [{si+1}/{len(mc_sel)}] {sample_name[:55]}...", end="", flush=True)
        data = rdf.AsNumpy(extract_cols)
        print(f" {len(data['w'])} events")

        dm = data["decayMode"]
        w = data["w"]
        is_sig = np.isin(dm, sig_modes)
        is_bkg = ~is_sig
        matrix = np.column_stack([data[f].astype(np.float64) for f in feature_names])

        if np.any(is_sig):
            sig_arrays.append(matrix[is_sig])
            sig_weights.append(w[is_sig].astype(np.float64))
        if np.any(is_bkg):
            bkg_arrays.append(matrix[is_bkg])
            bkg_weights.append(w[is_bkg].astype(np.float64))

    X_sig = np.vstack(sig_arrays) if sig_arrays else np.zeros((0, len(feature_names)))
    X_bkg = np.vstack(bkg_arrays) if bkg_arrays else np.zeros((0, len(feature_names)))
    w_sig = np.concatenate(sig_weights) if sig_weights else np.zeros(0)
    w_bkg = np.concatenate(bkg_weights) if bkg_weights else np.zeros(0)
    return X_sig, w_sig, X_bkg, w_bkg


# ═══════════════════════════════════════════════════════════════════════════════
#  BDT Training
# ═══════════════════════════════════════════════════════════════════════════════

def train_bdt(X_train, y_train, w_train, X_val, y_val, w_val,
              feature_names, use_gpu=True, params=None):
    """Train XGBoost BDT classifier.

    Parameters
    ----------
    X_train, y_train, w_train : training data (features, labels, physics weights)
    X_val, y_val, w_val : validation data (for early stopping + loss curves)
    feature_names : list[str]
    use_gpu : bool

    Returns
    -------
    (model, eval_results) — trained model + dict of loss curves per iteration
    """
    import xgboost as xgb

    # Normalize weights: scale so signal total = background total
    # Physics weights are extreme (sig~0.0004, bkg~19M per event)
    # BDT needs balanced classes for stable training
    abs_w_tr = np.abs(w_train)
    abs_w_val = np.abs(w_val)
    sig_mask_tr = y_train == 1
    sig_mask_val = y_val == 1

    sum_sig_tr = np.sum(abs_w_tr[sig_mask_tr])
    sum_bkg_tr = np.sum(abs_w_tr[~sig_mask_tr])
    ratio = sum_bkg_tr / sum_sig_tr if sum_sig_tr > 0 else 1.0

    # Scale signal weights up to match background total
    norm_w_tr = abs_w_tr.copy()
    norm_w_tr[sig_mask_tr] *= ratio
    norm_w_val = abs_w_val.copy()
    norm_w_val[sig_mask_val] *= ratio

    print(f"  Weight normalization: sig scaled by {ratio:.2f}")
    print(f"    Train: {np.sum(norm_w_tr[sig_mask_tr]):.4g} sig, "
          f"{np.sum(norm_w_tr[~sig_mask_tr]):.4g} bkg")

    defaults = {
        "max_depth": 4,
        "learning_rate": 0.05,
        "n_estimators": 1000,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "min_child_weight": 10,
        "gamma": 1.0,
        "reg_alpha": 0.1,
        "reg_lambda": 1.0,
        "objective": "binary:logistic",
        "eval_metric": ["logloss", "auc"],
        "tree_method": "hist",
        "random_state": 42,
        "verbosity": 0,
        "early_stopping_rounds": 50,
    }
    if params:
        defaults.update(params)
    params = defaults

    if use_gpu:
        params["device"] = "cuda"

    # Progress bar callback
    try:
        from tqdm import tqdm
        from xgboost.callback import TrainingCallback

        class TqdmCallback(TrainingCallback):
            def __init__(self, total):
                self.pbar = tqdm(total=total, desc="  BDT Training", unit="tree",
                                bar_format="{l_bar}{bar:30}{r_bar}")
            def after_iteration(self, model, epoch, evals_log):
                self.pbar.update(1)
                if evals_log:
                    last_set = list(evals_log.keys())[-1]
                    val_loss = evals_log[last_set]["logloss"][-1]
                    val_auc = evals_log[last_set]["auc"][-1]
                    self.pbar.set_postfix({"val_loss": f"{val_loss:.4f}",
                                          "val_auc": f"{val_auc:.4f}"})
                return False
            def after_training(self, model):
                self.pbar.close()
                return model

        callbacks = [TqdmCallback(params["n_estimators"])]
    except ImportError:
        callbacks = None

    print(f"  Training XGBoost (GPU={use_gpu}, {len(X_train)} train, {len(X_val)} val)...")
    t0 = time.time()

    model = xgb.XGBClassifier(**params, callbacks=callbacks)
    model.fit(
        X_train, y_train,
        sample_weight=norm_w_tr,
        eval_set=[(X_train, y_train), (X_val, y_val)],
        sample_weight_eval_set=[norm_w_tr, norm_w_val],
        verbose=False,
    )
    eval_results = model.evals_result()

    dt = time.time() - t0
    best_iter = getattr(model, 'best_iteration', params["n_estimators"])
    print(f"  Training done in {dt:.1f}s (best iteration: {best_iter}/{params['n_estimators']})")

    return model, eval_results


# ═══════════════════════════════════════════════════════════════════════════════
#  Significance scan
# ═══════════════════════════════════════════════════════════════════════════════

def significance_scan(scores_sig, w_sig, scores_bkg, w_bkg, lumi,
                      n_points=200, min_bkg_events=5):
    """Scan BDT score threshold to find optimal significance on TEST set."""
    thresholds = np.linspace(0, 1, n_points)
    S_arr = np.zeros(n_points)
    B_arr = np.zeros(n_points)
    Z_arr = np.zeros(n_points)

    for i, t in enumerate(thresholds):
        sig_pass = scores_sig >= t
        bkg_pass = scores_bkg >= t
        n_bkg_unw = int(np.sum(bkg_pass))
        if n_bkg_unw < min_bkg_events:
            continue
        S = float(np.sum(w_sig[sig_pass])) * lumi
        B = float(np.sum(w_bkg[bkg_pass])) * lumi
        Z = asimov_significance(S, B)
        S_arr[i] = S
        B_arr[i] = B
        Z_arr[i] = Z if not math.isnan(Z) else 0.0

    best_idx = np.argmax(Z_arr)
    return {
        "thresholds": thresholds, "S_arr": S_arr, "B_arr": B_arr, "Z_arr": Z_arr,
        "best_idx": best_idx,
        "best_threshold": thresholds[best_idx],
        "best_S": S_arr[best_idx], "best_B": B_arr[best_idx], "best_Z": Z_arr[best_idx],
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  Multi-class: Data extraction
# ═══════════════════════════════════════════════════════════════════════════════

def extract_features_multiclass(mc_sel, feature_names, class_map):
    """Extract per-event feature arrays with integer class labels.

    Parameters
    ----------
    class_map : dict[int, int]
        Maps decayMode → class index. Events with decayMode not in class_map
        are skipped.

    Returns
    -------
    (X, y, w) — feature matrix, class labels, physics weights
    """
    extract_cols = list(set(feature_names)) + ["w", "decayMode"]

    all_X, all_y, all_w = [], [], []

    for si, (sample_name, rdf) in enumerate(mc_sel.items()):
        print(f"  [{si+1}/{len(mc_sel)}] {sample_name[:55]}...", end="", flush=True)
        data = rdf.AsNumpy(extract_cols)
        print(f" {len(data['w'])} events")

        dm = data["decayMode"]
        w = data["w"].astype(np.float64)
        matrix = np.column_stack([data[f].astype(np.float64) for f in feature_names])

        # Map decayMode to class index, skip unmapped events
        mask = np.isin(dm, list(class_map.keys()))
        if not np.any(mask):
            continue
        y = np.array([class_map[int(d)] for d in dm[mask]])
        all_X.append(matrix[mask])
        all_y.append(y)
        all_w.append(w[mask])

    X = np.vstack(all_X) if all_X else np.zeros((0, len(feature_names)))
    y = np.concatenate(all_y) if all_y else np.zeros(0, dtype=int)
    w = np.concatenate(all_w) if all_w else np.zeros(0)

    # Print per-class counts and remap to contiguous labels
    unique, counts = np.unique(y, return_counts=True)
    all_classes = sorted(set(class_map.values()))
    missing = set(all_classes) - set(unique)

    print(f"\n  Total: {len(X)} events across {len(unique)} classes")
    for c, n in zip(unique, counts):
        print(f"    class {int(c):2d}: {n:>10d} events")
    if missing:
        print(f"  WARNING: classes with 0 events: {sorted(missing)} — will be dropped")

    # Remap to contiguous 0..N-1 (XGBoost requires this)
    present = sorted(unique)
    remap = {old: new for new, old in enumerate(present)}
    y = np.array([remap[int(c)] for c in y])

    return X, y, w, present


# ═══════════════════════════════════════════════════════════════════════════════
#  Multi-class: BDT Training
# ═══════════════════════════════════════════════════════════════════════════════

def train_bdt_multiclass(X_train, y_train, w_train, X_val, y_val, w_val,
                         feature_names, n_classes, class_names=None,
                         use_gpu=True, params=None):
    """Train multi-class XGBoost BDT classifier.

    Weight normalization: each class scaled to equal total weight.

    Returns
    -------
    (model, eval_results)
    """
    import xgboost as xgb

    # Normalize weights: equalize total weight per class
    abs_w_tr = np.abs(w_train)
    abs_w_val = np.abs(w_val)
    norm_w_tr = abs_w_tr.copy()
    norm_w_val = abs_w_val.copy()

    class_sums = []
    for c in range(n_classes):
        mask = y_train == c
        s = np.sum(abs_w_tr[mask])
        class_sums.append(s)

    max_sum = max(s for s in class_sums if s > 0)
    print(f"  Weight normalization (per-class, target={max_sum:.4g}):")
    for c in range(n_classes):
        mask_tr = y_train == c
        mask_val = y_val == c
        if class_sums[c] > 0:
            ratio = max_sum / class_sums[c]
            norm_w_tr[mask_tr] *= ratio
            norm_w_val[mask_val] *= ratio
            name = class_names[c] if class_names else f"class {c}"
            n_tr = int(np.sum(mask_tr))
            print(f"    {name:20s}: {n_tr:>8d} events, scale={ratio:.2f}")
        else:
            name = class_names[c] if class_names else f"class {c}"
            print(f"    {name:20s}: 0 events (skipped)")

    defaults = {
        "max_depth": 4,
        "learning_rate": 0.05,
        "n_estimators": 1000,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "min_child_weight": 10,
        "gamma": 1.0,
        "reg_alpha": 0.1,
        "reg_lambda": 1.0,
        "objective": "multi:softprob",
        "num_class": n_classes,
        "eval_metric": "mlogloss",
        "tree_method": "hist",
        "random_state": 42,
        "verbosity": 0,
        "early_stopping_rounds": 50,
    }
    if params:
        defaults.update(params)
    params = defaults

    if use_gpu:
        params["device"] = "cuda"

    # Progress bar callback
    try:
        from tqdm import tqdm
        from xgboost.callback import TrainingCallback

        class TqdmCallback(TrainingCallback):
            def __init__(self, total):
                self.pbar = tqdm(total=total, desc="  BDT Training", unit="tree",
                                bar_format="{l_bar}{bar:30}{r_bar}")
            def after_iteration(self, model, epoch, evals_log):
                self.pbar.update(1)
                if evals_log:
                    last_set = list(evals_log.keys())[-1]
                    val_loss = evals_log[last_set]["mlogloss"][-1]
                    self.pbar.set_postfix({"val_mlogloss": f"{val_loss:.4f}"})
                return False
            def after_training(self, model):
                self.pbar.close()
                return model

        callbacks = [TqdmCallback(params["n_estimators"])]
    except ImportError:
        callbacks = None

    print(f"  Training XGBoost multi-class (GPU={use_gpu}, {n_classes} classes, "
          f"{len(X_train)} train, {len(X_val)} val)...")
    t0 = time.time()

    model = xgb.XGBClassifier(**params, callbacks=callbacks)
    model.fit(
        X_train, y_train,
        sample_weight=norm_w_tr,
        eval_set=[(X_train, y_train), (X_val, y_val)],
        sample_weight_eval_set=[norm_w_tr, norm_w_val],
        verbose=False,
    )
    eval_results = model.evals_result()

    dt = time.time() - t0
    best_iter = getattr(model, 'best_iteration', params["n_estimators"])
    print(f"  Training done in {dt:.1f}s (best iteration: {best_iter}/{params['n_estimators']})")

    return model, eval_results


# ═══════════════════════════════════════════════════════════════════════════════
#  Multi-class: Significance scan
# ═══════════════════════════════════════════════════════════════════════════════

def significance_scan_multiclass(proba, w, y, lumi, sig_classes,
                                 n_points=200, min_bkg_events=5):
    """Scan signal score threshold for multi-class BDT.

    Signal score = sum of probabilities for signal classes.
    Uses original physics weights for yield calculation.
    Includes statistical uncertainties on S, B, Z_A.
    """
    sig_score = sum(proba[:, c] for c in sig_classes)
    is_sig = np.isin(y, sig_classes)
    is_bkg = ~is_sig

    scores_sig = sig_score[is_sig]
    scores_bkg = sig_score[is_bkg]
    w_sig = w[is_sig]
    w_bkg = w[is_bkg]

    thresholds = np.linspace(0, 1, n_points)
    S_arr = np.zeros(n_points)
    B_arr = np.zeros(n_points)
    Z_arr = np.zeros(n_points)
    dS_arr = np.zeros(n_points)
    dB_arr = np.zeros(n_points)
    dZ_arr = np.zeros(n_points)

    for i, t in enumerate(thresholds):
        sig_pass = scores_sig >= t
        bkg_pass = scores_bkg >= t
        n_bkg_unw = int(np.sum(bkg_pass))
        if n_bkg_unw < min_bkg_events:
            continue
        S = float(np.sum(w_sig[sig_pass])) * lumi
        B = float(np.sum(w_bkg[bkg_pass])) * lumi
        dS = float(np.sqrt(np.sum(w_sig[sig_pass]**2))) * lumi
        dB = float(np.sqrt(np.sum(w_bkg[bkg_pass]**2))) * lumi
        Z = asimov_significance(S, B)
        S_arr[i] = S
        B_arr[i] = B
        Z_arr[i] = Z if not math.isnan(Z) else 0.0
        dS_arr[i] = dS
        dB_arr[i] = dB
        # Propagate uncertainty to Z_A
        if S > 0 and B > 0 and Z > 0:
            # dZ/dS and dZ/dB from Asimov formula
            dZdS = math.log(1 + S / B) / Z if Z > 0 else 0
            dZdB = (math.log(1 + S / B) - S / B) / Z if Z > 0 else 0
            dZ_arr[i] = math.sqrt((dZdS * dS)**2 + (dZdB * dB)**2)

    best_idx = np.argmax(Z_arr)
    return {
        "thresholds": thresholds, "S_arr": S_arr, "B_arr": B_arr, "Z_arr": Z_arr,
        "dS_arr": dS_arr, "dB_arr": dB_arr, "dZ_arr": dZ_arr,
        "best_idx": best_idx,
        "best_threshold": thresholds[best_idx],
        "best_S": S_arr[best_idx], "best_B": B_arr[best_idx], "best_Z": Z_arr[best_idx],
        "best_dS": dS_arr[best_idx], "best_dB": dB_arr[best_idx], "best_dZ": dZ_arr[best_idx],
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  Multi-class: Plots
# ═══════════════════════════════════════════════════════════════════════════════

def plot_confusion_matrix_multiclass(y_true, y_pred, class_names, outpath):
    """Confusion matrix using unweighted event counts, normalized per true class."""
    plt = _setup_plot()
    from sklearn.metrics import confusion_matrix
    import itertools

    cm = confusion_matrix(y_true, y_pred, labels=range(len(class_names)))
    # Normalize per row (per true class)
    row_sums = cm.sum(axis=1, keepdims=True)
    cm_norm = np.where(row_sums > 0, cm / row_sums, 0)

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(cm_norm, interpolation='nearest', cmap='Blues', vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for i, j in itertools.product(range(len(class_names)), range(len(class_names))):
        val = cm_norm[i, j]
        color = "white" if val > 0.5 else "black"
        ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=7)

    ax.set_xticks(range(len(class_names)))
    ax.set_yticks(range(len(class_names)))
    ax.set_xticklabels(class_names, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(class_names, fontsize=8)
    ax.set_xlabel("Predicted Class")
    ax.set_ylabel("True Class")
    ax.set_title("Confusion Matrix (normalized per true class)")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_roc_multiclass(y_test, proba, class_names, outpath):
    """One-vs-rest ROC for all classes on a single plot."""
    plt = _setup_plot()
    from sklearn.metrics import roc_curve, roc_auc_score

    fig, ax = plt.subplots(figsize=(8, 7))
    for c in range(len(class_names)):
        y_binary = (y_test == c).astype(int)
        if np.sum(y_binary) == 0 or np.sum(y_binary) == len(y_binary):
            continue
        fpr, tpr, _ = roc_curve(y_binary, proba[:, c])
        auc_val = roc_auc_score(y_binary, proba[:, c])
        ax.plot(fpr, tpr, linewidth=1.5, label=f"{class_names[c]} (AUC={auc_val:.3f})")

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curves (one-vs-rest)")
    ax.legend(loc="lower right", fontsize=7)
    ax.grid(True, alpha=0.3)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_score_distribution_multiclass(sig_score, w, y, class_names, sig_classes, outpath):
    """Histogram of signal score (sum of signal probs), stacked by true class."""
    plt = _setup_plot()
    fig, ax = plt.subplots(figsize=(8, 6))
    bins = np.linspace(0, 1, 50)

    # Background classes first (stacked), then signal overlaid
    bkg_classes = [c for c in range(len(class_names)) if c not in sig_classes]
    bkg_data = [sig_score[y == c] for c in bkg_classes]
    bkg_weights = [np.abs(w[y == c]) for c in bkg_classes]
    bkg_labels = [class_names[c] for c in bkg_classes]

    if any(len(d) > 0 for d in bkg_data):
        ax.hist(bkg_data, bins=bins, weights=bkg_weights, density=True,
                histtype="stepfilled", alpha=0.5, stacked=True, label=bkg_labels)

    for c in sig_classes:
        mask = y == c
        if np.sum(mask) > 0:
            ax.hist(sig_score[mask], bins=bins, weights=np.abs(w[mask]), density=True,
                    histtype="step", linewidth=2, label=class_names[c])

    ax.set_xlabel("Signal Score (sum of signal class probabilities)")
    ax.set_ylabel("Normalized")
    ax.set_title("Multi-class BDT Signal Score")
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_yscale("log")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_overtraining_multiclass(model, X_train, y_train, X_test, y_test,
                                 sig_classes, outpath):
    """Overtraining check using signal score for multi-class BDT."""
    plt = _setup_plot()
    from scipy.stats import ks_2samp

    proba_train = model.predict_proba(X_train)
    proba_test = model.predict_proba(X_test)
    sig_score_train = sum(proba_train[:, c] for c in sig_classes)
    sig_score_test = sum(proba_test[:, c] for c in sig_classes)

    is_sig_train = np.isin(y_train, sig_classes)
    is_sig_test = np.isin(y_test, sig_classes)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    bins = np.linspace(0, 1, 50)
    for ax, label, mask_tr, mask_te in [
        (ax1, "Signal", is_sig_train, is_sig_test),
        (ax2, "Background", ~is_sig_train, ~is_sig_test),
    ]:
        tr = sig_score_train[mask_tr]
        te = sig_score_test[mask_te]
        if len(tr) == 0 or len(te) == 0:
            continue
        ks_stat, ks_pval = ks_2samp(tr, te)
        ax.hist(tr, bins=bins, density=True, histtype="stepfilled", alpha=0.4,
                color="tab:blue", label="Train")
        ax.hist(te, bins=bins, density=True, histtype="step", linewidth=2,
                color="tab:red", label="Test")
        ax.set_xlabel("Signal Score"); ax.set_ylabel("Normalized")
        ax.set_title(f"{label} (KS p={ks_pval:.3f})")
        ax.legend(); ax.grid(True, alpha=0.3)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Plots (binary + shared)
# ═══════════════════════════════════════════════════════════════════════════════

def _setup_plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_roc(fpr, tpr, auc_val, outpath):
    plt = _setup_plot()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(fpr, tpr, linewidth=2, label=f"BDT (AUC = {auc_val:.4f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate (Background efficiency)")
    ax.set_ylabel("True Positive Rate (Signal efficiency)")
    ax.set_title("ROC Curve")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_feature_importance(model, feature_names, outpath, top_n=23):
    plt = _setup_plot()
    importances = model.feature_importances_
    indices = np.argsort(importances)[-top_n:]
    fig, ax = plt.subplots(figsize=(8, max(4, len(indices) * 0.35)))
    ax.barh(range(len(indices)), importances[indices], align="center")
    ax.set_yticks(range(len(indices)))
    ax.set_yticklabels([feature_names[i] for i in indices])
    ax.set_xlabel("Feature Importance (gain)")
    ax.set_title("BDT Feature Importance")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_score_distribution(scores_sig, w_sig, scores_bkg, w_bkg, outpath):
    plt = _setup_plot()
    fig, ax = plt.subplots(figsize=(8, 6))
    bins = np.linspace(0, 1, 50)
    ax.hist(scores_bkg, bins=bins, weights=np.abs(w_bkg), density=True,
            histtype="stepfilled", alpha=0.4, color="tab:blue", label="Background")
    ax.hist(scores_sig, bins=bins, weights=np.abs(w_sig), density=True,
            histtype="step", linewidth=2, color="tab:red", label="Signal")
    ax.set_xlabel("BDT Score"); ax.set_ylabel("Normalized")
    ax.set_title("BDT Score Distribution")
    ax.legend(); ax.grid(True, alpha=0.3); ax.set_yscale("log")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_significance_scan(scan, outpath):
    plt = _setup_plot()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                     gridspec_kw={"height_ratios": [2, 1], "hspace": 0.05})
    t = scan["thresholds"]
    ax1.plot(t, scan["Z_arr"], linewidth=2, color="tab:red")
    ax1.axvline(scan["best_threshold"], color="gray", linestyle="--", alpha=0.7,
                label=f'Best: score > {scan["best_threshold"]:.3f}, Z_A = {scan["best_Z"]:.4f}')
    ax1.set_ylabel(r"$Z_A$ (Asimov significance)")
    ax1.legend(); ax1.grid(True, alpha=0.3); ax1.tick_params(labelbottom=False)
    ax2.plot(t, scan["S_arr"], linewidth=1.5, color="tab:red", label="S")
    ax2.plot(t, scan["B_arr"], linewidth=1.5, color="tab:blue", label="B")
    ax2.set_xlabel("BDT Score Threshold"); ax2.set_ylabel("Weighted Events")
    ax2.set_yscale("log"); ax2.legend(); ax2.grid(True, alpha=0.3)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_overtraining(model, X_train, y_train, X_test, y_test, outpath):
    plt = _setup_plot()
    from scipy.stats import ks_2samp
    scores_train = model.predict_proba(X_train)[:, 1]
    scores_test = model.predict_proba(X_test)[:, 1]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    bins = np.linspace(0, 1, 50)
    for ax, label, mv in [(ax1, "Signal", 1), (ax2, "Background", 0)]:
        tr = scores_train[y_train == mv]
        te = scores_test[y_test == mv]
        ks_stat, ks_pval = ks_2samp(tr, te)
        ax.hist(tr, bins=bins, density=True, histtype="stepfilled", alpha=0.4,
                color="tab:blue", label="Train")
        ax.hist(te, bins=bins, density=True, histtype="step", linewidth=2,
                color="tab:red", label="Test")
        ax.set_xlabel("BDT Score"); ax.set_ylabel("Normalized")
        ax.set_title(f"{label} (KS p={ks_pval:.3f})")
        ax.legend(); ax.grid(True, alpha=0.3)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_loss_curves(eval_results, best_iteration, outpath):
    """Plot train and validation loss curves vs boosting round.

    Handles both binary (logloss + auc) and multi-class (mlogloss only).
    """
    plt = _setup_plot()
    keys = list(eval_results.keys())  # ['validation_0', 'validation_1']
    train_key, val_key = keys[0], keys[1]

    # Detect metric name
    train_metrics = eval_results[train_key]
    loss_key = "mlogloss" if "mlogloss" in train_metrics else "logloss"
    has_auc = "auc" in train_metrics

    if has_auc:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 5))

    # Loss curve
    ax1.plot(train_metrics[loss_key], linewidth=1.5,
             color="tab:blue", label="Train", alpha=0.7)
    ax1.plot(eval_results[val_key][loss_key], linewidth=1.5,
             color="tab:red", label="Validation")
    if best_iteration:
        ax1.axvline(best_iteration, color="gray", linestyle="--", alpha=0.5,
                    label=f"Best iter: {best_iteration}")
    loss_label = "Multi-class Log Loss" if loss_key == "mlogloss" else "Log Loss"
    ax1.set_xlabel("Boosting Round"); ax1.set_ylabel(loss_label)
    ax1.set_title("Training Loss Curve"); ax1.legend(); ax1.grid(True, alpha=0.3)

    # AUC (binary only)
    if has_auc:
        ax2.plot(train_metrics["auc"], linewidth=1.5,
                 color="tab:blue", label="Train", alpha=0.7)
        ax2.plot(eval_results[val_key]["auc"], linewidth=1.5,
                 color="tab:red", label="Validation")
        if best_iteration:
            ax2.axvline(best_iteration, color="gray", linestyle="--", alpha=0.5,
                        label=f"Best iter: {best_iteration}")
        ax2.set_xlabel("Boosting Round"); ax2.set_ylabel("AUC")
        ax2.set_title("AUC Curve"); ax2.legend(); ax2.grid(True, alpha=0.3)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)
    print(f"  Saved: {outpath}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Save model + results
# ═══════════════════════════════════════════════════════════════════════════════

def save_model(model, feature_names, scan, auc_val, outdir):
    """Save BDT model and metadata."""
    models_dir = os.path.join(outdir, "models")
    os.makedirs(models_dir, exist_ok=True)

    model_path = os.path.join(models_dir, "bdt_model.json")
    model.save_model(model_path)
    print(f"  Saved: {model_path}")

    meta = {
        "features": feature_names,
        "auc": float(auc_val) if auc_val is not None else None,
        "best_threshold": float(scan["best_threshold"]),
        "best_S": float(scan["best_S"]),
        "best_B": float(scan["best_B"]),
        "best_Z_A": float(scan["best_Z"]),
        "n_estimators": int(model.n_estimators),
        "best_iteration": int(getattr(model, 'best_iteration', model.n_estimators)),
        "max_depth": int(model.max_depth),
    }
    meta_path = os.path.join(models_dir, "bdt_meta.yaml")
    with open(meta_path, "w") as f:
        yaml.dump(meta, f, default_flow_style=False)
    print(f"  Saved: {meta_path}")

    summary_path = os.path.join(outdir, "summary.txt")
    with open(summary_path, "w") as f:
        f.write("BDT Classifier Summary\n" + "=" * 50 + "\n\n")
        f.write(f"AUC: {auc_val:.4f}\n")
        f.write(f"Best BDT score threshold: {scan['best_threshold']:.4f}\n")
        f.write(f"S = {scan['best_S']:.4g}\n")
        f.write(f"B = {scan['best_B']:.4g}\n")
        f.write(f"Z_A = {scan['best_Z']:.4f}\n\n")
        f.write(f"Features ({len(feature_names)}):\n")
        for fn in feature_names:
            f.write(f"  {fn}\n")
    print(f"  Saved: {summary_path}")
