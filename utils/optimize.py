"""Cut optimization via RGS (Random Grid Search) and Grid Search.

Provides methods:
  run_rgs()            — signal events as candidate thresholds
  run_rgs_parallel()   — same, parallelized (CPU multiprocess or JIT/GPU)
  run_grid()           — evenly spaced grid
  run_grid_parallel()  — same, parallelized

Backends (auto-detected):
  "gpu"   — numba CUDA kernel on GPU (fastest, ~2-5s for 50k candidates)
  "jit"   — numba JIT with prange on CPU (fast, ~5-10s)
  "numpy" — pure numpy with ProcessPoolExecutor fallback (~30-60s)
"""

import math
import os
import time

import numpy as np

from analysis.cutflow import asimov_significance


# ═══════════════════════════════════════════════════════════════════════════════
#  Backend detection
# ═══════════════════════════════════════════════════════════════════════════════

def _get_backend():
    """Auto-detect best available compute backend."""
    try:
        from numba import cuda
        if cuda.is_available():
            return "gpu"
    except (ImportError, Exception):
        pass
    try:
        import numba
        return "jit"
    except ImportError:
        pass
    return "numpy"


# ═══════════════════════════════════════════════════════════════════════════════
#  Data extraction
# ═══════════════════════════════════════════════════════════════════════════════

def extract_arrays(mc_sel, opt_vars, sig_modes=(20, 21, 22)):
    """Extract per-event numpy arrays from RDataFrames."""
    col_names = list(set(v[0] for v in opt_vars))
    extract_cols = col_names + ["w", "decayMode"]

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
        matrix = np.column_stack([data[v[0]].astype(np.float64) for v in opt_vars])

        if np.any(is_sig):
            sig_arrays.append(matrix[is_sig])
            sig_weights.append(w[is_sig].astype(np.float64))
        if np.any(is_bkg):
            bkg_arrays.append(matrix[is_bkg])
            bkg_weights.append(w[is_bkg].astype(np.float64))

    sig = np.vstack(sig_arrays) if sig_arrays else np.zeros((0, len(opt_vars)))
    bkg = np.vstack(bkg_arrays) if bkg_arrays else np.zeros((0, len(opt_vars)))
    sig_w = np.concatenate(sig_weights) if sig_weights else np.zeros(0)
    bkg_w = np.concatenate(bkg_weights) if bkg_weights else np.zeros(0)
    return sig, sig_w, bkg, bkg_w


def train_test_split(sig, sig_w, bkg, bkg_w, test_fraction=0.5, seed=42):
    """Split signal and background arrays into train/test sets.

    Parameters
    ----------
    test_fraction : float
        Fraction of events to use for testing (default 0.5 = 50/50 split).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    (sig_train, sig_w_train, bkg_train, bkg_w_train,
     sig_test, sig_w_test, bkg_test, bkg_w_test)
    """
    rng = np.random.RandomState(seed)

    n_sig = len(sig)
    n_bkg = len(bkg)
    sig_idx = rng.permutation(n_sig)
    bkg_idx = rng.permutation(n_bkg)

    n_sig_test = int(n_sig * test_fraction)
    n_bkg_test = int(n_bkg * test_fraction)

    sig_train = sig[sig_idx[n_sig_test:]]
    sig_w_train = sig_w[sig_idx[n_sig_test:]]
    sig_test = sig[sig_idx[:n_sig_test]]
    sig_w_test = sig_w[sig_idx[:n_sig_test]]

    bkg_train = bkg[bkg_idx[n_bkg_test:]]
    bkg_w_train = bkg_w[bkg_idx[n_bkg_test:]]
    bkg_test = bkg[bkg_idx[:n_bkg_test]]
    bkg_w_test = bkg_w[bkg_idx[:n_bkg_test]]

    print(f"  Train: {len(sig_train)} sig, {len(bkg_train)} bkg")
    print(f"  Test:  {len(sig_test)} sig, {len(bkg_test)} bkg")

    return (sig_train, sig_w_train, bkg_train, bkg_w_train,
            sig_test, sig_w_test, bkg_test, bkg_w_test)


# ═══════════════════════════════════════════════════════════════════════════════
#  Backend: numpy (fallback)
# ═══════════════════════════════════════════════════════════════════════════════

def _eval_numpy(candidates, sig, bkg, sig_w, bkg_w, lumi,
                show_progress=False, min_bkg_events=5):
    """Pure numpy evaluation — loop over candidates, vectorized over events."""
    K = len(candidates)
    results = np.zeros((K, 3))

    iterator = range(K)
    if show_progress:
        try:
            from tqdm import tqdm
            iterator = tqdm(iterator, desc="  numpy", unit="pt",
                           bar_format="{l_bar}{bar:30}{r_bar}")
        except ImportError:
            pass

    for i in iterator:
        thresh = candidates[i]
        sig_pass = np.all(sig >= thresh, axis=1)
        bkg_pass = np.all(bkg >= thresh, axis=1)
        n_bkg_unw = int(np.sum(bkg_pass))
        if n_bkg_unw < min_bkg_events:
            continue
        S = float(np.sum(sig_w[sig_pass])) * lumi
        B = float(np.sum(bkg_w[bkg_pass])) * lumi
        Z = asimov_significance(S, B)
        results[i] = [S, B, Z if not math.isnan(Z) else 0.0]
    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Backend: numba JIT + prange (CPU)
# ═══════════════════════════════════════════════════════════════════════════════

def _build_jit_functions():
    """Build JIT-compiled evaluation function (lazy, only when needed)."""
    import numba

    @numba.njit
    def _asimov_jit(S, B):
        if B <= 0.0 or S <= 0.0:
            return 0.0
        val = 2.0 * ((S + B) * math.log(1.0 + S / B) - S)
        if val < 0.0:
            return 0.0
        return math.sqrt(val)

    @numba.njit(parallel=True)
    def _eval_jit_inner(candidates, sig, bkg, sig_w, bkg_w, lumi, min_bkg):
        K = candidates.shape[0]
        D = candidates.shape[1]
        results = np.zeros((K, 3))
        for i in numba.prange(K):
            S = 0.0
            for j in range(sig.shape[0]):
                passes = True
                for d in range(D):
                    if sig[j, d] < candidates[i, d]:
                        passes = False
                        break
                if passes:
                    S += sig_w[j]
            B = 0.0
            n_bkg = 0
            for j in range(bkg.shape[0]):
                passes = True
                for d in range(D):
                    if bkg[j, d] < candidates[i, d]:
                        passes = False
                        break
                if passes:
                    B += bkg_w[j]
                    n_bkg += 1
            S *= lumi
            B *= lumi
            if n_bkg >= min_bkg:
                results[i, 0] = S
                results[i, 1] = B
                results[i, 2] = _asimov_jit(S, B)
        return results

    return _eval_jit_inner


_jit_func = None


def _eval_jit(candidates, sig, bkg, sig_w, bkg_w, lumi, min_bkg_events=5):
    """numba JIT evaluation with prange parallelism over candidates."""
    global _jit_func
    if _jit_func is None:
        print("  [JIT] Compiling (first call only)...", flush=True)
        _jit_func = _build_jit_functions()
        # Warm up with tiny arrays
        _jit_func(candidates[:2].copy(), sig[:10].copy(), bkg[:10].copy(),
                  sig_w[:10].copy(), bkg_w[:10].copy(), lumi, min_bkg_events)
        print("  [JIT] Compiled.")
    return _jit_func(candidates, sig, bkg, sig_w, bkg_w, lumi, min_bkg_events)


# ═══════════════════════════════════════════════════════════════════════════════
#  Backend: numba CUDA (GPU)
# ═══════════════════════════════════════════════════════════════════════════════

def _build_cuda_kernel():
    """Build CUDA kernel (lazy, only when needed)."""
    from numba import cuda

    @cuda.jit
    def _eval_cuda_kernel(candidates, sig, bkg, sig_w, bkg_w, results, lumi):
        i = cuda.grid(1)
        if i >= candidates.shape[0]:
            return
        D = candidates.shape[1]
        S = 0.0
        for j in range(sig.shape[0]):
            passes = True
            for d in range(D):
                if sig[j, d] < candidates[i, d]:
                    passes = False
                    break
            if passes:
                S += sig_w[j]
        B = 0.0
        n_bkg = 0
        for j in range(bkg.shape[0]):
            passes = True
            for d in range(D):
                if bkg[j, d] < candidates[i, d]:
                    passes = False
                    break
            if passes:
                B += bkg_w[j]
                n_bkg += 1
        results[i, 0] = S * lumi
        results[i, 1] = B * lumi
        results[i, 2] = n_bkg  # Z_A computed on CPU after

    return _eval_cuda_kernel


_cuda_kernel = None


def _eval_cuda(candidates, sig, bkg, sig_w, bkg_w, lumi, min_bkg_events=5):
    """CUDA GPU evaluation — one thread per candidate."""
    from numba import cuda

    global _cuda_kernel
    if _cuda_kernel is None:
        print("  [CUDA] Compiling kernel...", flush=True)
        _cuda_kernel = _build_cuda_kernel()

    K = len(candidates)

    # Copy to GPU
    d_candidates = cuda.to_device(candidates)
    d_sig = cuda.to_device(sig)
    d_bkg = cuda.to_device(bkg)
    d_sig_w = cuda.to_device(sig_w)
    d_bkg_w = cuda.to_device(bkg_w)
    d_results = cuda.device_array((K, 3), dtype=np.float64)

    # Launch kernel
    threads_per_block = 256
    blocks = (K + threads_per_block - 1) // threads_per_block
    _cuda_kernel[blocks, threads_per_block](
        d_candidates, d_sig, d_bkg, d_sig_w, d_bkg_w, d_results, lumi)

    # Copy results back
    results = d_results.copy_to_host()

    # Compute Z_A on CPU (CUDA can't call Python math functions)
    for i in range(K):
        S, B, n_bkg = results[i]
        if int(n_bkg) < min_bkg_events:
            results[i] = [0.0, 0.0, 0.0]
        else:
            Z = asimov_significance(S, B)
            results[i, 2] = Z if not math.isnan(Z) else 0.0

    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Dispatch: select backend
# ═══════════════════════════════════════════════════════════════════════════════

def _evaluate_candidates(candidates, sig, bkg, sig_w, bkg_w, lumi,
                         show_progress=False, min_bkg_events=5, backend="auto"):
    """Evaluate Z_A for each candidate — dispatches to best available backend."""
    if backend == "auto":
        backend = _get_backend()

    if backend == "gpu":
        return _eval_cuda(candidates, sig, bkg, sig_w, bkg_w, lumi, min_bkg_events)
    elif backend == "jit":
        return _eval_jit(candidates, sig, bkg, sig_w, bkg_w, lumi, min_bkg_events)
    else:
        return _eval_numpy(candidates, sig, bkg, sig_w, bkg_w, lumi,
                           show_progress, min_bkg_events)


# ═══════════════════════════════════════════════════════════════════════════════
#  Methods
# ═══════════════════════════════════════════════════════════════════════════════

def run_rgs(sig, sig_w, bkg, bkg_w, opt_vars, lumi,
            n_samples=5000, seed=42, backend="auto"):
    """Random Grid Search: use signal events as candidate thresholds."""
    rng = np.random.RandomState(seed)
    n = min(n_samples, len(sig))
    idx = rng.choice(len(sig), n, replace=False)
    candidates = sig[idx]

    be = backend if backend != "auto" else _get_backend()
    print(f"[RGS] {len(candidates)} candidates, {len(sig)} sig, "
          f"{len(bkg)} bkg, backend={be}")
    t0 = time.time()
    results = _evaluate_candidates(candidates, sig, bkg, sig_w, bkg_w, lumi,
                                    show_progress=(be == "numpy"), backend=backend)
    dt = time.time() - t0
    print(f"[RGS] Done in {dt:.1f}s")
    return _pack_result(candidates, results, sig_w, lumi, f"RGS ({be})")


def run_rgs_parallel(sig, sig_w, bkg, bkg_w, opt_vars, lumi,
                     n_samples=5000, seed=42, n_workers=8, backend="auto"):
    """RGS — parallel. Uses JIT/GPU internally if available, else ProcessPoolExecutor."""
    rng = np.random.RandomState(seed)
    n = min(n_samples, len(sig))
    idx = rng.choice(len(sig), n, replace=False)
    candidates = sig[idx]

    be = backend if backend != "auto" else _get_backend()
    print(f"[RGS Parallel] {len(candidates)} candidates, backend={be}")
    t0 = time.time()

    if be in ("gpu", "jit"):
        # JIT/GPU handle parallelism internally
        results = _evaluate_candidates(candidates, sig, bkg, sig_w, bkg_w, lumi,
                                        backend=be)
    else:
        # numpy fallback: use ProcessPoolExecutor
        results = _run_parallel_numpy(candidates, sig, bkg, sig_w, bkg_w, lumi, n_workers)

    dt = time.time() - t0
    print(f"[RGS Parallel] Done in {dt:.1f}s")
    return _pack_result(candidates, results, sig_w, lumi, f"RGS Parallel ({be})")


def run_grid(sig, sig_w, bkg, bkg_w, opt_vars, lumi, n_steps=10, backend="auto"):
    """Grid search: evaluate evenly spaced grid."""
    candidates = _build_grid(opt_vars, n_steps)
    be = backend if backend != "auto" else _get_backend()
    print(f"[Grid] {n_steps} steps/var, {len(candidates)} points, backend={be}")
    t0 = time.time()
    results = _evaluate_candidates(candidates, sig, bkg, sig_w, bkg_w, lumi,
                                    show_progress=(be == "numpy"), backend=backend)
    dt = time.time() - t0
    print(f"[Grid] Done in {dt:.1f}s")
    return _pack_result(candidates, results, sig_w, lumi, f"Grid ({be})")


def run_grid_parallel(sig, sig_w, bkg, bkg_w, opt_vars, lumi,
                      n_steps=10, n_workers=8, backend="auto"):
    """Grid search — parallel."""
    candidates = _build_grid(opt_vars, n_steps)
    be = backend if backend != "auto" else _get_backend()
    print(f"[Grid Parallel] {len(candidates)} points, backend={be}")
    t0 = time.time()

    if be in ("gpu", "jit"):
        results = _evaluate_candidates(candidates, sig, bkg, sig_w, bkg_w, lumi,
                                        backend=be)
    else:
        results = _run_parallel_numpy(candidates, sig, bkg, sig_w, bkg_w, lumi, n_workers)

    dt = time.time() - t0
    print(f"[Grid Parallel] Done in {dt:.1f}s")
    return _pack_result(candidates, results, sig_w, lumi, f"Grid Parallel ({be})")


# ═══════════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _build_grid(opt_vars, n_steps):
    """Build evenly spaced grid from opt_vars scan ranges."""
    grids = [np.linspace(lo, hi, n_steps) for _, _, (lo, hi) in opt_vars]
    mesh = np.meshgrid(*grids, indexing="ij")
    return np.column_stack([m.ravel() for m in mesh])


def _run_parallel_numpy(candidates, sig, bkg, sig_w, bkg_w, lumi, n_workers):
    """numpy fallback: split candidates across ProcessPoolExecutor."""
    from concurrent.futures import ProcessPoolExecutor, as_completed

    n_chunks = min(n_workers * 4, len(candidates))
    chunks = np.array_split(candidates, n_chunks)

    results_list = [None] * n_chunks
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        future_to_idx = {
            pool.submit(_eval_numpy, c, sig, bkg, sig_w, bkg_w, lumi): i
            for i, c in enumerate(chunks)
        }
        try:
            from tqdm import tqdm
            with tqdm(total=n_chunks, desc="  Parallel", unit="chunk",
                     bar_format="{l_bar}{bar:30}{r_bar}") as pbar:
                for future in as_completed(future_to_idx):
                    idx = future_to_idx[future]
                    results_list[idx] = future.result()
                    pbar.update(1)
        except ImportError:
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                results_list[idx] = future.result()

    return np.vstack(results_list)


def _pack_result(candidates, results, sig_w, lumi, method):
    """Pack optimization results into a dict."""
    best_idx = np.argmax(results[:, 2])
    S, B, Z = results[best_idx]
    total_sig = float(np.sum(sig_w)) * lumi
    sig_eff = S / total_sig * 100 if total_sig > 0 else 0
    return {
        "thresholds": candidates[best_idx],
        "S": S, "B": B, "Z_A": Z,
        "sig_eff": sig_eff, "all_results": results, "method": method,
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  Output
# ═══════════════════════════════════════════════════════════════════════════════

def evaluate_on_test(result, sig_test, sig_w_test, bkg_test, bkg_w_test, lumi):
    """Evaluate the optimized thresholds on the test set.

    Updates result dict in-place with test_S, test_B, test_Z_A, test_sig_eff.
    """
    thresh = result["thresholds"]
    sig_pass = np.all(sig_test >= thresh, axis=1)
    bkg_pass = np.all(bkg_test >= thresh, axis=1)
    S = float(np.sum(sig_w_test[sig_pass])) * lumi
    B = float(np.sum(bkg_w_test[bkg_pass])) * lumi
    Z = asimov_significance(S, B)
    total_sig = float(np.sum(sig_w_test)) * lumi
    sig_eff = S / total_sig * 100 if total_sig > 0 else 0
    result["test_S"] = S
    result["test_B"] = B
    result["test_Z_A"] = Z if not math.isnan(Z) else 0.0
    result["test_sig_eff"] = sig_eff
    result["test_n_bkg"] = int(np.sum(bkg_pass))
    return result


def print_results(result, opt_vars):
    """Print formatted table of optimized cuts with train/test performance."""
    print(f"\n{'=' * 70}")
    print(f"  {result['method']}")
    print(f"{'=' * 70}")
    print(f"\n{'Variable':<25s} {'Dir':>3s} {'Threshold':>12s}")
    print("-" * 45)
    for d, (col, direction, _) in enumerate(opt_vars):
        print(f"{col:<25s} {direction:>3s} {result['thresholds'][d]:>12.4f}")

    print(f"\n  {'':15s} {'S':>10s} {'B':>10s} {'Z_A':>10s} {'Eff(%)':>10s}")
    print(f"  {'Train:':15s} {result['S']:>10.4g} {result['B']:>10.4g} "
          f"{result['Z_A']:>10.4f} {result['sig_eff']:>10.1f}")
    if "test_S" in result:
        print(f"  {'Test:':15s} {result['test_S']:>10.4g} {result['test_B']:>10.4g} "
              f"{result['test_Z_A']:>10.4f} {result['test_sig_eff']:>10.1f}")
    if result['B'] > 0 and result['B'] < 10:
        print(f"  WARNING: Train B = {result['B']:.2f} — fewer than 10 weighted bkg events")
    if "test_B" in result and result['test_B'] > 0 and result['test_B'] < 10:
        print(f"  WARNING: Test B = {result['test_B']:.2f} — fewer than 10 weighted bkg events")


def save_results(results_dict, opt_vars, out_dir):
    """Save optimization results to files."""
    import yaml
    os.makedirs(out_dir, exist_ok=True)

    yaml_data = {}
    for method_name, result in results_dict.items():
        cuts = {}
        for d, (col, direction, _) in enumerate(opt_vars):
            key = col
            val = float(result['thresholds'][d])
            cuts[key] = f"{col} > {val:.4f}" if direction == ">" else f"{col} < {val:.4f}"
        yaml_data[method_name] = {
            "cuts": cuts, "S": float(result['S']), "B": float(result['B']),
            "Z_A": float(result['Z_A']), "signal_eff_pct": float(result['sig_eff']),
        }

    yaml_path = os.path.join(out_dir, "results.yaml")
    with open(yaml_path, "w") as f:
        f.write("# Cut optimization results\n# Copy 'cuts' into config/cuts.yaml\n\n")
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
    print(f"  Saved: {yaml_path}")

    summary_path = os.path.join(out_dir, "summary.txt")
    with open(summary_path, "w") as f:
        f.write("Cut Optimization Summary\n" + "=" * 70 + "\n\n")
        for method_name, result in results_dict.items():
            f.write(f"{method_name}\n" + "-" * 45 + "\n")
            f.write(f"{'Variable':<25s} {'Dir':>3s} {'Threshold':>12s}\n")
            for d, (col, direction, _) in enumerate(opt_vars):
                f.write(f"{col:<25s} {direction:>3s} {result['thresholds'][d]:>12.4f}\n")
            f.write(f"\nS = {result['S']:.4g},  B = {result['B']:.4g},  "
                    f"Z_A = {result['Z_A']:.4f},  Signal eff = {result['sig_eff']:.1f}%\n\n")
    print(f"  Saved: {summary_path}")

    for method_name, result in results_dict.items():
        tag = method_name.lower().replace(" ", "_").replace("(", "").replace(")", "")
        npz_path = os.path.join(out_dir, f"{tag}_scan.npz")
        np.savez_compressed(npz_path, all_results=result['all_results'],
                            best_thresholds=result['thresholds'],
                            opt_var_names=[v[0] for v in opt_vars],
                            opt_var_dirs=[v[1] for v in opt_vars])
        print(f"  Saved: {npz_path}")


def plot_fom_ranking(result, outpath):
    """Plot Z_A vs rank (sorted descending)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z_vals = result["all_results"][:, 2]
    z_sorted = np.sort(z_vals[z_vals > 0])[::-1]
    if len(z_sorted) == 0:
        print(f"  [plot] No valid candidates — skipping {outpath}")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(len(z_sorted)), z_sorted, linewidth=1.5)
    ax.set_xlabel("Candidate rank")
    ax.set_ylabel(r"$Z_A$ (Asimov significance)")
    ax.set_title(f"{result['method']} — FOM ranking")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, min(len(z_sorted), 500))

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Saved: {outpath}")
