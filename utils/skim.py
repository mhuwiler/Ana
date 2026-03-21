"""
Snapshot/skim utilities for I/O optimization.

Writes slim ROOT files with only the branches used by the analysis,
eliminating the I/O bottleneck of reading 200+ branches over EOS.
Slim files are stored on /depot/ for fast local access.

Usage from cutflow_TrigEff.py:
    from utils.skim import run_skim, load_slim_or_eos

    # Generate slim files (one-time):
    python cutflow_TrigEff.py --skim --max-mc-files 10 --no-data

    # Force regeneration after adding new branches:
    python cutflow_TrigEff.py --reslim --max-mc-files 10 --no-data
"""

import os
import time

import ROOT

SLIM_DIR = "/depot/cms/users/das214/tmp/slim"


def get_used_branches(df, original_file):
    """Auto-detect tree branches actually used by the DataFrame chain.

    Compares columns in the Define/Filter chain against the original
    tree branches to find which file-level branches are referenced.
    """
    all_cols = set(str(c) for c in df.GetColumnNames())
    df_raw = ROOT.RDataFrame("Events", original_file)
    tree_branches = set(str(c) for c in df_raw.GetColumnNames())
    used = all_cols & tree_branches
    for extra in ("run", "luminosityBlock", "genWeight"):
        if extra in tree_branches:
            used.add(extra)
    return sorted(used)


def run_skim(mc_dict, mc_files, force=False):
    """Write slim ROOT files for all MC samples.

    Parameters
    ----------
    mc_dict : dict[str, RDataFrame]
        Sample name -> RDataFrame with all Define/Filter chains applied.
    mc_files : dict[str, list[str]]
        Sample name -> list of original EOS file paths.
    force : bool
        If True, overwrite existing slim files.
    """
    os.makedirs(SLIM_DIR, exist_ok=True)
    print(f"\n[skim] Output directory: {SLIM_DIR}")

    first_name = next(iter(mc_dict))
    branches = get_used_branches(mc_dict[first_name], mc_files[first_name][0])
    mc_branches = [b for b in branches if b not in ("run", "luminosityBlock")]
    print(f"[skim] Auto-detected {len(mc_branches)} branches")

    for name, df in mc_dict.items():
        out = os.path.join(SLIM_DIR, f"{name}.root")
        if os.path.exists(out) and not force:
            sz = os.path.getsize(out) / 1e6
            print(f"  {name:55s} exists ({sz:.1f} MB), skip")
            continue
        t0 = time.time()
        df.Snapshot("Events", out, mc_branches)
        sz = os.path.getsize(out) / 1e6
        print(f"  {name:55s} saved ({sz:.1f} MB, {time.time() - t0:.1f}s)")

    total = sum(
        os.path.getsize(os.path.join(SLIM_DIR, f)) / 1e6
        for f in os.listdir(SLIM_DIR)
        if f.endswith(".root")
    )
    print(f"\n[skim] Total: {total:.0f} MB in {SLIM_DIR}")


def load_slim_or_eos(name, eos_files, max_files=None):
    """Return RDataFrame from slim file if available, else from EOS.

    Parameters
    ----------
    name : str
        Sample name (used as slim filename).
    eos_files : list[str]
        Original EOS file paths.
    max_files : int or None
        Limit number of EOS files (ignored when slim exists).

    Returns
    -------
    ROOT.RDataFrame
    """
    slim = os.path.join(SLIM_DIR, f"{name}.root")
    if os.path.exists(slim):
        sz = os.path.getsize(slim) / 1e6
        print(f"  {name:55s} (slim, {sz:.1f} MB)")
        return ROOT.RDataFrame("Events", slim)
    files = eos_files[:max_files] if max_files else eos_files
    print(f"  {name:55s} ({len(files)} files)")
    return ROOT.RDataFrame("Events", files)
