"""
Snapshot/skim utilities for I/O optimization.

Writes slim ROOT files with only the branches used by the analysis.
Supports incremental branch addition via TTree::AddFriend — new branches
are stored in a separate _friends.root file without rebuilding the main slim.

Slim files are named with file count: {sample}_{N}files.root
so different scripts with different MAX_MC_FILES don't conflict.

Usage:
    # Explicit skim (all branches, all files):
    python skim.py

    # Auto-skim happens automatically in load_and_run() on first run.
"""

import json
import os
import re
import time

import ROOT

SLIM_DIR = "/depot/cms/users/das214/tmp/slim"

# Collects all C++ expression strings from r_define/r_filter calls
rdf_exprs = []


# ── Expression tracking wrappers ─────────────────────────────────────────────

def r_define(df, col, expr):
    """Define a column and record the expression for branch auto-detection."""
    rdf_exprs.append(expr)
    return df.Define(col, expr)


def r_filter(df, expr):
    """Filter rows and record the expression for branch auto-detection."""
    rdf_exprs.append(expr)
    return df.Filter(expr)


# ── Slim file naming (with file count) ──────────────────────────────────────

def _slim_stem(name, n_files):
    """Build the base name for slim files: {sample}_{N}files."""
    return f"{name}_{n_files}files"


def _slim_path(name, n_files):
    return os.path.join(SLIM_DIR, f"{_slim_stem(name, n_files)}.root")


def _friends_path(name, n_files):
    return os.path.join(SLIM_DIR, f"{_slim_stem(name, n_files)}_friends.root")


def _branches_file(name, n_files):
    return os.path.join(SLIM_DIR, f"{_slim_stem(name, n_files)}.branches.json")


# ── Branch tracking (JSON sidecar files) ─────────────────────────────────────

def _load_known_branches(name, n_files):
    path = _branches_file(name, n_files)
    if os.path.exists(path):
        with open(path) as f:
            return set(json.load(f))
    return set()


def _save_known_branches(name, n_files, branches):
    os.makedirs(SLIM_DIR, exist_ok=True)
    with open(_branches_file(name, n_files), "w") as f:
        json.dump(sorted(branches), f, indent=2)


# ── Auto-detection ───────────────────────────────────────────────────────────

def detect_used_branches(expressions, tree_file):
    """Auto-detect tree branches referenced in C++ expression strings.

    Extracts all C-identifier tokens from the collected expressions,
    then intersects with actual tree branch names.
    """
    df_raw = ROOT.RDataFrame("Events", tree_file)
    tree_branches = set(str(c) for c in df_raw.GetColumnNames())

    all_exprs = " ".join(expressions)
    tokens = set(re.findall(r'[A-Za-z_][A-Za-z0-9_]*', all_exprs))

    used = tokens & tree_branches

    for e in ("genWeight", "run", "luminosityBlock"):
        if e in tree_branches:
            used.add(e)

    return sorted(used)


# ── Skim (with incremental friend-tree support) ─────────────────────────────

def run_skim(mc_dict, branches, n_files_map, mc_files_map=None, force=False):
    """Write slim ROOT files using detected branch list.

    Parameters
    ----------
    mc_dict : dict[str, RDataFrame]
        Sample name -> RDataFrame (with all Define/Filter chains).
    branches : list[str]
        All branch names needed by the analysis.
    n_files_map : dict[str, int]
        Sample name -> number of files used.
    mc_files_map : dict[str, list[str]] or None
        Sample name -> original EOS file paths (needed for incremental).
    force : bool
        If True, full rebuild (overwrite existing slim files).
    """
    os.makedirs(SLIM_DIR, exist_ok=True)
    print(f"\n[skim] Output directory: {SLIM_DIR}")
    print(f"[skim] Total needed: {len(branches)} branches")

    for name, df in mc_dict.items():
        n_files = n_files_map.get(name, 0)
        slim = _slim_path(name, n_files)
        friends = _friends_path(name, n_files)
        known = _load_known_branches(name, n_files)

        if force or not os.path.exists(slim):
            t0 = time.time()
            df.Snapshot("Events", slim, branches)
            sz = os.path.getsize(slim) / 1e6
            print(f"  {name:55s} saved ({sz:.1f} MB, {time.time() - t0:.1f}s) [{n_files} files]")
            _save_known_branches(name, n_files, branches)
            if os.path.exists(friends):
                os.remove(friends)
            continue

        # Incremental: find new branches
        new_branches = sorted(set(branches) - known)
        if not new_branches:
            sz = os.path.getsize(slim) / 1e6
            print(f"  {name:55s} up to date ({sz:.1f} MB) [{n_files} files]")
            continue

        print(f"  {name:55s} adding {len(new_branches)} new branches:")
        for b in new_branches:
            print(f"    + {b}")

        if mc_files_map and name in mc_files_map:
            eos_df = ROOT.RDataFrame("Events", mc_files_map[name])
        else:
            print(f"    ERROR: no EOS files for {name}, use force=True for full rebuild")
            continue

        t0 = time.time()
        eos_df.Snapshot("Events", friends, new_branches)
        sz = os.path.getsize(friends) / 1e6
        print(f"    friends saved ({sz:.1f} MB, {time.time() - t0:.1f}s)")
        _save_known_branches(name, n_files, known | set(new_branches))

    total = sum(
        os.path.getsize(os.path.join(SLIM_DIR, f)) / 1e6
        for f in os.listdir(SLIM_DIR)
        if f.endswith(".root")
    )
    print(f"\n[skim] Total disk: {total:.0f} MB in {SLIM_DIR}")


def run_skim_pipeline(mc, mc_files_map, exprs, args):
    """Detect branches from an explicit expression list and run the skim.

    Parameters
    ----------
    mc : dict[str, RDataFrame]
        Sample name -> weighted RDataFrame with all analysis columns defined.
    mc_files_map : dict[str, list[str]]
        Sample name -> original EOS file paths (for incremental skim).
    exprs : list[str]
        All C++ expression strings used by the analysis.
    args : argparse.Namespace
        Must have: force (bool).
    """
    import sys
    first_name = next(iter(mc_files_map))
    branches = detect_used_branches(exprs, mc_files_map[first_name][0])
    mc_branches = [b for b in branches if b not in ("run", "luminosityBlock")]
    print(f"[skim] Auto-detected {len(mc_branches)} branches from {len(exprs)} expressions")
    for b in mc_branches:
        print(f"  {b}")
    n_files_map = {name: len(files) for name, files in mc_files_map.items()}
    run_skim(mc, mc_branches, n_files_map, mc_files_map=mc_files_map,
             force=getattr(args, "force", False))
    print("Skim complete. Re-run without --skim to use slim files.")
    sys.exit(0)


# ── Slim file completeness check ──────────────────────────────────────────

def has_complete_slim(mc, mc_files_map):
    """Check if ALL samples have slim files with derived columns."""
    for name in mc:
        n_files = len(mc_files_map.get(name, []))
        slim = _slim_path(name, n_files)
        if not os.path.exists(slim) or not has_complete_slim_sample(slim):
            return False
    return True


def has_complete_slim_sample(slim_path):
    """Check if a single slim file contains derived columns (not just raw branches)."""
    if not os.path.exists(slim_path):
        return False
    try:
        df_check = ROOT.RDataFrame("Events", slim_path)
        cols = set(str(c) for c in df_check.GetColumnNames())
        return "ak4_pt0" in cols and "w" in cols
    except Exception:
        return False


# ── Auto-skim (called from load_and_run) ──────────────────────────────────

def collect_all_expressions():
    """Collect all C++ expressions by scanning analysis source + C++ files."""
    import importlib
    import inspect
    exprs = list(rdf_exprs)  # weight expressions from load_mc_samples

    # Scan Python source code for branch references
    for mod_path in ["analysis.definitions", "analysis.histograms",
                     "utils.triggers", "analysis.event_loop"]:
        try:
            mod = importlib.import_module(mod_path)
            src = inspect.getsource(mod)
            exprs.append(src)
        except Exception:
            pass

    # Scan C++ macro files for branch references
    ana_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for cpp_file in ["elements/GenMatching.C", "elements/RecoObjects.C",
                     "elements/common.h"]:
        cpp_path = os.path.join(ana_dir, cpp_file)
        if os.path.exists(cpp_path):
            with open(cpp_path) as f:
                exprs.append(f.read())

    # Add resolved branch prefixes (f-strings in definitions.py construct these at runtime)
    from analysis.constants import AK4, AK8, AK8_SGP, AK4_UPART
    for prefix in [AK4, AK8, AK8_SGP]:
        for suffix in ["pt", "eta", "phi", "mass", "msoftdrop"]:
            exprs.append(f"{prefix}_{suffix}")
        exprs.append(f"n{prefix}")
    # UParT branches
    for prob in ["probb", "probc", "probg", "probuds", "problepb", "probtaup", "probtaum"]:
        exprs.append(f"{AK4}_scoutUParT_{prob}")
    # AK8 ScoutGlobalParT branches
    from analysis.constants import AK8_SGP_PROB_NAMES
    for prob in AK8_SGP_PROB_NAMES:
        exprs.append(f"{AK8_SGP}_prob_{prob}")
    exprs.append(f"{AK8_SGP}_massCorrGeneric")
    exprs.append(f"{AK8_SGP}_massCorrResonance")

    return exprs


def ensure_slim(mc, mc_files_map, expressions=None):
    """Auto-create slim files with ALL columns (raw + derived + weights).

    Called after define_kinematics/define_gen_columns/define_gen_matched_ak4.
    Snapshots the fully-defined RDataFrame so subsequent runs skip Define chains.
    """
    if not mc_files_map:
        return

    n_files_map = {name: len(files) for name, files in mc_files_map.items()}

    samples_to_skim = {}
    for name in mc:
        n_files = n_files_map.get(name, 0)
        slim = _slim_path(name, n_files)
        if not os.path.exists(slim):
            samples_to_skim[name] = "missing"

    if not samples_to_skim:
        return

    print(f"\n[auto-skim] {len(samples_to_skim)} sample(s) need skimming:")
    for name, reason in samples_to_skim.items():
        print(f"  {name}: {reason}")

    os.makedirs(SLIM_DIR, exist_ok=True)

    for name in samples_to_skim:
        n_files = n_files_map[name]
        slim = _slim_path(name, n_files)
        df = mc[name]
        # Get all columns (raw + derived), excluding non-serializable struct types
        all_cols = []
        for c in sorted(str(c) for c in df.GetColumnNames()):
            try:
                ct = df.GetColumnType(c)
                # Skip struct/class types (GenHHResult, TLorentzVector, etc.)
                if any(x in ct for x in ["Result", "Ana::", "TLorentz", "vector<", "RVec"]):
                    continue
                all_cols.append(c)
            except Exception:
                continue
        t0 = time.time()
        df.Snapshot("Events", slim, all_cols)
        sz = os.path.getsize(slim) / 1e6
        print(f"  {name[:55]:55s} saved ({sz:.1f} MB, {time.time() - t0:.1f}s) [{n_files} files, {len(all_cols)} columns]")
        _save_known_branches(name, n_files, all_cols)

    total = sum(
        os.path.getsize(os.path.join(SLIM_DIR, f)) / 1e6
        for f in os.listdir(SLIM_DIR)
        if f.endswith(".root")
    )
    print(f"\n[auto-skim] Total disk: {total:.0f} MB in {SLIM_DIR}")

    # Reload from slim
    for name in samples_to_skim:
        n_files = n_files_map[name]
        mc[name] = _load_from_slim(name, n_files)


def _load_from_slim(name, n_files):
    """Load RDataFrame from slim file (+ friends if present)."""
    slim = _slim_path(name, n_files)
    friends = _friends_path(name, n_files)

    chain = ROOT.TChain("Events")
    chain.Add(slim)
    sz = os.path.getsize(slim) / 1e6

    if os.path.exists(friends):
        fsz = os.path.getsize(friends) / 1e6
        friend_chain = ROOT.TChain("Events")
        friend_chain.Add(friends)
        chain.AddFriend(friend_chain)
        print(f"  {name:55s} (slim {sz:.1f}MB + friends {fsz:.1f}MB) [{n_files} files]")
    else:
        print(f"  {name:55s} (slim, {sz:.1f} MB) [{n_files} files]")

    return ROOT.RDataFrame(chain)


# ── Loading (slim + friends, or fallback to EOS) ────────────────────────────

def load_slim_or_eos(name, eos_files, max_files=None):
    """Return RDataFrame from slim file if available, else from EOS.

    Checks for slim file with _{N}files suffix matching the actual file count.
    If a _friends.root file exists, it is attached via AddFriend.
    """
    files = eos_files[:max_files] if max_files else eos_files
    n_files = len(files)
    slim = _slim_path(name, n_files)

    if os.path.exists(slim):
        return _load_from_slim(name, n_files)

    print(f"  {name:55s} ({n_files} files, from EOS)")
    return ROOT.RDataFrame("Events", files)
