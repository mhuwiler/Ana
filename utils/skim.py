"""
Snapshot/skim utilities for I/O optimization.

Writes slim ROOT files with only the branches used by the analysis.
Supports incremental branch addition via TTree::AddFriend — new branches
are stored in a separate _friends.root file without rebuilding the main slim.

Usage:
    # First skim (full):
    python cutflow_TrigEff.py --skim --max-mc-files 2 --no-data

    # After adding new Define() that uses a new branch:
    python cutflow_TrigEff.py --reslim --max-mc-files 2 --no-data
    # → only reads the NEW branch from EOS, adds to friends file

    # Full rebuild (merge friends into main):
    python cutflow_TrigEff.py --reslim --force --max-mc-files 2 --no-data
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


# ── Branch tracking (JSON sidecar files) ─────────────────────────────────────

def _branches_file(sample_name):
    return os.path.join(SLIM_DIR, f"{sample_name}.branches.json")


def _load_known_branches(sample_name):
    path = _branches_file(sample_name)
    if os.path.exists(path):
        with open(path) as f:
            return set(json.load(f))
    return set()


def _save_known_branches(sample_name, branches):
    with open(_branches_file(sample_name), "w") as f:
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

def run_skim(mc_dict, branches, mc_files_map=None, force=False):
    """Write slim ROOT files using detected branch list.

    On first run, writes all branches to slim/<sample>.root.
    On subsequent runs (--reslim), only writes NEW branches to
    slim/<sample>_friends.root via TTree::AddFriend pattern.
    With force=True, does a full rebuild merging everything.

    Parameters
    ----------
    mc_dict : dict[str, RDataFrame]
        Sample name -> RDataFrame (with all Define/Filter chains).
    branches : list[str]
        All branch names needed by the analysis.
    mc_files_map : dict[str, list[str]] or None
        Sample name -> original EOS file paths (needed for incremental).
    force : bool
        If True, full rebuild (overwrite existing slim files).
    """
    os.makedirs(SLIM_DIR, exist_ok=True)
    print(f"\n[skim] Output directory: {SLIM_DIR}")
    print(f"[skim] Total needed: {len(branches)} branches")

    for name, df in mc_dict.items():
        slim_path = os.path.join(SLIM_DIR, f"{name}.root")
        friends_path = os.path.join(SLIM_DIR, f"{name}_friends.root")
        known = _load_known_branches(name)

        if force or not os.path.exists(slim_path):
            # Full skim
            t0 = time.time()
            df.Snapshot("Events", slim_path, branches)
            sz = os.path.getsize(slim_path) / 1e6
            print(f"  {name:55s} saved ({sz:.1f} MB, {time.time() - t0:.1f}s)")
            _save_known_branches(name, branches)
            if os.path.exists(friends_path):
                os.remove(friends_path)
            continue

        # Incremental: find new branches
        new_branches = sorted(set(branches) - known)
        if not new_branches:
            sz = os.path.getsize(slim_path) / 1e6
            print(f"  {name:55s} up to date ({sz:.1f} MB)")
            continue

        print(f"  {name:55s} adding {len(new_branches)} new branches:")
        for b in new_branches:
            print(f"    + {b}")

        if mc_files_map and name in mc_files_map:
            eos_df = ROOT.RDataFrame("Events", mc_files_map[name])
        else:
            print(f"    ERROR: no EOS files for {name}, use --force for full rebuild")
            continue

        t0 = time.time()
        eos_df.Snapshot("Events", friends_path, new_branches)
        sz = os.path.getsize(friends_path) / 1e6
        print(f"    friends saved ({sz:.1f} MB, {time.time() - t0:.1f}s)")
        _save_known_branches(name, known | set(new_branches))

    total = sum(
        os.path.getsize(os.path.join(SLIM_DIR, f)) / 1e6
        for f in os.listdir(SLIM_DIR)
        if f.endswith(".root")
    )
    print(f"\n[skim] Total disk: {total:.0f} MB in {SLIM_DIR}")


# ── Loading (slim + friends, or fallback to EOS) ────────────────────────────

def load_slim_or_eos(name, eos_files, max_files=None):
    """Return RDataFrame from slim file if available, else from EOS.

    If a _friends.root file exists, it is attached via AddFriend
    so all branches (original + incremental) are available.
    """
    slim = os.path.join(SLIM_DIR, f"{name}.root")
    friends = os.path.join(SLIM_DIR, f"{name}_friends.root")

    if os.path.exists(slim):
        sz = os.path.getsize(slim) / 1e6
        chain = ROOT.TChain("Events")
        chain.Add(slim)
        if os.path.exists(friends):
            fsz = os.path.getsize(friends) / 1e6
            friend_chain = ROOT.TChain("Events")
            friend_chain.Add(friends)
            chain.AddFriend(friend_chain)
            print(f"  {name:55s} (slim {sz:.1f}MB + friends {fsz:.1f}MB)")
        else:
            print(f"  {name:55s} (slim, {sz:.1f} MB)")
        return ROOT.RDataFrame(chain)

    files = eos_files[:max_files] if max_files else eos_files
    print(f"  {name:55s} ({len(files)} files)")
    return ROOT.RDataFrame("Events", files)
