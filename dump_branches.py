#!/usr/bin/env python3
"""
Dump branch names from NanoAOD ROOT files for each MC sample group and data.

Outputs one text file per source into output/:
  - output/branches_data.txt
  - output/branches_DY.txt
  - output/branches_TT.txt
  - output/branches_HHbbtt.txt
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ROOT
from scouting_utils.data import (
    ls_nanoaod_files_groups,
    load_scouting_data,
    BASE, GROUPS,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUT_DIR, exist_ok=True)


def get_branches(root_file_path):
    """Return sorted list of branch names from the Events tree."""
    tf = ROOT.TFile.Open(root_file_path)
    tree = tf.Get("Events")
    branches = sorted(b.GetName() for b in tree.GetListOfBranches())
    tf.Close()
    return branches


def write_branches(branches, outpath):
    with open(outpath, "w") as f:
        for b in branches:
            f.write(b + "\n")
    print(f"  Wrote {len(branches)} branches -> {outpath}")


# --- MC samples (one file per group) ---
_, group_files_by_sample = ls_nanoaod_files_groups(
    base_dir=BASE, groups=GROUPS, verbose=0,
)

for group_name, samples_dict in group_files_by_sample.items():
    # Pick the first file from the first sample in this group
    first_sample = next(iter(samples_dict))
    first_file = samples_dict[first_sample][0]
    print(f"[{group_name}] Using: {first_file}")
    branches = get_branches(first_file)
    write_branches(branches, os.path.join(OUT_DIR, f"branches_{group_name}.txt"))

# --- Data ---
print("\n[Data]")
data_files, _ = load_scouting_data(years=["2024"])
print(f"  Using: {data_files[0]}")
branches = get_branches(data_files[0])
write_branches(branches, os.path.join(OUT_DIR, "branches_data.txt"))

print("\nDone!")
