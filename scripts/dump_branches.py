#!/usr/bin/env python3
"""List branches from ROOT files specified in config/samples.yaml.

Saves one file per sample to output/branches/ with metadata header.
Filename: {group}_{sample_short}_{production_id}.txt

Usage:
    python scripts/list_branches.py                            # all samples
    python scripts/list_branches.py --group HHbbtt              # only signal
    python scripts/list_branches.py --filter "Scouting*"        # grep pattern
    python scripts/list_branches.py --config config/samples.yaml
    python scripts/list_branches.py --diff config/samples_copy.yaml
"""

import argparse, os, fnmatch, datetime, yaml, ROOT
ROOT.gROOT.SetBatch(True)

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "output", "branches")


def get_first_file(base_dir, sample_name, dataset_tag, production_id):
    """Find the first .root file for a sample."""
    sample_dir = os.path.join(base_dir, sample_name)
    if not os.path.isdir(sample_dir):
        return None, None, None

    if dataset_tag:
        tag_dirs = [os.path.join(sample_dir, dataset_tag)]
    else:
        tag_dirs = sorted([os.path.join(sample_dir, d)
                          for d in os.listdir(sample_dir)
                          if os.path.isdir(os.path.join(sample_dir, d))])

    for tag_dir in reversed(tag_dirs):
        if not os.path.isdir(tag_dir):
            continue
        found_tag = os.path.basename(tag_dir)
        if production_id:
            prod_dirs = [os.path.join(tag_dir, production_id)]
        else:
            prod_dirs = sorted([os.path.join(tag_dir, d)
                               for d in os.listdir(tag_dir)
                               if os.path.isdir(os.path.join(tag_dir, d))])

        for prod_dir in reversed(prod_dirs):
            if not os.path.isdir(prod_dir):
                continue
            found_prod = os.path.basename(prod_dir)
            for sub in sorted(os.listdir(prod_dir)):
                sub_path = os.path.join(prod_dir, sub)
                if os.path.isdir(sub_path):
                    roots = sorted([f for f in os.listdir(sub_path) if f.endswith(".root")])
                    if roots:
                        return os.path.join(sub_path, roots[0]), found_tag, found_prod
    return None, None, None


def list_branches(filepath, tree_name="Events"):
    """Get sorted list of branch names from a ROOT file."""
    f = ROOT.TFile.Open(filepath)
    if not f or f.IsZombie():
        return []
    t = f.Get(tree_name)
    if not t:
        return []
    branches = sorted([b.GetName() for b in t.GetListOfBranches()])
    f.Close()
    return branches


def nanoaod_version(base_dir):
    """Extract NanoAOD version from base_dir path."""
    for part in base_dir.split("/"):
        if "NanoAOD" in part:
            return part
    return "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/samples.yaml")
    parser.add_argument("--group", default=None,
                        help="Only process this sample group (e.g. HHbbtt, DY, TT, QCD)")
    parser.add_argument("--filter", default=None,
                        help="Only list branches matching this pattern (fnmatch)")
    parser.add_argument("--diff", default=None,
                        help="Compare branches against samples in another config")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    base_dir = cfg["settings"]["base_dir"]
    nano_ver = nanoaod_version(base_dir)
    samples = cfg.get("samples", {})

    os.makedirs(OUT_DIR, exist_ok=True)

    all_results = {}  # sample_name -> (branches, filepath, group, tag, prod_id)

    for group_name, group_samples in samples.items():
        if args.group and args.group != group_name:
            continue
        for sample_name, sample_cfg in group_samples.items():
            dataset_tag = sample_cfg.get("dataset_tag")
            production_id = sample_cfg.get("production_id")
            filepath, found_tag, found_prod = get_first_file(
                base_dir, sample_name, dataset_tag, production_id)

            if not filepath:
                print(f"  [{group_name}] {sample_name}: NO FILES FOUND")
                continue

            branches = list_branches(filepath)
            all_results[sample_name] = (branches, filepath, group_name, found_tag, found_prod)

            # Build output filename: {group}_{sample_short}_{production_id}.txt
            prod_str = found_prod or "unknown"
            out_name = f"{group_name}_{sample_name}_{prod_str}.txt"
            out_path = os.path.join(OUT_DIR, out_name)

            filtered = branches
            if args.filter:
                filtered = [b for b in branches if fnmatch.fnmatch(b, args.filter)]

            with open(out_path, "w") as fout:
                fout.write(f"# Group: {group_name}\n")
                fout.write(f"# Sample: {sample_name}\n")
                fout.write(f"# Dataset tag: {found_tag}\n")
                fout.write(f"# Production ID: {prod_str}\n")
                fout.write(f"# NanoAOD version: {nano_ver}\n")
                fout.write(f"# File: {filepath}\n")
                fout.write(f"# Total branches: {len(branches)}\n")
                if args.filter:
                    fout.write(f"# Filter: {args.filter} ({len(filtered)} matched)\n")
                fout.write(f"# Date: {datetime.date.today()}\n")
                fout.write("#\n")
                for b in filtered:
                    fout.write(b + "\n")

            print(f"  [{group_name}] {sample_name} -> {out_path} ({len(filtered)} branches)")

    # Diff mode
    if args.diff and all_results:
        with open(args.diff) as f:
            cfg2 = yaml.safe_load(f)
        base_dir2 = cfg2["settings"]["base_dir"]
        nano_ver2 = nanoaod_version(base_dir2)
        print(f"\n=== Diff: {nano_ver} vs {nano_ver2} ===")

        for gn, gs in cfg2.get("samples", {}).items():
            for sn, sc in gs.items():
                fp2, _, _ = get_first_file(base_dir2, sn, sc.get("dataset_tag"), sc.get("production_id"))
                if not fp2:
                    continue
                branches2 = set(list_branches(fp2))
                # Compare against first matching sample in current config
                for sn1, (br1, _, _, _, _) in all_results.items():
                    branches1 = set(br1)
                    added = sorted(branches1 - branches2)
                    removed = sorted(branches2 - branches1)
                    diff_path = os.path.join(OUT_DIR, f"diff_{nano_ver}_vs_{nano_ver2}.txt")
                    with open(diff_path, "w") as fd:
                        fd.write(f"# Diff: {nano_ver} vs {nano_ver2}\n")
                        fd.write(f"# Sample A: {sn1}\n")
                        fd.write(f"# Sample B: {sn}\n")
                        fd.write(f"# Date: {datetime.date.today()}\n")
                        fd.write(f"#\n")
                        fd.write(f"# + Added in {nano_ver} ({len(added)}):\n")
                        for b in added:
                            fd.write(f"+ {b}\n")
                        fd.write(f"#\n# - Removed from {nano_ver2} ({len(removed)}):\n")
                        for b in removed:
                            fd.write(f"- {b}\n")
                    print(f"  Diff saved: {diff_path}")
                    print(f"    + {len(added)} added, - {len(removed)} removed")
                    break
                break

    print(f"\nAll branch files saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
