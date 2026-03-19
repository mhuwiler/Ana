#!/usr/bin/env python3
"""
Quick XCache connectivity and performance test.

Usage:
    python test_xcache.py                    # default: 1 file, 1 thread
    python test_xcache.py --nthreads 32      # test with ImplicitMT
    python test_xcache.py --nfiles 10        # test with 10 files
    python test_xcache.py --nfiles 10 --nthreads 32  # full test
"""
import argparse
import time
import sys
import subprocess
import json

XCACHE_PREFIX = "root://xcache.cms.rcac.purdue.edu/"
DASGOCLIENT = "/cvmfs/cms.cern.ch/common/dasgoclient"
SCOUTING_DATA_JSON = "/home/das214/HHtobbtautau/Run3_nano_submission/datasets/Scouting_DATA.json"


def get_data_files(max_files=1):
    """Get scouting data file paths from DAS via XCache."""
    with open(SCOUTING_DATA_JSON) as f:
        datasets = json.load(f)

    # Use Run2024C (smallest era) for testing
    das_path = datasets.get("2024", {}).get("Run2024C")
    if not das_path:
        print("ERROR: Could not find Run2024C in Scouting_DATA.json")
        sys.exit(1)

    print(f"Querying DAS for Run2024C: {das_path} ...")
    t0 = time.time()
    result = subprocess.run(
        [DASGOCLIENT, "-query", f"file dataset={das_path}"],
        capture_output=True, text=True, check=True,
    )
    files = sorted([XCACHE_PREFIX + f.strip()
                    for f in result.stdout.strip().split("\n") if f.strip()])
    dt = time.time() - t0
    print(f"  DAS returned {len(files)} files in {dt:.1f}s")
    return files[:max_files]


def test_xcache(files, nthreads=1):
    """Run a series of tests on XCache files using ROOT RDataFrame."""
    import ROOT

    if nthreads > 1:
        ROOT.ROOT.EnableImplicitMT(nthreads)
        print(f"ImplicitMT enabled with {nthreads} threads")
    else:
        print("ImplicitMT disabled (single-threaded)")

    print(f"\nTesting {len(files)} file(s):")
    for f in files:
        print(f"  {f}")
    print()

    # ── Test 1: Create RDataFrame (no I/O yet) ──
    print("=" * 60)
    print("Test 1: Create RDataFrame")
    t0 = time.time()
    df = ROOT.RDataFrame("Events", files)
    dt = time.time() - t0
    print(f"  RDataFrame created in {dt:.1f}s")

    # ── Test 2: Simple Count (triggers full event loop) ──
    print("=" * 60)
    print("Test 2: Count() — full event loop, minimal column reads")
    t0 = time.time()
    n = df.Count().GetValue()
    dt = time.time() - t0
    print(f"  {n:,} events in {dt:.1f}s  "
          f"({n/dt:,.0f} events/s)" if dt > 0 else f"  {n:,} events")

    # ── Test 3: 1D histogram of a single column ──
    print("=" * 60)
    print("Test 3: Histo1D('ScoutingPFJetRecluster_pt') — reads one branch")
    t0 = time.time()
    h = df.Histo1D(("h_pt", ";pT;", 100, 0, 500), "ScoutingPFJetRecluster_pt")
    h.GetValue()  # trigger evaluation
    dt = time.time() - t0
    entries = h.GetValue().GetEntries()
    print(f"  {entries:,.0f} entries in {dt:.1f}s")

    # ── Test 4: Take two columns (what lumi extraction does) ──
    print("=" * 60)
    print("Test 4: Take('run') + Take('luminosityBlock') — lumi-style")
    t0 = time.time()
    run_take = df.Take["unsigned int"]("run")
    ls_take = df.Take["unsigned int"]("luminosityBlock")
    ROOT.RDF.RunGraphs([run_take, ls_take])
    dt_rg = time.time() - t0
    print(f"  RunGraphs done in {dt_rg:.1f}s")

    t0 = time.time()
    runs = run_take.GetValue()
    ls = ls_take.GetValue()
    dt_get = time.time() - t0
    print(f"  GetValue: {len(runs):,} events in {dt_get:.1f}s")

    import numpy as np
    t0 = time.time()
    r = np.asarray(runs, dtype=np.int64)
    l = np.asarray(ls, dtype=np.int64)
    keys = r * 100000 + l
    unique = np.unique(keys)
    u_runs = unique // 100000
    unique_runs, counts = np.unique(u_runs, return_counts=True)
    dt_np = time.time() - t0
    print(f"  numpy: {len(unique):,} unique (run,LS) pairs "
          f"across {len(unique_runs)} runs in {dt_np:.2f}s")

    # ── Test 5: Multi-column read (simulates Phase 1) ──
    print("=" * 60)
    print("Test 5: Multiple histograms + counts — simulates Phase 1")
    cols = ["ScoutingPFJetRecluster_pt", "ScoutingPFJetRecluster_eta",
            "ScoutingPFJetRecluster_mass", "run", "luminosityBlock"]
    ptrs = []
    ptrs.append(df.Count())
    ptrs.append(df.Histo1D(("h1", "", 50, 0, 500), "ScoutingPFJetRecluster_pt"))
    ptrs.append(df.Histo1D(("h2", "", 50, -5, 5), "ScoutingPFJetRecluster_eta"))
    ptrs.append(df.Histo1D(("h3", "", 50, 0, 100), "ScoutingPFJetRecluster_mass"))
    ptrs.append(df.Histo1D(("h4", "", 10000, 378000, 388000), "run"))

    t0 = time.time()
    ROOT.RDF.RunGraphs(ptrs)
    dt = time.time() - t0
    print(f"  RunGraphs({len(ptrs)} actions) done in {dt:.1f}s")

    # ── Summary ──
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Files:   {len(files)}")
    print(f"  Threads: {nthreads}")
    print(f"  Events:  {n:,}")
    print(f"  XCache seems to be working." if n > 0 else "  WARNING: 0 events read!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test XCache performance")
    parser.add_argument("--nfiles", type=int, default=1,
                        help="Number of data files to test (default: 1)")
    parser.add_argument("--nthreads", type=int, default=1,
                        help="Number of ImplicitMT threads (default: 1, i.e. off)")
    args = parser.parse_args()

    files = get_data_files(max_files=args.nfiles)
    test_xcache(files, nthreads=args.nthreads)
