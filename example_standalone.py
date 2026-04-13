#!/usr/bin/env python3
"""
Standalone example: Access MC samples and ROOT branches directly.

This script shows how to use the sample config (config/samples.yaml)
to load ROOT files into RDataFrames without the full analysis framework.
Useful for users who want to develop their own tools or just explore the data.

Usage:
    python example_standalone.py
"""

import ROOT
import numpy as np

# ═══════════════════════════════════════════════════════════════════════════════
#  1. Load sample configuration
# ═══════════════════════════════════════════════════════════════════════════════

# The sample config lives in config/samples.yaml and is parsed by utils/data.py.
# It provides: BASE (EOS path), GROUPS (sample groups), XSEC (cross sections).
from utils.data import ls_nanoaod_files_groups, BASE, GROUPS, XSEC

# List available sample groups
print("Available sample groups:")
for group_name, samples in GROUPS.items():
    print(f"  {group_name}: {list(samples.keys())}")
print()

# ═══════════════════════════════════════════════════════════════════════════════
#  2. Get file paths for specific samples
# ═══════════════════════════════════════════════════════════════════════════════

# Get files for selected groups
selected_groups = {k: v for k, v in GROUPS.items() if k in ["HHbbtt", "DY"]}
group_files, group_files_by_sample = ls_nanoaod_files_groups(
    base_dir=BASE, groups=selected_groups, verbose=1)

# Print file paths
for group, samples in group_files_by_sample.items():
    print(f"\n{group}:")
    for sample_name, files in samples.items():
        print(f"  {sample_name}: {len(files)} files")
        if files:
            print(f"    First file: {files[0]}")

# Create RDataFrames directly

# Pick one sample — limit to 1 file for speed
signal_files = group_files_by_sample["HHbbtt"]
sample_name = list(signal_files.keys())[0]
files = signal_files[sample_name][:1]  # just 1 file

print(f"\nLoading {sample_name} ({len(files)} file(s))...")
df = ROOT.RDataFrame("Events", files)

# ═══════════════════════════════════════════════════════════════════════════════
#  4. Explore available branches
# ═══════════════════════════════════════════════════════════════════════════════

all_cols = sorted(str(c) for c in df.GetColumnNames())
print(f"\nTotal branches: {len(all_cols)}")

# Filter by category
jet_branches = [c for c in all_cols if "ScoutingPFJetRecluster2" in c]
muon_branches = [c for c in all_cols if "ScoutingMuonVtx_" in c]
electron_branches = [c for c in all_cols if "ScoutingElectron_" in c]
gen_branches = [c for c in all_cols if "GenPart_" in c]

print(f"  AK4 jet branches: {len(jet_branches)}")
print(f"  Muon branches: {len(muon_branches)}")
print(f"  Electron branches: {len(electron_branches)}")
print(f"  GenPart branches: {len(gen_branches)}")

# ═══════════════════════════════════════════════════════════════════════════════
#  5. Define columns and make histograms (standard RDataFrame usage)
# ═══════════════════════════════════════════════════════════════════════════════

AK4 = "ScoutingPFJetRecluster2"

# Define some basic variables
df = (df
    .Define("lead_jet_pt", f"{AK4}_pt[0]")
    .Define("sublead_jet_pt", f"{AK4}_pt[1]")
    .Define("nJets", f"n{AK4}")
    .Define("HT", f"Sum({AK4}_pt)")
)

# Book histograms
h_pt = df.Histo1D(("h_pt", ";Leading jet p_{T} [GeV];Events", 50, 0, 250), "lead_jet_pt")
h_ht = df.Histo1D(("h_ht", ";H_{T} [GeV];Events", 100, 0, 2000), "HT")
h_njets = df.Histo1D(("h_njets", ";Number of jets;Events", 16, 0, 16), "nJets")

# Run event loop (lazy — executes only when you access results)
print(f"\nLeading jet pT: mean = {h_pt.GetMean():.1f} GeV, entries = {h_pt.GetEntries():.0f}")
print(f"HT: mean = {h_ht.GetMean():.1f} GeV")
print(f"nJets: mean = {h_njets.GetMean():.1f}")

# ═══════════════════════════════════════════════════════════════════════════════
#  6. Export to numpy (for ML, plotting, etc.)
# ═══════════════════════════════════════════════════════════════════════════════

# Export specific columns to numpy arrays
data = df.AsNumpy(["lead_jet_pt", "HT", "nJets"])
print(f"\nNumpy arrays:")
print(f"  lead_jet_pt: shape={data['lead_jet_pt'].shape}, "
      f"mean={data['lead_jet_pt'].mean():.1f}")
print(f"  HT: shape={data['HT'].shape}, mean={data['HT'].mean():.1f}")

# ═══════════════════════════════════════════════════════════════════════════════
#  7. Apply weights for MC normalisation
# ═══════════════════════════════════════════════════════════════════════════════

xsec_pb = XSEC[sample_name]
lumi_fb = 103.965  # fb^-1

# Compute sum of genWeights
sum_genw = df.Sum("genWeight").GetValue()
scale = xsec_pb * 1000.0 / sum_genw  # convert pb to fb

# Define per-event weight
df_weighted = df.Define("w", f"genWeight * {scale:.10e}")

# Weighted histogram
h_pt_w = df_weighted.Histo1D(
    ("h_pt_w", ";Leading jet p_{T} [GeV];Events", 50, 0, 250),
    "lead_jet_pt", "w")

# Scale by luminosity
h_pt_w.Scale(lumi_fb)
print(f"\nWeighted yield (×lumi): {h_pt_w.Integral():.2f} events")
print(f"  xsec = {xsec_pb} pb, sum_genw = {sum_genw:.4g}, lumi = {lumi_fb} fb^-1")

# ═══════════════════════════════════════════════════════════════════════════════
#  8. Access slim files (pre-processed with derived columns)
# ═══════════════════════════════════════════════════════════════════════════════

# If you've run `python skim.py`, slim files with pre-computed columns
# (ak4_pt0, HT, b_coi0_score, tau_coi0_score, etc.) are at:
SLIM_DIR = "/depot/cms/users/das214/tmp/slim"

import os
slim_files = [f for f in os.listdir(SLIM_DIR) if f.endswith(".root")]
if slim_files:
    print(f"\nSlim files available in {SLIM_DIR}:")
    for f in sorted(slim_files)[:5]:
        print(f"  {f}")
    print(f"  ... ({len(slim_files)} total)")

    # Load a slim file — all derived columns are pre-computed
    slim_path = os.path.join(SLIM_DIR, slim_files[0])
    df_slim = ROOT.RDataFrame("Events", slim_path)
    slim_cols = sorted(str(c) for c in df_slim.GetColumnNames())
    derived = [c for c in slim_cols if c.startswith(("ak4_", "b_coi", "tau_coi", "mbb", "HT"))]
    print(f"  Derived columns available: {derived[:10]}...")
else:
    print(f"\nNo slim files found. Run `python skim.py` first.")

print("\nDone!")
