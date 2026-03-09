#!/usr/bin/env python3
"""
Generator-level analysis: match gen b-quarks to reco scouting jets
and evaluate ParticleNet b-tagging performance.

Only runs on signal MC samples (HH -> bb tau tau).

Usage:
    python genAna.py                        # default: light theme, MAX_EVENTS limit
    python genAna.py --theme dark           # dark theme
    python genAna.py --max-files 2          # quick test
    python genAna.py --all-events           # process all events
    python genAna.py --overwrite            # regenerate existing plots
"""

import os
import sys
import argparse
import resource

import numpy as np

parser = argparse.ArgumentParser(description="Gen-level b-jet matching analysis")
parser.add_argument("--overwrite", action="store_true",
                    help="Overwrite existing plots")
parser.add_argument("--theme", choices=["light", "dark"], default="light",
                    help="Plot colour theme (default: light)")
parser.add_argument("--max-files", type=int, default=0,
                    help="Max files per sample (0 = use MAX_EVENTS limit)")
parser.add_argument("--all-events", action="store_true",
                    help="Process all events (ignore MAX_EVENTS limit)")
parser.add_argument("--dr-threshold", type=float, default=0.4,
                    help="DeltaR threshold for gen-reco matching (default: 0.4)")
ARGS = parser.parse_args()

# Increase stack size
STACK_SIZE = 64 * 1024 * 1024
resource.setrlimit(resource.RLIMIT_STACK, (STACK_SIZE, resource.RLIM_INFINITY))

import ROOT
import anaConfig

from scouting_utils.data import (
    ls_nanoaod_files_groups,
    BASE, GROUPS, XSEC, MAX_EVENTS,
)
from scouting_utils.plotting import setup_style

import matplotlib.pyplot as plt
import mplhep as hep

# ──────────────────────────────────────────────────────────────────────────────
# ROOT setup
# ──────────────────────────────────────────────────────────────────────────────

ANA_DIR = os.path.expanduser(
    "/home/das214/HHtobbtautau/Scouting/CMSSW_15_0_15/src/Ana"
)
os.chdir(ANA_DIR)
sys.path.insert(0, ANA_DIR)

ROOT.gInterpreter.AddIncludePath(ANA_DIR)
ROOT.gROOT.LoadMacro("Particle.h+")
ROOT.gROOT.LoadMacro("HHbbtautauAnaElements.C+")

ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.ROOT.EnableImplicitMT(32)

setup_style(dark=(ARGS.theme == "dark"))

PLOT_DIR = os.path.join(ANA_DIR, "plots", ARGS.theme, "genAna")
os.makedirs(PLOT_DIR, exist_ok=True)


# ──────────────────────────────────────────────────────────────────────────────
# Load signal MC only
# ──────────────────────────────────────────────────────────────────────────────

# Only need the signal group
sig_groups = {"HHbbtt": GROUPS["HHbbtt"]}
_, group_files_by_sample = ls_nanoaod_files_groups(
    base_dir=BASE, groups=sig_groups, verbose=0,
)


def _limit_files(files, sample_name):
    good_files = [f for f in files if os.path.isfile(f)]
    if not good_files:
        raise FileNotFoundError(f"No valid files for {sample_name}")
    if ARGS.max_files > 0:
        good_files = good_files[:ARGS.max_files]
    elif not ARGS.all_events and MAX_EVENTS > 0:
        n_files = max(1, MAX_EVENTS // 10_000)
        good_files = good_files[:n_files]
    return good_files


# Build RDataFrames for signal samples
sig_dfs = {}
for samples_dict in group_files_by_sample.values():
    for sample_name, files in samples_dict.items():
        good_files = _limit_files(files, sample_name)
        df = ROOT.RDataFrame("Events", good_files)
        sig_dfs[sample_name] = df
        print(f"  {sample_name[:60]:60s}  ({len(good_files)} files)")


# ──────────────────────────────────────────────────────────────────────────────
# Gen-level matching + reco jet definitions
# ──────────────────────────────────────────────────────────────────────────────

gen_decay_expr = (
    "Ana::DecayGenMatching(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)"
)

DR_THRESH = ARGS.dr_threshold

# C++ helper: match a gen particle to the closest reco scouting jet within dR threshold.
# Returns the reco jet index, or -1 if no match found.
closest_match_cpp = f"""
int genMatchToRecoJet(int genIdx,
                      const ROOT::VecOps::RVec<float>& genPt,
                      const ROOT::VecOps::RVec<float>& genEta,
                      const ROOT::VecOps::RVec<float>& genPhi,
                      const ROOT::VecOps::RVec<float>& genMass,
                      const ROOT::VecOps::RVec<float>& recoPt,
                      const ROOT::VecOps::RVec<float>& recoEta,
                      const ROOT::VecOps::RVec<float>& recoPhi,
                      const ROOT::VecOps::RVec<float>& recoMass) {{
    if (genIdx < 0 || genIdx >= (int)genPt.size()) return -1;
    TLorentzVector gen;
    gen.SetPtEtaPhiM(genPt[genIdx], genEta[genIdx], genPhi[genIdx], genMass[genIdx]);
    double bestDR = {DR_THRESH};
    int best = -1;
    for (int i = 0; i < (int)recoPt.size(); i++) {{
        TLorentzVector reco;
        reco.SetPtEtaPhiM(recoPt[i], recoEta[i], recoPhi[i], recoMass[i]);
        double dR = gen.DeltaR(reco);
        if (dR < bestDR) {{
            bestDR = dR;
            best = i;
        }}
    }}
    return best;
}}
"""
ROOT.gInterpreter.Declare(closest_match_cpp)

# ParticleNet BvsAll discriminator
pnet_bvsall_expr = (
    "(ScoutingPFJetRecluster_particleNet_prob_b + ScoutingPFJetRecluster_particleNet_prob_bb)"
    " / (ScoutingPFJetRecluster_particleNet_prob_b + ScoutingPFJetRecluster_particleNet_prob_bb"
    " + ScoutingPFJetRecluster_particleNet_prob_c + ScoutingPFJetRecluster_particleNet_prob_cc"
    " + ScoutingPFJetRecluster_particleNet_prob_g + ScoutingPFJetRecluster_particleNet_prob_uds"
    " + ScoutingPFJetRecluster_particleNet_prob_undef)"
)

for name, df in sig_dfs.items():
    df = (df
        .Define("GenDecay", gen_decay_expr)
        .Define("pnet_BvsAll", pnet_bvsall_expr)
        # Match gen b-quarks (b1, b2) to reco scouting jets
        .Define("reco_b1_idx",
                "genMatchToRecoJet(GenDecay.b1, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass,"
                " ScoutingPFJetRecluster_pt, ScoutingPFJetRecluster_eta,"
                " ScoutingPFJetRecluster_phi, ScoutingPFJetRecluster_mass)")
        .Define("reco_b2_idx",
                "genMatchToRecoJet(GenDecay.b2, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass,"
                " ScoutingPFJetRecluster_pt, ScoutingPFJetRecluster_eta,"
                " ScoutingPFJetRecluster_phi, ScoutingPFJetRecluster_mass)")
        # b-scores for matched jets (-1 if no match)
        .Define("matched_b1_score", "reco_b1_idx >= 0 ? pnet_BvsAll[reco_b1_idx] : -1.0f")
        .Define("matched_b2_score", "reco_b2_idx >= 0 ? pnet_BvsAll[reco_b2_idx] : -1.0f")
        # pT of matched jets
        .Define("matched_b1_pt", "reco_b1_idx >= 0 ? ScoutingPFJetRecluster_pt[reco_b1_idx] : -1.0f")
        .Define("matched_b2_pt", "reco_b2_idx >= 0 ? ScoutingPFJetRecluster_pt[reco_b2_idx] : -1.0f")
        # dR between gen and matched reco
        .Define("matched_b1_dR",
                "reco_b1_idx >= 0 ? ROOT::VecOps::DeltaR("
                "GenPart_eta[GenDecay.b1], ScoutingPFJetRecluster_eta[reco_b1_idx],"
                "GenPart_phi[GenDecay.b1], ScoutingPFJetRecluster_phi[reco_b1_idx]) : -1.0f")
        .Define("matched_b2_dR",
                "reco_b2_idx >= 0 ? ROOT::VecOps::DeltaR("
                "GenPart_eta[GenDecay.b2], ScoutingPFJetRecluster_eta[reco_b2_idx],"
                "GenPart_phi[GenDecay.b2], ScoutingPFJetRecluster_phi[reco_b2_idx]) : -1.0f")
        # Gen b-quark pT
        .Define("gen_b1_pt", "GenDecay.b1 >= 0 ? GenPart_pt[GenDecay.b1] : -1.0f")
        .Define("gen_b2_pt", "GenDecay.b2 >= 0 ? GenPart_pt[GenDecay.b2] : -1.0f")
        # Decay type
        .Define("decayType", "(int)GenDecay.decayType")
        # Number of reco jets
        .Define("nJets", "(int)nScoutingPFJetRecluster")
        # All jet b-scores (for comparison)
        .Define("bsort_idx", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(pnet_BvsAll))")
        .Define("best_bscore", "pnet_BvsAll[bsort_idx[0]]")
        .Define("second_bscore", "bsort_idx.size() > 1 ? pnet_BvsAll[bsort_idx[1]] : -1.0f")
    )
    sig_dfs[name] = df


# ──────────────────────────────────────────────────────────────────────────────
# Book histograms and run
# ──────────────────────────────────────────────────────────────────────────────

print("\nBooking histograms...")

hists = {}
for name, df in sig_dfs.items():
    # Require at least 4 reco jets and valid gen decay
    df_cut = df.Filter("nJets >= 4 && GenDecay.b1 >= 0 && GenDecay.b2 >= 0")

    # Matched b-jet b-scores (only events with successful match)
    df_b1_matched = df_cut.Filter("reco_b1_idx >= 0")
    df_b2_matched = df_cut.Filter("reco_b2_idx >= 0")
    df_both_matched = df_cut.Filter("reco_b1_idx >= 0 && reco_b2_idx >= 0")

    h = {}
    # b-score distributions
    h["matched_b1_score"] = df_b1_matched.Histo1D(
        ("h_matched_b1_score", "Gen-matched b1 jet;ParticleNet BvsAll;Events", 50, 0, 1),
        "matched_b1_score")
    h["matched_b2_score"] = df_b2_matched.Histo1D(
        ("h_matched_b2_score", "Gen-matched b2 jet;ParticleNet BvsAll;Events", 50, 0, 1),
        "matched_b2_score")

    # Best reco b-score (regardless of gen matching)
    h["best_bscore"] = df_cut.Histo1D(
        ("h_best_bscore", "Best reco b-score;ParticleNet BvsAll;Events", 50, 0, 1),
        "best_bscore")
    h["second_bscore"] = df_cut.Histo1D(
        ("h_second_bscore", "2nd best reco b-score;ParticleNet BvsAll;Events", 50, 0, 1),
        "second_bscore")

    # All jets b-score (flattened)
    h["all_bscore"] = df_cut.Histo1D(
        ("h_all_bscore", "All jets b-score;ParticleNet BvsAll;Events", 50, 0, 1),
        "pnet_BvsAll")

    # dR between gen and matched reco
    h["matched_b1_dR"] = df_b1_matched.Histo1D(
        ("h_matched_b1_dR", "dR(gen b1, reco jet);#DeltaR;Events", 50, 0, 0.5),
        "matched_b1_dR")
    h["matched_b2_dR"] = df_b2_matched.Histo1D(
        ("h_matched_b2_dR", "dR(gen b2, reco jet);#DeltaR;Events", 50, 0, 0.5),
        "matched_b2_dR")

    # Gen b-quark pT
    h["gen_b1_pt"] = df_cut.Histo1D(
        ("h_gen_b1_pt", "Gen b1 pT;p_{T} [GeV];Events", 50, 0, 500),
        "gen_b1_pt")
    h["gen_b2_pt"] = df_cut.Histo1D(
        ("h_gen_b2_pt", "Gen b2 pT;p_{T} [GeV];Events", 50, 0, 500),
        "gen_b2_pt")

    # Matched reco jet pT
    h["matched_b1_pt"] = df_b1_matched.Histo1D(
        ("h_matched_b1_pt", "Matched reco b1 jet pT;p_{T} [GeV];Events", 50, 0, 500),
        "matched_b1_pt")
    h["matched_b2_pt"] = df_b2_matched.Histo1D(
        ("h_matched_b2_pt", "Matched reco b2 jet pT;p_{T} [GeV];Events", 50, 0, 500),
        "matched_b2_pt")

    # Matching efficiency: count events
    h["n_total"] = df_cut.Count()
    h["n_b1_matched"] = df_b1_matched.Count()
    h["n_b2_matched"] = df_b2_matched.Count()
    h["n_both_matched"] = df_both_matched.Count()

    hists[name] = h

# Run all histograms in one go
all_ptrs = []
for h in hists.values():
    for v in h.values():
        all_ptrs.append(v)

print(f"Running {len(all_ptrs)} histogram actions in parallel...")
ROOT.RDF.RunGraphs(all_ptrs)
print("Done.")


# ──────────────────────────────────────────────────────────────────────────────
# Print matching statistics
# ──────────────────────────────────────────────────────────────────────────────

for name, h in hists.items():
    n_total = h["n_total"].GetValue()
    n_b1 = h["n_b1_matched"].GetValue()
    n_b2 = h["n_b2_matched"].GetValue()
    n_both = h["n_both_matched"].GetValue()

    print(f"\n{'='*60}")
    print(f"Sample: {name}")
    print(f"{'='*60}")
    print(f"Events with >=4 jets and valid gen b1,b2: {n_total}")
    print(f"  b1 matched to reco jet (dR < {DR_THRESH}): {n_b1} ({100*n_b1/max(1,n_total):.1f}%)")
    print(f"  b2 matched to reco jet (dR < {DR_THRESH}): {n_b2} ({100*n_b2/max(1,n_total):.1f}%)")
    print(f"  Both matched:                              {n_both} ({100*n_both/max(1,n_total):.1f}%)")


# ──────────────────────────────────────────────────────────────────────────────
# Helper to extract numpy arrays from TH1
# ──────────────────────────────────────────────────────────────────────────────

def th1_to_np(h):
    nbins = h.GetNbinsX()
    edges = np.array([h.GetBinLowEdge(i+1) for i in range(nbins)] + [h.GetBinLowEdge(nbins+1)])
    vals = np.array([h.GetBinContent(i+1) for i in range(nbins)])
    errs = np.array([h.GetBinError(i+1) for i in range(nbins)])
    return edges, vals, errs


# ──────────────────────────────────────────────────────────────────────────────
# Plot
# ──────────────────────────────────────────────────────────────────────────────

print("\nMaking plots...")

for name, h in hists.items():
    short_name = name.split("_")[0]  # e.g. "GluGluHHto2B2Tau"

    # ── 1. b-score comparison: gen-matched b-jets vs all jets ──
    outpath = os.path.join(PLOT_DIR, "bscore_genmatched_vs_all.png")
    if os.path.exists(outpath) and not ARGS.overwrite:
        print(f"Skipping (exists): {outpath}")
    else:
        fig, ax = plt.subplots(figsize=(8, 6))

        edges, vals_all, _ = th1_to_np(h["all_bscore"].GetPtr())
        edges, vals_b1, _ = th1_to_np(h["matched_b1_score"].GetPtr())
        edges, vals_b2, _ = th1_to_np(h["matched_b2_score"].GetPtr())

        # Normalize to unity
        def norm(v):
            s = v.sum()
            return v / s if s > 0 else v

        ax.step(edges, np.r_[norm(vals_all), norm(vals_all)[-1]], where="post",
                linewidth=2, color="gray", label="All jets")
        ax.step(edges, np.r_[norm(vals_b1), norm(vals_b1)[-1]], where="post",
                linewidth=2, color="tab:red", label="Gen-matched b1")
        ax.step(edges, np.r_[norm(vals_b2), norm(vals_b2)[-1]], where="post",
                linewidth=2, color="tab:blue", label="Gen-matched b2")

        ax.set_xlabel("ParticleNet BvsAll score")
        ax.set_ylabel("Fraction of jets (normalised)")
        ax.set_title(f"Gen-matched b-jets vs all jets (dR < {DR_THRESH})")
        ax.legend()
        ax.grid(True, alpha=0.3)
        hep.cms.label("Simulation Preliminary", ax=ax, data=False)
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
        plt.close(fig)

    # ── 2. Best reco b-score vs gen-matched b-score ──
    outpath = os.path.join(PLOT_DIR, "bscore_best_reco_vs_genmatched.png")
    if os.path.exists(outpath) and not ARGS.overwrite:
        print(f"Skipping (exists): {outpath}")
    else:
        fig, ax = plt.subplots(figsize=(8, 6))

        edges, vals_best, _ = th1_to_np(h["best_bscore"].GetPtr())
        edges, vals_2nd, _ = th1_to_np(h["second_bscore"].GetPtr())
        edges, vals_b1, _ = th1_to_np(h["matched_b1_score"].GetPtr())
        edges, vals_b2, _ = th1_to_np(h["matched_b2_score"].GetPtr())

        ax.step(edges, np.r_[norm(vals_best), norm(vals_best)[-1]], where="post",
                linewidth=2, color="tab:orange", label="Best reco b-score")
        ax.step(edges, np.r_[norm(vals_2nd), norm(vals_2nd)[-1]], where="post",
                linewidth=2, color="tab:green", label="2nd best reco b-score")
        ax.step(edges, np.r_[norm(vals_b1), norm(vals_b1)[-1]], where="post",
                linewidth=2, color="tab:red", linestyle="--", label="Gen-matched b1")
        ax.step(edges, np.r_[norm(vals_b2), norm(vals_b2)[-1]], where="post",
                linewidth=2, color="tab:blue", linestyle="--", label="Gen-matched b2")

        ax.set_xlabel("ParticleNet BvsAll score")
        ax.set_ylabel("Fraction of events (normalised)")
        ax.set_title("Reco b-tag ranking vs gen-matched b-jets")
        ax.legend()
        ax.grid(True, alpha=0.3)
        hep.cms.label("Simulation Preliminary", ax=ax, data=False)
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
        plt.close(fig)

    # ── 3. dR distribution for gen-reco matching ──
    outpath = os.path.join(PLOT_DIR, "dR_gen_reco_match.png")
    if os.path.exists(outpath) and not ARGS.overwrite:
        print(f"Skipping (exists): {outpath}")
    else:
        fig, ax = plt.subplots(figsize=(8, 6))

        edges, vals_dr1, _ = th1_to_np(h["matched_b1_dR"].GetPtr())
        edges, vals_dr2, _ = th1_to_np(h["matched_b2_dR"].GetPtr())

        ax.step(edges, np.r_[vals_dr1, vals_dr1[-1]], where="post",
                linewidth=2, color="tab:red", label="b1 match dR")
        ax.step(edges, np.r_[vals_dr2, vals_dr2[-1]], where="post",
                linewidth=2, color="tab:blue", label="b2 match dR")

        ax.set_xlabel(r"$\Delta R$(gen b, reco jet)")
        ax.set_ylabel("Events")
        ax.set_title(f"Gen-reco matching quality (threshold dR < {DR_THRESH})")
        ax.legend()
        ax.grid(True, alpha=0.3)
        hep.cms.label("Simulation Preliminary", ax=ax, data=False)
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
        plt.close(fig)

    # ── 4. Gen b-quark pT vs matched reco jet pT ──
    outpath = os.path.join(PLOT_DIR, "pt_gen_vs_reco_matched.png")
    if os.path.exists(outpath) and not ARGS.overwrite:
        print(f"Skipping (exists): {outpath}")
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        for ax, label, gen_key, reco_key, color in [
            (ax1, "b1", "gen_b1_pt", "matched_b1_pt", "tab:red"),
            (ax2, "b2", "gen_b2_pt", "matched_b2_pt", "tab:blue"),
        ]:
            edges_g, vals_g, _ = th1_to_np(h[gen_key].GetPtr())
            edges_r, vals_r, _ = th1_to_np(h[reco_key].GetPtr())

            ax.step(edges_g, np.r_[vals_g, vals_g[-1]], where="post",
                    linewidth=2, color=color, label=f"Gen {label} pT")
            ax.step(edges_r, np.r_[vals_r, vals_r[-1]], where="post",
                    linewidth=2, color=color, linestyle="--", label=f"Matched reco {label} pT")

            ax.set_xlabel(r"$p_T$ [GeV]")
            ax.set_ylabel("Events")
            ax.set_title(f"{label} quark")
            ax.legend()
            ax.grid(True, alpha=0.3)

        fig.suptitle("Gen b-quark pT vs matched reco jet pT")
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
        plt.close(fig)

    # ── 5. b-score vs gen b pT (profile-like: binned average) ──
    # Use the matched b1 for this study
    outpath = os.path.join(PLOT_DIR, "bscore_vs_genpt.png")
    if os.path.exists(outpath) and not ARGS.overwrite:
        print(f"Skipping (exists): {outpath}")
    else:
        # Get b1 score and gen pT as numpy arrays via TProfile
        # Re-book a TProfile for b-score vs gen pT
        df_b1m = sig_dfs[name].Filter(
            "nJets >= 4 && GenDecay.b1 >= 0 && GenDecay.b2 >= 0 && reco_b1_idx >= 0"
        )
        prof_ptr = df_b1m.Profile1D(
            ("prof_bscore_vs_pt", "Mean b-score vs gen b pT;Gen b pT [GeV];Mean BvsAll",
             20, 0, 300),
            "gen_b1_pt", "matched_b1_score"
        )
        prof = prof_ptr.GetValue()  # trigger computation, get TProfile

        # Extract profile values
        nbins = prof.GetNbinsX()
        centers = np.array([prof.GetBinCenter(i+1) for i in range(nbins)])
        means = np.array([prof.GetBinContent(i+1) for i in range(nbins)])
        errs = np.array([prof.GetBinError(i+1) for i in range(nbins)])
        # Mask empty bins
        mask = np.array([prof.GetBinEntries(i+1) > 0 for i in range(nbins)])

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.errorbar(centers[mask], means[mask], yerr=errs[mask],
                    fmt="o-", color="tab:red", capsize=3, label="Gen-matched b1")
        ax.set_xlabel(r"Gen b-quark $p_T$ [GeV]")
        ax.set_ylabel("Mean ParticleNet BvsAll score")
        ax.set_title("b-tag score vs gen b-quark pT")
        ax.set_xlim(0, 300)
        ax.set_ylim(0, 1)
        ax.legend()
        ax.grid(True, alpha=0.3)
        hep.cms.label("Simulation Preliminary", ax=ax, data=False)
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
        plt.close(fig)

print("\nDone! Plots saved to:", PLOT_DIR)
