// elements/RecoObjects.C
//
// Reco-level object selection for all processes in the HH->bb tautau
// scouting analysis.
//
// Functions provided (all in namespace Ana):
//
//   BvsAll()            — compute BvsAll b-tag discriminant per jet
//   findDijetInWindow() — find jet pair in mass window closest to target
//
// Load in cutflow_TrigEff.py after elements/GenMatching.C:
//   ROOT.gROOT.LoadMacro("elements/RecoObjects.C+")
//

#include "common.h"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Ana {


// ═════════════════════════════════════════════════════════════════════════════
//  BvsAll discriminant  (Signal + all processes)
// ═════════════════════════════════════════════════════════════════════════════
//
// Computes the BvsAll b-tag discriminant for each jet:
//   BvsAll = prob_b / (prob_b + prob_c + prob_cc + prob_g + prob_uds + prob_undef)
//
// Note: prob_bb is intentionally excluded — it is an AK8 fat-jet variable
// and is not meaningful for AK4 scouting jets.
//
// This mirrors the Define() expressions currently inlined in cutflow_TrigEff.py,
// moved here for reuse and clarity.
//
// Usage:
//   df.Define("pnet_BvsAll", "Ana::BvsAll("
//       "ScoutingPFJetRecluster_particleNet_prob_b, "
//       "ScoutingPFJetRecluster_particleNet_prob_c, "
//       "ScoutingPFJetRecluster_particleNet_prob_cc, "
//       "ScoutingPFJetRecluster_particleNet_prob_g, "
//       "ScoutingPFJetRecluster_particleNet_prob_uds, "
//       "ScoutingPFJetRecluster_particleNet_prob_undef)")
//
ROOT::RVec<float> BvsAll(
    const ROOT::RVec<float>& prob_b,
    const ROOT::RVec<float>& prob_c,
    const ROOT::RVec<float>& prob_cc,
    const ROOT::RVec<float>& prob_g,
    const ROOT::RVec<float>& prob_uds,
    const ROOT::RVec<float>& prob_undef)
{
    ROOT::RVec<float> result(prob_b.size(), 0.f);
    for (int i = 0; i < (int)prob_b.size(); ++i) {
        float denom = prob_b[i] + prob_c[i] + prob_cc[i]
                    + prob_g[i] + prob_uds[i] + prob_undef[i];
        result[i] = (denom > 0.f) ? prob_b[i] / denom : 0.f;
    }
    return result;
}


// ═════════════════════════════════════════════════════════════════════════════
//  Dijet mass-window pairing
// ═════════════════════════════════════════════════════════════════════════════
//
// Finds the jet pair with invariant mass inside [mass_lo, mass_hi] that is
// closest to `target`. Optionally excludes jet indices (e.g. jets already
// used by the H→bb candidate when searching for H→ττ).
//
// Used for H→bb (mass window [100,150], target 125) and H→ττ (mass window
// [40,150], target 80).
//
// Cut values come from config/regions.yaml via Python interpolation:
//   reg = CUTS_REG["hbb"]
//   df.Define("hbb_pair",
//       f"Ana::findDijetInWindow(ScoutingPFJetRecluster_pt, "
//       f"ScoutingPFJetRecluster_eta, ScoutingPFJetRecluster_phi, "
//       f"ScoutingPFJetRecluster_mass, "
//       f"{reg['mass_lo']}f, {reg['mass_hi']}f, {reg['target']}f)")
//
DijetPair findDijetInWindow(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi,
    const ROOT::RVec<float>& mass,
    float mass_lo, float mass_hi, float target,
    const std::vector<int>& exclude = {})
{
    DijetPair best;
    float bestDist = 1e9f;
    int n = (int)pt.size();

    for (int i = 0; i < n; ++i) {
        bool skip_i = false;
        for (int ex : exclude) { if (i == ex) { skip_i = true; break; } }
        if (skip_i) continue;

        for (int j = i + 1; j < n; ++j) {
            bool skip_j = false;
            for (int ex : exclude) { if (j == ex) { skip_j = true; break; } }
            if (skip_j) continue;

            TLorentzVector v1, v2;
            v1.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
            v2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
            float mjj = (float)(v1 + v2).M();

            if (mjj >= mass_lo && mjj <= mass_hi) {
                float dist = std::abs(mjj - target);
                if (dist < bestDist) {
                    bestDist = dist;
                    best.i1 = i;
                    best.i2 = j;
                    best.mass = mjj;
                }
            }
        }
    }
    return best;
}


// ═════════════════════════════════════════════════════════════════════════════
//  Muon selection  (ScoutingMuonVtx collection)
// ═════════════════════════════════════════════════════════════════════════════
//
// Selects the highest-pT ScoutingMuonVtx passing quality cuts.
// Returns collection index or -1 if none pass.
//
// Uses combined relative isolation: (trackIso + ecalIso + hcalIso) / pT
//
// Sources:
//   AN-2025/103 Table 6.7 (pT, eta, dxy, dz)
//   CMS SWGuideMuonIsolation TWiki (combined isolation)
//   Legacy code: Ana/elements/RecoObjects.C (SelectMuon)
//
// TODO: Add normchi2 < 10         [CMS Tight ID, EXO-19-018]
// TODO: Add nValidPixelHits > 0   [CMS Tight ID, EXO-19-018]
// TODO: Add nTrackerLayers > 5    [CMS Tight ID, EXO-19-018]
//
int SelectMuon(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& trackIso,
    const ROOT::RVec<float>& ecalIso,
    const ROOT::RVec<float>& hcalIso,
    const ROOT::RVec<float>& trk_dxy,
    const ROOT::RVec<float>& trk_dz,
    float pt_min,
    float eta_max,
    float iso_max,
    float dxy_max,
    float dz_max)
{
    int best = -1;
    float best_pt = -1.f;
    for (int i = 0; i < (int)pt.size(); ++i) {
        if (pt[i]                  < pt_min)  continue;
        if (std::abs(eta[i])       > eta_max) continue;
        if (pt[i] > 0 && (trackIso[i] + ecalIso[i] + hcalIso[i]) / pt[i] > iso_max) continue;
        if (std::abs(trk_dxy[i])   > dxy_max) continue;
        if (std::abs(trk_dz[i])    > dz_max)  continue;
        if (pt[i] > best_pt) { best_pt = pt[i]; best = i; }
    }
    return best;
}


// ═════════════════════════════════════════════════════════════════════════════
//  Electron selection  (ScoutingElectron collection)
// ═════════════════════════════════════════════════════════════════════════════
//
// Selects the highest-pT ScoutingElectron passing quality cuts.
// ECAL crack (crack_lo < |eta| < crack_hi) is excluded.
// Uses combined relative isolation: (trackIso + ecalIso + hcalIso) / pT
//
// Sources:
//   AN-2025/103 Table 6.8 (pT, eta, d0, dz)
//   CMS EgammaPublicData TWiki (crack veto)
//   Legacy code: Ana/elements/RecoObjects.C (SelectElectron)
//
// TODO: Add hOverE < 0.15                          [CMS cut-based ID]
// TODO: Add sigmaIetaIeta barrel < 0.011 / endcap < 0.031  [CMS Loose ID]
// TODO: Add |dEtaIn| barrel < 0.007 / endcap < 0.009       [CMS Loose ID]
// TODO: Add |dPhiIn| barrel < 0.15 / endcap < 0.10         [CMS Loose ID]
//
int SelectElectron(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& trackIso,
    const ROOT::RVec<float>& ecalIso,
    const ROOT::RVec<float>& hcalIso,
    const ROOT::RVec<float>& bestTrack_d0,
    const ROOT::RVec<float>& bestTrack_dz,
    float pt_min,
    float eta_max,
    float crack_lo,
    float crack_hi,
    float iso_max,
    float d0_max,
    float dz_max)
{
    int best = -1;
    float best_pt = -1.f;
    for (int i = 0; i < (int)pt.size(); ++i) {
        float aeta = std::abs(eta[i]);
        if (pt[i]  < pt_min)  continue;
        if (aeta   > eta_max) continue;
        if (aeta > crack_lo && aeta < crack_hi) continue;  // ECAL crack veto
        if (pt[i] > 0 && (trackIso[i] + ecalIso[i] + hcalIso[i]) / pt[i] > iso_max) continue;
        if (std::abs(bestTrack_d0[i]) > d0_max) continue;
        if (std::abs(bestTrack_dz[i]) > dz_max) continue;
        if (pt[i] > best_pt) { best_pt = pt[i]; best = i; }
    }
    return best;
}


} // namespace Ana
