// elements/common.h
//
// Shared types and inline utility functions for the elements/ C++ module.
// Included by GenMatching.C and RecoObjects.C.
//
// Design rules:
//   - ONLY inline functions and struct/class definitions here.
//   - NO non-inline function definitions (avoids linker conflicts with
//     HHbbtautauAnaElements.C when both .so files are loaded via ACLiC).
//   - All functions live in namespace Ana.
//
// PDG IDs used throughout:
//   e=11, mu=13, tau=15, b=5, t=6, W=24, Z=23, H=25
//
#pragma once

#include <cmath>
#include <vector>
#include <limits>
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"

namespace Ana {

// ─────────────────────────────────────────────────────────────────────────────
// PDG mass constants [GeV]
// ─────────────────────────────────────────────────────────────────────────────

constexpr float Z_MASS   = 91.1876f;
constexpr float W_MASS   = 80.377f;
constexpr float TOP_MASS = 172.69f;
constexpr float H_MASS   = 125.25f;
constexpr float TAU_MASS = 1.77686f;

// ─────────────────────────────────────────────────────────────────────────────
// Result structs
// ─────────────────────────────────────────────────────────────────────────────

// Shared: dijet pair from mass-window pairing
struct DijetPair {
    int   i1   = -1;   // index of jet 1
    int   i2   = -1;   // index of jet 2
    float mass = -1.f; // invariant mass [GeV]
    bool  found() const { return i1 >= 0 && i2 >= 0; }
};

// DY: gen-level Z boson decay classification
// decayMode:
//   0 = unknown / not found
//   1 = Z -> e+ e-
//   2 = Z -> mu+ mu-
//   3 = Z -> tau tau (both taus hadronic)
//   4 = Z -> tau tau (one tau -> mu)
//   5 = Z -> tau tau (one tau -> e)
struct GenZResult {
    int z_idx    = -1;
    int lep1     = -1;   // index of first Z daughter
    int lep2     = -1;   // index of second Z daughter
    int decayMode = 0;
    bool found() const { return z_idx >= 0; }
};

// TT: gen-level top quark decay classification
// decayMode:
//   0 = unknown / not found
//   1 = fully hadronic  (both W -> qq)
//   2 = semi-leptonic   (one W -> lv, one W -> qq)
//   3 = fully leptonic  (both W -> lv)
struct GenTopResult {
    int top_idx        = -1;
    int antitop_idx    = -1;
    int b_from_top     = -1;   // b quark index from t -> W b
    int b_from_antitop = -1;   // b quark index from tbar -> W bbar
    int decayMode      = 0;
    bool found() const { return top_idx >= 0; }
};

// Signal: gen-level HH -> bb tautau decay classification
// Mirrors GenMatchingResult from HHbbtautauAnaElements.C.
// decayType (matches GenDecay enum in that file):
//   0 = None / not found
//   1 = tau_h tau_h  (both taus hadronic)
//   2 = tau_mu tau_h (one tau -> mu)
//   3 = tau_e  tau_h (one tau -> e)
struct GenHHResult {
    int Htob      = -1;   // gen index of H -> bb Higgs
    int Htotau    = -1;   // gen index of H -> tautau Higgs
    int b1        = -1;   // first b quark
    int b2        = -1;   // second b quark
    int tau1      = -1;   // first tau
    int tau2      = -1;   // second tau
    int mu        = -1;   // muon from tau (if tau_mu tau_h)
    int e         = -1;   // electron from tau (if tau_e tau_h)
    int decayType = 0;
    bool found() const { return Htob >= 0 && Htotau >= 0; }
};

// ─────────────────────────────────────────────────────────────────────────────
// Angular utilities and gen status flag helpers
//
// NOTE: These are defined in HHbbtautauAnaElements.C as well (in namespace Ana).
// When both .so files are loaded, Cling sees duplicate definitions.
// We OMIT them here to avoid conflicts — the legacy definitions are used.
// The .C files that need them include local copies in their anonymous namespace.
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// Gen particle tree traversal
// ─────────────────────────────────────────────────────────────────────────────

// Returns true if `descendant` has `ancestor` anywhere in its mother chain.
inline bool isDescendantOf(int descendant, int ancestor,
                            const ROOT::RVec<int>& mothers) {
    int current = descendant;
    int safety  = 0;
    while (current >= 0 && current < (int)mothers.size() && safety < 200) {
        if (current == ancestor) return true;
        int next = mothers[current];
        if (next == current) break;   // self-loop guard
        current = next;
        ++safety;
    }
    return false;
}

// ─────────────────────────────────────────────────────────────────────────────
// 4-vector helpers
// ─────────────────────────────────────────────────────────────────────────────

// Build TLorentzVector from a collection at index `idx`.
// Returns zero-vector if idx is out of range.
inline TLorentzVector getP4(int idx,
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi,
    const ROOT::RVec<float>& mass)
{
    TLorentzVector v;
    if (idx >= 0 && idx < (int)pt.size()) {
        v.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], mass[idx]);
    }
    return v;
}

} // namespace Ana
