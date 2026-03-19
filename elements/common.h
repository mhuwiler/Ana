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

// Shared: opposite-sign same-flavour dilepton pair
struct DileptonPair {
    int   i1         = -1;   // index of lepton 1 in collection
    int   i2         = -1;   // index of lepton 2 in collection
    float mass       = -1.f; // invariant mass [GeV]
    int   charge_sum = 99;   // i1.charge + i2.charge (0 for OS pair)
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
// Angular utilities
// ─────────────────────────────────────────────────────────────────────────────

inline double deltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    return dphi;
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);
    return std::sqrt(deta * deta + dphi * dphi);
}

// ─────────────────────────────────────────────────────────────────────────────
// Gen particle status flag helpers
// (same bit definitions as HHbbtautauAnaElements.C)
// ─────────────────────────────────────────────────────────────────────────────

inline bool isLastCopy(int flags)      { return (flags >> 13) & 1; }
inline bool isHardProcess(int flags)   { return (flags >>  7) & 1; }
inline bool fromHardProcess(int flags) { return (flags >>  8) & 1; }

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

// ─────────────────────────────────────────────────────────────────────────────
// Safe array access
// ─────────────────────────────────────────────────────────────────────────────

template<typename T>
inline T safeAt(const ROOT::RVec<T>& v, int idx, T fallback = T{}) {
    if (idx >= 0 && idx < (int)v.size()) return v[idx];
    return fallback;
}

} // namespace Ana
