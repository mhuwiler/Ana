"""
Trigger bit definitions for scouting and parking analyses.

Contains:
  - DST scouting trigger bits (Muon, Electron, JetHT)
  - Parking trigger bits (HH, Muon, Electron)
  - Combined OR expressions
"""


def or_expr(bits):
    """Build a C++ OR expression string from a list of branch names."""
    return "(" + " || ".join(bits) + ")"


# ---------------------------------------------------------------------------
# DST Scouting triggers
# ---------------------------------------------------------------------------

DST_MU_BIT = [
    "DST_PFScouting_DatasetMuon",
    "DST_PFScouting_DoubleMuon",
    "DST_PFScouting_SingleMuon",
]

DST_EL_BITS = [
    "DST_PFScouting_DoubleEG",
    "DST_PFScouting_SinglePhotonEB",
]

DST_JetHT_BITS = [
    "DST_PFScouting_JetHT",
]

DST_MU_expr = or_expr(DST_MU_BIT)
DST_EL_expr = or_expr(DST_EL_BITS)
DST_JetHT_expr = or_expr(DST_JetHT_BITS)


# ---------------------------------------------------------------------------
# Parking HH triggers
# ---------------------------------------------------------------------------

PARKING_HH_BITS = [
    "HLT_PFHT280_QuadPFJet30",
    "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55",
    "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60",
    "HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60",
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40",
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5",
    "HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70",
    "HLT_PFHT400_SixPFJet32",
    "HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50",
    "HLT_PFHT450_SixPFJet36",
    "HLT_PFHT450_SixPFJet36_PNetBTag0p35",
    # TODO: Add VBF parking triggers too
]


# ---------------------------------------------------------------------------
# Parking Muon triggers (2023 menu)
# ---------------------------------------------------------------------------

PARKING_MUON_BITS_2023 = [
    "HLT_Dimuon0_Jpsi3p5_Muon2",
    "HLT_Dimuon0_Jpsi_NoVertexing",
    "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
    "HLT_Dimuon0_LowMass",
    "HLT_Dimuon0_Upsilon_L1_4p5er2p0",
    "HLT_Dimuon0_Upsilon_NoVertexing",
    "HLT_Dimuon10_Upsilon_y1p4",
    "HLT_Dimuon12_Upsilon_y1p4",
    "HLT_Dimuon14_Phi_Barrel_Seagulls",
    "HLT_Dimuon14_PsiPrime",
    "HLT_Dimuon18_PsiPrime",
    "HLT_Dimuon18_PsiPrime_noCorrL1",
    "HLT_Dimuon24_Phi_noCorrL1",
    "HLT_Dimuon24_Upsilon_noCorrL1",
    "HLT_Dimuon25_Jpsi",
    "HLT_Dimuon25_Jpsi_noCorrL1",
    "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",
    "HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon",
    "HLT_DoubleMu3_TkMu_DsTau3Mu",
    "HLT_DoubleMu3_Trk_Tau3mu",
    "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass",
    "HLT_DoubleMu4_3_Bs",
    "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG",
    "HLT_DoubleMu4_3_Jpsi",
    "HLT_DoubleMu4_3_LowMass",
    "HLT_DoubleMu4_3_Photon4_BsToMMG",
    "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
    "HLT_DoubleMu4_JpsiTrk_Bc",
    "HLT_DoubleMu4_Jpsi_Displaced",
    "HLT_DoubleMu4_Jpsi_NoVertexing",
    "HLT_DoubleMu4_LowMass_Displaced",
    "HLT_DoubleMu4_MuMuTrk_Displaced",
    "HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL",
    "HLT_Mu25_TkMu0_Phi",
    "HLT_Mu30_TkMu0_Psi",
    "HLT_Mu30_TkMu0_Upsilon",
    "HLT_Mu7p5_L2Mu2_Jpsi",
    "HLT_Mu7p5_L2Mu2_Upsilon",
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15",
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1",
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15",
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1",
    "HLT_Trimuon5_3p5_2_Upsilon_Muon",
    "HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon",
]


# ---------------------------------------------------------------------------
# Parking Electron triggers (2023 menu)
# ---------------------------------------------------------------------------

PARKING_ELECTRON_BITS_2023 = [
    "HLT_DoubleEle10_eta1p22_mMax6",
    "HLT_DoubleEle6p5_eta1p22_mMax6",
    "HLT_DoubleEle8_eta1p22_mMax6",
    "HLT_SingleEle8",
    "HLT_SingleEle8_SingleEGL1",
]


# ---------------------------------------------------------------------------
# Combined expressions
# ---------------------------------------------------------------------------

PARKING_HH_expr = or_expr(PARKING_HH_BITS)
PARKING_MU_expr = or_expr(PARKING_MUON_BITS_2023)
PARKING_EG_expr = or_expr(PARKING_ELECTRON_BITS_2023)
PARKING_LEPT_expr = f"({PARKING_MU_expr} || {PARKING_EG_expr})"
