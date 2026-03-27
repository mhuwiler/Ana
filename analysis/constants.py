"""Physics constants and branch name prefixes for the scouting analysis."""

# ── Branch prefixes ──
AK4 = "ScoutingPFJetRecluster2"
AK8 = "ScoutingFatPFJetRecluster"
AK4_UPART = f"{AK4}_scoutUParT"
AK8_SGP = f"{AK8}_scoutGlobalParT"

# ── UParT probability names (AK4 BvsAll / TauVsAll denominators) ──
UPART_PROB_NAMES = ["probb", "probc", "probg", "probuds", "problepb", "probtaup", "probtaum"]

# ── AK8 ScoutGlobalParT probability names (XvsAll discriminants) ──
AK8_SGP_PROB_NAMES = [
    "Xbb", "Xbc", "Xbs", "Xcc", "Xcs", "Xss", "Xud", "Xgg", "Xqq",
    "Xtauhtauh", "Xtauhtaum", "Xtauhtaue", "QCD",
]

# ── Physics thresholds ──
AK8_PT_MIN = 150.0
AK8_DR_MATCH = 0.8
JET_PT_MIN = 20.0
N_LEAD_JETS = 4
STACK_SIZE_MB = 64
