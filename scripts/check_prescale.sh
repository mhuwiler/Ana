#!/bin/bash
# Check trigger prescales for a given run and HLT path.
# Prescale=1 means unprescaled (good). Prescale>1 means only 1/N events kept.
#
# Usage:
#   ./check_prescale.sh <hltpath> <run_number>
#
# Examples:
#   ./check_prescale.sh "DST_PFScouting_JetHT*" 379416
#   ./check_prescale.sh "HLT_PFHT280_QuadPFJet30*" 379416
#
# Prerequisites:
#   export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda310/bin:$PATH
#   pip install --user brilws

set -euo pipefail

HLT_PATH="${1:?Usage: $0 <hltpath> <run_number>}"
RUN="${2:?Usage: $0 <hltpath> <run_number>}"

export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda310/bin:$PATH

echo "HLT path: $HLT_PATH"
echo "Run:      $RUN"
echo ""

brilcalc trg --prescale -c web -r "$RUN" --hltpath "$HLT_PATH"
