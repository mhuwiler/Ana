#!/bin/bash
# Compute integrated luminosity for a given HLT trigger path.
#
# Usage:
#   ./run_brilcalc.sh <hltpath> [golden_json]
#
# Examples:
#   ./run_brilcalc.sh "DST_PFScouting_JetHT_v*"
#   ./run_brilcalc.sh "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55_v*"
#   ./run_brilcalc.sh "DST_PFScouting_JetHT_v*" /path/to/golden.json
#
# Prerequisites:
#   export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda310/bin:$PATH
#   pip install --user brilws

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANA_DIR="$(dirname "$SCRIPT_DIR")"

HLT_PATH="${1:?Usage: $0 <hltpath> [golden_json]}"
GOLDEN_JSON="${2:-${ANA_DIR}/data/Cert_Collisions2024_378981_386951_Golden.json}"
OUTPUT_DIR="${ANA_DIR}/output/brilcalc"

mkdir -p "$OUTPUT_DIR"

export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda310/bin:$PATH

# Derive output filename from HLT path (strip version glob)
OUT_NAME=$(echo "$HLT_PATH" | sed 's/_v\*$//' | sed 's/_v[0-9]*$//')
OUTFILE="${OUTPUT_DIR}/brilcalc_${OUT_NAME}.txt"

echo "HLT path:    $HLT_PATH"
echo "Golden JSON: $GOLDEN_JSON"
echo "Output:      $OUTFILE"
echo ""

brilcalc lumi \
  -b "STABLE BEAMS" \
  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json \
  -u /fb \
  --hltpath "$HLT_PATH" \
  -i "$GOLDEN_JSON" \
  -c web \
  | tee "$OUTFILE"
