#!/bin/bash
# ============================================================================
# Run GenXSecAnalyzer on TT and DYJetsNLO samples
#
# Prerequisites:
#   1. cmsenv must be set (source /cvmfs/... && cd CMSSW_*/src && cmsenv)
#   2. Grid proxy: voms-proxy-init -voms cms
#   3. genXsec_cfg.py must exist in this directory (auto-downloaded below)
#
# Usage:
#   cd scripts/
#   source run_xsec.sh          # runs all datasets
#   source run_xsec.sh 3        # runs only dataset #3 (0-indexed)
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANA_DIR="$(dirname "$SCRIPT_DIR")"
cd "$SCRIPT_DIR"

DATASET_FILE="${ANA_DIR}/output/xsec/datasets_xsec.txt"
LOG_DIR="${ANA_DIR}/logs/xsec"
NEVENTS=1000000   # 1M events -> ~0.1% precision, ~3 min per dataset
DASGOCLIENT="/cvmfs/cms.cern.ch/common/dasgoclient"
XROOTD_PREFIX="root://cms-xrd-global.cern.ch/"

# Download genXsec_cfg.py if not present
if [ ! -f genXsec_cfg.py ]; then
    echo ">>> Downloading genXsec_cfg.py ..."
    curl -sL https://raw.githubusercontent.com/cms-sw/genproductions/master/Utilities/calculateXSectionAndFilterEfficiency/genXsec_cfg.py -o genXsec_cfg.py
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to download genXsec_cfg.py"
        return 1 2>/dev/null || exit 1
    fi
fi

mkdir -p "$LOG_DIR"

# Check prerequisites
if ! command -v cmsRun &>/dev/null; then
    echo "ERROR: cmsRun not found. Did you run 'cmsenv'?"
    return 1 2>/dev/null || exit 1
fi

if ! voms-proxy-info -exists &>/dev/null; then
    echo "WARNING: No valid grid proxy found. Run: voms-proxy-init -voms cms"
fi

# Read datasets into array
mapfile -t DATASETS < <(grep -v '^#' "$DATASET_FILE" | grep -v '^$')

echo "============================================"
echo " GenXSecAnalyzer - ${#DATASETS[@]} datasets"
echo " Events per dataset: $NEVENTS"
echo "============================================"
echo ""

run_one() {
    local idx=$1
    local dataset="${DATASETS[$idx]}"
    # Short name for the log file (primary dataset name)
    local short_name
    short_name=$(echo "$dataset" | cut -d'/' -f2)
    local logfile="${LOG_DIR}/${short_name}.log"

    echo ">>> [$((idx+1))/${#DATASETS[@]}] $short_name"
    echo "    Dataset: $dataset"

    # Query DAS for one file
    local das_file
    das_file=$($DASGOCLIENT -query "file dataset=${dataset}" -limit 1 2>/dev/null)
    if [ -z "$das_file" ]; then
        echo "    ERROR: No files found via DAS for $dataset"
        echo "ERROR: No files found via DAS" > "$logfile"
        return 1
    fi

    local input_file="${XROOTD_PREFIX}${das_file}"
    echo "    Input: $input_file"
    echo "    Log:   $logfile"

    cmsRun genXsec_cfg.py \
        inputFiles="$input_file" \
        maxEvents=$NEVENTS \
        &> "$logfile"

    echo ""
}

# Run a single dataset or all
if [ -n "$1" ]; then
    run_one "$1"
else
    for i in "${!DATASETS[@]}"; do
        run_one "$i"
    done

    echo "============================================"
    echo " Summary"
    echo "============================================"
    for logfile in "$LOG_DIR"/*.log; do
        local_name=$(basename "$logfile" .log)
        result=$(grep "After filter: final cross section" "$logfile" 2>/dev/null)
        if [ -n "$result" ]; then
            printf "  %-70s %s\n" "$local_name" "$result"
        else
            printf "  %-70s %s\n" "$local_name" "FAILED (check log)"
        fi
    done
fi
