#!/usr/bin/env bash
# ==============================================================================
# Script: run_pipeline_6to8.sh
# Purpose: Run R DSS analysis pipeline (Steps 6-8)
#          Import BSseq → DSS Modeling → Visualization
# Input: .cov files from results/methylation/
# Output: DMR results and visualization plots
# ==============================================================================
set -euo pipefail

# 프로젝트 루트 디렉토리로 이동
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "=============================================="
echo " R DSS Analysis Pipeline (Steps 6-8)"
echo " Started: $(date)"
echo "=============================================="

# .cov 파일 존재 확인
COV_DIR="results/methylation"
if [ ! -d "$COV_DIR" ] || [ -z "$(ls -A "$COV_DIR"/*.cov* 2>/dev/null)" ]; then
    echo "[ERROR] No .cov files found in $COV_DIR"
    echo "        Please run 'run_pipeline_1to5.sh' first."
    exit 1
fi

echo "[INFO] Found .cov files in $COV_DIR"
ls -la "$COV_DIR"/*.cov* 2>/dev/null | head -5
echo ""

echo "[STEP 6/8] Import BSseq - Create HDF5-backed BSseq object"
echo "----------------------------------------------"
Rscript scripts/06_import_bsseq.R

echo ""
echo "[STEP 7/8] DSS Modeling - Differential Methylation Analysis"
echo "----------------------------------------------"
Rscript scripts/07_dss_modeling.R

echo ""
echo "[STEP 8/8] Visualization - Generate Plots"
echo "----------------------------------------------"
Rscript scripts/08_visualization.R

echo ""
echo "=============================================="
echo " Pipeline 6-8 COMPLETE!"
echo " Finished: $(date)"
echo ""
echo " Results:"
echo "   - HDF5 BSseq: results/r_objects/dss_h5_store/"
echo "   - DMR Results: results/dss_results/"
echo "   - Plots: results/plots/"
echo "=============================================="
