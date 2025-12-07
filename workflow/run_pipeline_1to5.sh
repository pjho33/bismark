#!/usr/bin/env bash
# ==============================================================================
# Script: run_pipeline_1to5.sh
# Purpose: Run Bismark preprocessing pipeline (Steps 1-5)
#          QC → Trim → Align → Dedup → Methylation Extraction
# Output: .cov files in results/methylation/
# ==============================================================================
set -euo pipefail

# 프로젝트 루트 디렉토리로 이동
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "=============================================="
echo " Bismark Preprocessing Pipeline (Steps 1-5)"
echo " Started: $(date)"
echo "=============================================="

echo ""
echo "[STEP 1/5] FastQC - Quality Control"
echo "----------------------------------------------"
python scripts/01_qc.py

echo ""
echo "[STEP 2/5] Trim Galore - Adapter Trimming"
echo "----------------------------------------------"
python scripts/02_trim.py

echo ""
echo "[STEP 3/5] Bismark Align - Read Alignment"
echo "----------------------------------------------"
python scripts/03_align.py

echo ""
echo "[STEP 4/5] Deduplicate - Remove PCR Duplicates"
echo "----------------------------------------------"
python scripts/04_dedup.py

echo ""
echo "[STEP 5/5] Methylation Extractor - Generate .cov files"
echo "----------------------------------------------"
python scripts/05_methylation.py

echo ""
echo "=============================================="
echo " Pipeline 1-5 COMPLETE!"
echo " Finished: $(date)"
echo " Output: results/methylation/*.cov.gz"
echo ""
echo " Next: Run 'run_pipeline_6to8.sh' for R analysis"
echo "=============================================="
