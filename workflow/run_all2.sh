#!/usr/bin/env bash
set -euo pipefail

# 필요한 conda env 활성화 (이름은 각자 맞게 바꿔야 함)
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate bismark_env

echo "[RUN_ALL] Step 1: FastQC"
python scripts/01_qc.py

echo "[RUN_ALL] Step 2: Trim Galore"
python scripts/02_trim.py

echo "[RUN_ALL] Step 3: Bismark Align"
python scripts/03_align.py

echo "[RUN_ALL] Step 4: Deduplicate"
python scripts/04_dedup.py

echo "[RUN_ALL] Step 5: Methylation Extractor"
python scripts/05_methylation.py

echo "[RUN_ALL] DONE"

# ... python pipeline ends ...

echo "Running R DSS Analysis..."
Rscript scripts/06_import_bsseq.R
Rscript scripts/07_dss_modeling.R
Rscript scripts/08_visualization.R
