# scripts/06_import_bsseq.R
suppressPackageStartupMessages({
  library(yaml)
  library(bsseq)
  library(tidyverse)
  library(parallel)
})

# 1. 설정 로드
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

# 2. 경로 설정
# Python script output: results/methylation
input_dir <- file.path(paths$output_dir, "methylation")
output_obj_dir <- file.path(paths$output_dir, "r_objects")

if(!dir.exists(output_obj_dir)) dir.create(output_obj_dir, recursive = TRUE)

message("[R-DSS] Loading Bismark coverage files...")

# 3. 파일 리스트 확보 (gzip된 coverage 파일)
# Python: *.deduplicated.bam -> Bismark -> *.bismark.cov.gz
cov_files <- list.files(input_dir, pattern = "cov.gz$", full.names = TRUE)

if(length(cov_files) == 0) {
  stop(paste("[Error] No coverage files found in:", input_dir))
}

# 샘플 이름 추출 (파일명 기반)
# 예: sampleA_val_1...cov.gz -> sampleA
sample_names <- gsub("\\..*", "", basename(cov_files))

message(paste("[R-DSS] Found", length(cov_files), "samples:", paste(sample_names, collapse=", ")))

# 4. bsseq 객체 생성
# params$dss$threads가 있다면 사용, 없으면 기본값 4
n_cores <- if(!is.null(params$dss$threads)) params$dss$threads else 4

BSobj <- read.bismark(
  files = cov_files,
  sampleNames = sample_names,
  rmZeroCov = TRUE,
  strandCollapse = TRUE,
  verbose = TRUE
)

# 5. 저장 (Serialization)
saveRDS(BSobj, file = file.path(output_obj_dir, "bsseq_raw.rds"))

message(paste("[R-DSS] BSseq object saved successfully to:", output_obj_dir))