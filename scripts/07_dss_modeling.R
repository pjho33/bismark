# scripts/07_dss_modeling.R
suppressPackageStartupMessages({
  library(yaml)
  library(DSS)
  library(bsseq)
  library(rlang)
})

paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

obj_dir <- file.path(paths$output_dir, "r_objects")
bs_file <- file.path(obj_dir, "bsseq_raw.rds")

if(!file.exists(bs_file)) stop("Run 06_import_bsseq.R first!")

message("[R-DSS] Loading BSseq object...")
BSobj <- readRDS(bs_file)

# ---------------------------------------------------------
# [중요] 실험 디자인 설정 (Design Matrix)
# 선생님의 샘플 명명 규칙에 따라 수정이 필요한 부분입니다.
# 예: 파일명이 'Tumor_01', 'Normal_01' 식이라면 grepl로 자동화 가능
# ---------------------------------------------------------
sample_names <- sampleNames(BSobj)

# (예시) 'Control'이라는 단어가 들어가면 Group1, 아니면 Group2
group1 <- sample_names[grepl("Control|WT", sample_names, ignore.case = TRUE)]
group2 <- sample_names[!sample_names %in% group1]

message(paste("[R-DSS] Group 1 (Control):", paste(group1, collapse=", ")))
message(paste("[R-DSS] Group 2 (Case):   ", paste(group2, collapse=", ")))

if(length(group1) == 0 || length(group2) == 0) {
  stop("[Error] Groups are not defined correctly. Check sample naming logic.")
}

# 1. DML Test (Smoothing 포함)
# Smoothing은 메모리를 많이 먹으므로 필요시 파라미터 조정
message("[R-DSS] Performing DML test (this may take a while)...")
dmlTest_res <- DMLtest(
  BSobj, 
  group1 = group1, 
  group2 = group2, 
  smoothing = params$dss$smoothing %||% TRUE,
  smoothing.span = params$dss$smoothing_span %||% 500
)

# 2. Call DMRs
delta_val <- params$dss$delta %||% 0.1
p_val <- params$dss$p_threshold %||% 0.001

message(paste("[R-DSS] Calling DMRs (Delta:", delta_val, ", P-val:", p_val, ")..."))
dmrs <- callDMR(dmlTest_res, delta = delta_val, p.threshold = p_val)

# 3. 결과 저장
saveRDS(dmlTest_res, file = file.path(obj_dir, "dml_test_result.rds"))
saveRDS(dmrs, file = file.path(obj_dir, "dmrs_final.rds"))

# 간단한 CSV 내보내기
output_csv_dir <- file.path(paths$output_dir, "dss_tables")
if(!dir.exists(output_csv_dir)) dir.create(output_csv_dir)
write.csv(dmrs, file.path(output_csv_dir, "DMRs.csv"), row.names = FALSE)

message("[R-DSS] Modeling finished. Results saved.")