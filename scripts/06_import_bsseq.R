# ==============================================================================
# Script: 06_import_bsseq.R (HDF5 Version)
# Purpose: Create HDF5-backed BSseq object for memory-efficient DSS analysis
# Output: 'r_objects/dss_h5_store/se.h5' and 'r_objects/dss_h5_store/bsseq.rds'
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(bsseq)
  library(data.table)
  library(rhdf5)
  library(HDF5Array)
  library(GenomicRanges)
  library(gtools)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

# 입력/출력 경로 설정
input_dir <- file.path(paths$output_dir, "methylation")
H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
OUT_LOG <- file.path(paths$output_dir, "r_objects", "data_loading_log.csv")

# [중요] 표준 염색체 필터 (5.4억개 방지용)
basic_chrs <- c(1:22, "X", "Y", "M", "MT")
STD_CHRS <- unique(c(paste0("chr", basic_chrs), as.character(basic_chrs)))

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# 안전장치: 스크립트 종료 시 HDF5 핸들 닫기
on.exit({ try(H5close(), silent = TRUE) }, add = TRUE)

# ------------------------------------------------------------------------------
# 2. Setup & Validation
# ------------------------------------------------------------------------------
logit("=== [Step 1] Initializing HDF5 Data Pipeline ===")

if(!dir.exists(input_dir)) stop(paste("Directory not found:", input_dir))

# 기존 HDF5가 있다면 삭제하고 새로 시작 (Clean Start)
if(dir.exists(H5_DIR)) {
  logit("Removing existing HDF5 store to ensure fresh build...")
  unlink(H5_DIR, recursive = TRUE)
}
dir.create(H5_DIR, recursive = TRUE)

# 파일 리스트 확보 (Bismark coverage 파일)
cov_files <- list.files(input_dir, full.names = TRUE, 
                        pattern = "\\.cov$|\\.cov\\.gz$|\\.txt$|\\.tsv$")
if(length(cov_files) == 0) stop(paste("No coverage files found in:", input_dir))

# 샘플 이름 추출 및 그룹 할당
sample_ids <- gsub("\\..*", "", basename(cov_files))

# 그룹 패턴 (params에서 읽거나 기본값 사용)
control_pattern <- params$dss$control_pattern %||% "[0-9]+N|Control|WT|Normal"
groups <- ifelse(grepl(control_pattern, sample_ids, ignore.case = TRUE), "Control", "Case")

targets <- data.table(SampleID = sample_ids, Group = groups, FilePath = cov_files)

logit(paste("Detected Samples:", nrow(targets)))
print(table(targets$Group))

# ------------------------------------------------------------------------------
# 3. Build CpG Universe (Map Phase)
# ------------------------------------------------------------------------------
logit("=== [Step 2] Building Standard CpG Universe ===")
logit("Scanning files for standard chromosomes (chr1-22, X, Y)...")

pos_list <- list()
pb <- txtProgressBar(min=0, max=nrow(targets), style=3)

for(i in 1:nrow(targets)){
  # 1,2번 컬럼(좌표)만 읽기
  dt <- fread(targets$FilePath[i], select=1:2, header=FALSE, col.names=c("chr", "pos"))
  
  # 필터링 (Junk Data 제거)
  dt[, chr := as.character(chr)]
  dt <- dt[chr %in% STD_CHRS]
  
  pos_list[[i]] <- dt
  setTxtProgressBar(pb, i)
}
close(pb)

logit("Merging and Sorting Universe...")
# 유니크 좌표 추출
unique_pos <- unique(rbindlist(pos_list))
rm(pos_list); gc()

# 정렬 (Natural Sort: 1, 2, ..., 10)
unique_pos[, chr_rank := match(chr, mixedsort(unique(chr)))]
setorder(unique_pos, chr_rank, pos)

# 인덱스 부여
unique_pos[, idx := 1:.N]
setkey(unique_pos, chr, pos)

n_cpgs <- nrow(unique_pos)
logit(paste("Final Cleaned CpGs:", format(n_cpgs, big.mark=",")))

# ------------------------------------------------------------------------------
# 4. Initialize HDF5 Container
# ------------------------------------------------------------------------------
logit("=== [Step 3] Creating HDF5 Container ===")

h5_file <- file.path(H5_DIR, "se.h5")
h5createFile(h5_file)

# Chunk 설정 (속도 최적화)
chunk_dim <- c(min(n_cpgs, 10000), 1)

logit("Allocating Disk Space...")
h5createDataset(h5_file, "M", c(n_cpgs, nrow(targets)), storage.mode="integer", chunk=chunk_dim, level=1)
h5createDataset(h5_file, "Cov", c(n_cpgs, nrow(targets)), storage.mode="integer", chunk=chunk_dim, level=1)

# ------------------------------------------------------------------------------
# 5. Fill Data (Load Phase)
# ------------------------------------------------------------------------------
logit("=== [Step 4] Filling Data (Sequential Write) ===")

stats_list <- vector("list", nrow(targets))

for(i in 1:nrow(targets)){
  st_time <- Sys.time()
  
  # 진행 상황 출력
  cat(sprintf("[%d/%d] %s ... ", i, nrow(targets), targets$SampleID[i]))
  
  # 데이터 읽기
  raw <- fread(targets$FilePath[i], select=c(1,2,5,6), header=FALSE, 
               col.names=c("chr", "pos", "n_meth", "n_unmeth"))
  raw[, chr := as.character(chr)]
  raw <- raw[chr %in% STD_CHRS]
  
  # 유니버스 매칭
  merged <- unique_pos[raw, on=.(chr, pos), nomatch=NULL]
  
  # 데이터 준비 (M, Cov)
  merged[, n_cov := n_meth + n_unmeth]
  setorder(merged, idx)
  
  # HDF5 쓰기
  if(nrow(merged) > 0){
    h5write(as.integer(merged$n_meth), h5_file, "M", index=list(merged$idx, i))
    h5write(as.integer(merged$n_cov),  h5_file, "Cov", index=list(merged$idx, i))
  }
  
  # 통계 기록
  dur <- round(as.numeric(difftime(Sys.time(), st_time, units="secs")), 1)
  stats_list[[i]] <- data.frame(Sample=targets$SampleID[i], TotalRaw=nrow(raw), 
                                 Matched=nrow(merged), TimeSec=dur)
  
  cat(sprintf("Done! (%.1f sec)\n", dur))
  rm(raw, merged); gc(verbose=FALSE)
}

# 로딩 통계 저장
write.csv(do.call(rbind, stats_list), OUT_LOG, row.names=FALSE)
logit("Data Loading Complete.")

# ------------------------------------------------------------------------------
# 6. Finalize & Save Object
# ------------------------------------------------------------------------------
logit("=== [Step 5] Saving Final BSseq Object ===")

# Universe GRanges 생성
gr_universe <- GRanges(seqnames=unique_pos$chr, ranges=IRanges(start=unique_pos$pos, width=1))

# HDF5Array 연결
HDF5_M   <- HDF5Array(h5_file, "M")
HDF5_Cov <- HDF5Array(h5_file, "Cov")

# BSseq 객체 생성 (HDF5-backed, 메모리 효율적)
bs_obj <- BSseq(M=HDF5_M, Cov=HDF5_Cov, gr=gr_universe, sampleNames=targets$SampleID)
pData(bs_obj)$Group <- targets$Group

# 최종 저장 (껍데기만 저장하므로 용량이 작음)
saveRDS(bs_obj, file.path(H5_DIR, "bsseq.rds"))

# targets 정보도 저장 (07번 스크립트에서 사용)
saveRDS(targets, file.path(H5_DIR, "targets.rds"))

logit("======================================================")
logit(" SUCCESS! HDF5 Database created at: ", H5_DIR)
logit(paste(" Total CpGs:", format(n_cpgs, big.mark=",")))
logit(paste(" Total Samples:", nrow(targets)))
logit(" You can now run 07_dss_modeling.R for analysis.")
logit("======================================================")