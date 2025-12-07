# ==============================================================================
# Script: 07_dmrseq.R (HDF5 Version)
# Purpose: Run dmrseq DMR analysis using HDF5-backed BSseq object
# Input: 'r_objects/dss_h5_store/bsseq.rds' (HDF5-backed)
# Output: DMR results in 'dmrseq_results/'
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dmrseq)
  library(bsseq)
  library(data.table)
  library(rhdf5)
  library(HDF5Array)
  library(GenomicRanges)
  library(gtools)
  library(rlang)
  library(BiocParallel)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
OUT_DIR <- file.path(paths$output_dir, "dmrseq_results")
h5_file <- file.path(H5_DIR, "se.h5")

# dmrseq Parameters (from config or defaults)
DMRSEQ_CUTOFF <- params$dmrseq$cutoff %||% 0.05
DMRSEQ_MIN_NUMREGION <- params$dmrseq$min_num_region %||% 3
DMRSEQ_BPSPAN <- params$dmrseq$bp_span %||% 1000
DMRSEQ_MAXGAP <- params$dmrseq$max_gap %||% 1000
DMRSEQ_MAXPERMS <- params$dmrseq$max_perms %||% 10
DMRSEQ_THREADS <- params$dmrseq$threads %||% 4

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# ------------------------------------------------------------------------------
# 2. Validation
# ------------------------------------------------------------------------------
logit("=== dmrseq Analysis Pipeline (HDF5 Mode) ===")

if(!file.exists(h5_file)) stop("se.h5 not found! Run 06_import_bsseq.R first.")

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. Load or Rebuild BSseq Object
# ------------------------------------------------------------------------------
bs_file <- file.path(H5_DIR, "bsseq.rds")
targets_file <- file.path(H5_DIR, "targets.rds")

if(file.exists(bs_file) && file.exists(targets_file)) {
  logit("Loading pre-built BSseq object...")
  bs_obj <- readRDS(bs_file)
  targets <- readRDS(targets_file)
} else {
  # targets 정보 재생성 (fallback)
  logit("WARNING: Re-building Universe to map HDF5 data...")
  
  input_dir <- file.path(paths$output_dir, "methylation")
  cov_files <- list.files(input_dir, full.names = TRUE, 
                          pattern = "\\.cov$|\\.cov\\.gz$|\\.txt$|\\.tsv$")
  sample_ids <- gsub("\\..*", "", basename(cov_files))
  control_pattern <- params$dss$control_pattern %||% "[0-9]+N|Control|WT|Normal"
  groups <- ifelse(grepl(control_pattern, sample_ids, ignore.case = TRUE), "Control", "Case")
  targets <- data.table(SampleID = sample_ids, Group = groups, FilePath = cov_files)
  
  # Universe 재생성
  basic_chrs <- c(1:22, "X", "Y", "M", "MT")
  STD_CHRS <- unique(c(paste0("chr", basic_chrs), as.character(basic_chrs)))
  
  pos_list <- list()
  pb <- txtProgressBar(min=0, max=nrow(targets), style=3)
  for(i in 1:nrow(targets)){
    dt <- fread(targets$FilePath[i], select=1:2, header=FALSE, col.names=c("chr", "pos"))
    dt[, chr := as.character(chr)]
    dt <- dt[chr %in% STD_CHRS]
    pos_list[[i]] <- dt
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  logit("Merging and calculating Union...")
  unique_pos <- unique(rbindlist(pos_list))
  rm(pos_list); gc()
  
  logit("Sorting Universe...")
  unique_pos[, chr_rank := match(chr, mixedsort(unique(chr)))]
  setorder(unique_pos, chr_rank, pos)
  
  gr_universe <- GRanges(seqnames=unique_pos$chr, ranges=IRanges(start=unique_pos$pos, width=1))
  logit(paste("Universe Re-built:", length(gr_universe), "sites"))
  
  # HDF5 연결
  HDF5_M <- HDF5Array(h5_file, "M")
  HDF5_Cov <- HDF5Array(h5_file, "Cov")
  
  bs_obj <- BSseq(M=HDF5_M, Cov=HDF5_Cov, gr=gr_universe, sampleNames=targets$SampleID)
  pData(bs_obj)$Group <- targets$Group
}

logit(paste("BSseq object loaded:", length(bs_obj), "CpGs,", ncol(bs_obj), "samples"))

# ------------------------------------------------------------------------------
# 4. Prepare for dmrseq
# ------------------------------------------------------------------------------
# dmrseq는 pData에 condition 정보가 필요
pData(bs_obj)$Condition <- factor(pData(bs_obj)$Group, levels = c("Control", "Case"))

group1_n <- sum(pData(bs_obj)$Condition == "Control")
group2_n <- sum(pData(bs_obj)$Condition == "Case")

logit(paste("Control samples:", group1_n))
logit(paste("Case samples:", group2_n))

if(group1_n == 0 || group2_n == 0) {
  stop("[Error] Groups are not defined correctly. Check sample naming or control_pattern in params.yaml")
}

# ------------------------------------------------------------------------------
# 5. Filter Low Coverage CpGs (Memory Optimization)
# ------------------------------------------------------------------------------
logit("=== Filtering Low Coverage CpGs ===")

# 최소 coverage 필터링 (메모리 절약 + 품질 향상)
MIN_COV <- params$dmrseq$min_coverage %||% 2
MIN_SAMPLES <- params$dmrseq$min_samples %||% ceiling(ncol(bs_obj) * 0.5)

logit(paste("Filtering: min coverage =", MIN_COV, ", min samples =", MIN_SAMPLES))

# Coverage 기반 필터링
cov_mat <- getCoverage(bs_obj, type = "Cov")
keep_loci <- rowSums(cov_mat >= MIN_COV) >= MIN_SAMPLES

logit(paste("Original CpGs:", length(bs_obj)))
logit(paste("After filtering:", sum(keep_loci), "(", round(sum(keep_loci)/length(bs_obj)*100, 1), "%)"))

bs_filtered <- bs_obj[keep_loci, ]
rm(cov_mat, keep_loci); gc()

# ------------------------------------------------------------------------------
# 6. Set up Parallel Processing
# ------------------------------------------------------------------------------
logit(paste("Setting up parallel processing with", DMRSEQ_THREADS, "threads"))

# BiocParallel 설정
if(.Platform$OS.type == "unix") {
  register(MulticoreParam(workers = DMRSEQ_THREADS))
} else {
  register(SnowParam(workers = DMRSEQ_THREADS))
}

# ------------------------------------------------------------------------------
# 7. Run dmrseq Analysis
# ------------------------------------------------------------------------------
logit("=== Running dmrseq Analysis ===")
logit(paste("Parameters: cutoff=", DMRSEQ_CUTOFF, 
            ", minNumRegion=", DMRSEQ_MIN_NUMREGION,
            ", bpSpan=", DMRSEQ_BPSPAN,
            ", maxGap=", DMRSEQ_MAXGAP,
            ", maxPerms=", DMRSEQ_MAXPERMS))

logit("Running dmrseq (this may take a while)...")

# dmrseq 실행
dmrs <- tryCatch({
  dmrseq(bs = bs_filtered,
         testCovariate = "Condition",
         cutoff = DMRSEQ_CUTOFF,
         minNumRegion = DMRSEQ_MIN_NUMREGION,
         bpSpan = DMRSEQ_BPSPAN,
         maxGap = DMRSEQ_MAXGAP,
         maxPerms = DMRSEQ_MAXPERMS,
         verbose = TRUE)
}, error = function(e) {
  logit(paste("ERROR in dmrseq:", e$message))
  return(NULL)
})

# ------------------------------------------------------------------------------
# 8. Save Results
# ------------------------------------------------------------------------------
if(!is.null(dmrs) && length(dmrs) > 0) {
  logit(paste("DMRs found:", length(dmrs)))
  
  # GRanges를 data.frame으로 변환
  dmr_df <- as.data.frame(dmrs)
  
  # TSV 저장
  fwrite(dmr_df, file.path(OUT_DIR, "dmrseq_DMRs.tsv"), sep="\t")
  logit("   -> TSV file saved")
  
  # Significant DMRs (qval < 0.05)
  sig_dmrs <- dmr_df[dmr_df$qval < 0.05, ]
  if(nrow(sig_dmrs) > 0) {
    fwrite(sig_dmrs, file.path(OUT_DIR, "dmrseq_Significant_DMRs.tsv"), sep="\t")
    logit(paste("   -> Significant DMRs (qval < 0.05):", nrow(sig_dmrs)))
  }
  
  # BED 파일 생성 (시각화용)
  bed_df <- data.frame(
    chr = dmr_df$seqnames,
    start = dmr_df$start,
    end = dmr_df$end,
    name = paste0("DMR_", 1:nrow(dmr_df)),
    score = round(dmr_df$stat * 100),
    strand = "."
  )
  write.table(bed_df, file.path(OUT_DIR, "dmrseq_DMRs.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  logit("   -> BED file created for visualization")
  
  # RDS 저장 (GRanges 객체 보존)
  saveRDS(dmrs, file.path(OUT_DIR, "dmrseq_result.rds"))
  logit("   -> RDS file saved")
  
  # Summary 통계
  summary_stats <- data.frame(
    Metric = c("Total DMRs", "Significant (qval<0.05)", "Hyper-methylated", "Hypo-methylated",
               "Mean Width (bp)", "Median Width (bp)", "Mean CpGs per DMR"),
    Value = c(
      nrow(dmr_df),
      sum(dmr_df$qval < 0.05, na.rm=TRUE),
      sum(dmr_df$stat > 0, na.rm=TRUE),
      sum(dmr_df$stat < 0, na.rm=TRUE),
      round(mean(dmr_df$width)),
      round(median(dmr_df$width)),
      round(mean(dmr_df$L))
    )
  )
  fwrite(summary_stats, file.path(OUT_DIR, "dmrseq_summary.tsv"), sep="\t")
  logit("   -> Summary statistics saved")
  
} else {
  logit("   -> No DMRs found with current parameters")
  logit("   -> Try adjusting: cutoff, minNumRegion, or coverage filters")
}

logit("======================================================")
logit(" SUCCESS! dmrseq Analysis Complete")
logit(paste(" Results saved to:", OUT_DIR))
logit("======================================================")
