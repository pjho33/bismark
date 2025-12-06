# ==============================================================================
# Script: 07_dss_modeling.R (HDF5 Version)
# Purpose: Run DSS DML/DMR analysis using HDF5-backed BSseq object
# Input: 'r_objects/dss_h5_store/bsseq.rds' (HDF5-backed)
# Output: DML/DMR results in 'dss_results/'
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(DSS)
  library(bsseq)
  library(data.table)
  library(rhdf5)
  library(HDF5Array)
  library(GenomicRanges)
  library(gtools)
  library(rlang)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
OUT_DIR <- file.path(paths$output_dir, "dss_results")
h5_file <- file.path(H5_DIR, "se.h5")

# DSS Parameters (from config or defaults)
DSS_SMOOTHING <- params$dss$smoothing %||% TRUE
DSS_SPAN <- params$dss$smoothing_span %||% 500
DSS_DELTA <- params$dss$delta %||% 0.1
DSS_P_VAL <- params$dss$p_threshold %||% 0.001
DMR_MIN_LEN <- params$dss$dmr_min_len %||% 50
DMR_MIN_CG <- params$dss$dmr_min_cg %||% 3
DMR_DIS_MERGE <- params$dss$dmr_dis_merge %||% 100

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# ------------------------------------------------------------------------------
# 2. Validation
# ------------------------------------------------------------------------------
logit("=== DSS Analysis Pipeline (HDF5 Mode) ===")

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
# 4. Define Groups
# ------------------------------------------------------------------------------
group1 <- targets$SampleID[targets$Group == "Control"]
group2 <- targets$SampleID[targets$Group == "Case"]

logit(paste("Group 1 (Control):", length(group1), "samples"))
logit(paste("Group 2 (Case):", length(group2), "samples"))

if(length(group1) == 0 || length(group2) == 0) {
  stop("[Error] Groups are not defined correctly. Check sample naming or control_pattern in params.yaml")
}

# ------------------------------------------------------------------------------
# 5. Run DSS Analysis
# ------------------------------------------------------------------------------
logit("=== Running DSS Analysis ===")
logit(paste("Parameters: smoothing=", DSS_SMOOTHING, ", delta=", DSS_DELTA, ", p.threshold=", DSS_P_VAL))

logit("1. Running DML Test (This is the heavy calculation)...")
dml_res <- DMLtest(bs_obj, group1=group1, group2=group2, smoothing=DSS_SMOOTHING)

# 전체 DML 결과 저장
fwrite(as.data.table(dml_res), file.path(OUT_DIR, "All_DMLs.tsv.gz"), sep="\t", compress="gzip")
logit(paste("   -> Total Sites tested:", nrow(dml_res)))

logit("2. Filtering Significant DMLs...")
sig_dmls <- dml_res[!is.na(dml_res$pval) & dml_res$pval < DSS_P_VAL & abs(dml_res$diff) > DSS_DELTA, ]
fwrite(sig_dmls, file.path(OUT_DIR, "Significant_DMLs.tsv"), sep="\t")
logit(paste("   -> Significant DMLs:", nrow(sig_dmls)))

logit("3. Calling DMRs...")
dmrs <- callDMR(dml_res, p.threshold=DSS_P_VAL, delta=DSS_DELTA,
                minlen=DMR_MIN_LEN, minCG=DMR_MIN_CG, dis.merge=DMR_DIS_MERGE)

if(!is.null(dmrs) && is.data.frame(dmrs) && nrow(dmrs) > 0){
  fwrite(dmrs, file.path(OUT_DIR, "Final_DMRs.tsv"), sep="\t")
  logit(paste("   -> DMRs found:", nrow(dmrs)))
  
  # BED 파일 생성 (시각화용)
  bed_df <- as.data.frame(dmrs[, c("chr", "start", "end", "diff.Methy")])
  bed_df$name <- paste0("DMR_", 1:nrow(bed_df))
  bed_df$score <- round(bed_df$diff.Methy * 1000)
  write.table(bed_df[, c(1,2,3,5,6)], file.path(OUT_DIR, "DMRs_visualization.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  logit("   -> BED file created for visualization")
} else {
  logit("   -> No significant DMRs found with current parameters")
}

# RDS 저장 (선택적 - 용량이 클 수 있음)
# saveRDS(dml_res, file.path(OUT_DIR, "dml_test_result.rds"))
# saveRDS(dmrs, file.path(OUT_DIR, "dmrs_final.rds"))

logit("======================================================")
logit(" SUCCESS! DSS Analysis Complete")
logit(paste(" Results saved to:", OUT_DIR))
logit("======================================================")