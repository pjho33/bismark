# ==============================================================================
# Script: 08_visualization.R
# Purpose: Comprehensive visualization of DSS DML/DMR results
# Output: Multiple plots in 'results/plots/'
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(DSS)
  library(bsseq)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(rlang)
})

logit <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep="")

# ------------------------------------------------------------------------------
# 1. Configuration & Data Loading
# ------------------------------------------------------------------------------
paths <- read_yaml("config/paths.yaml")
params <- read_yaml("config/params.yaml")

H5_DIR <- file.path(paths$output_dir, "r_objects", "dss_h5_store")
OUT_DIR <- file.path(paths$output_dir, "dss_results")
plot_dir <- file.path(paths$output_dir, "plots")

if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

logit("=== DSS Visualization Pipeline ===")
logit("Loading data...")

# BSseq 객체 로드
bs_file <- file.path(H5_DIR, "bsseq.rds")
targets_file <- file.path(H5_DIR, "targets.rds")

if(file.exists(bs_file)) {
  BSobj <- readRDS(bs_file)
  targets <- readRDS(targets_file)
} else {
  stop("BSseq object not found! Run 06_import_bsseq.R first.")
}

# DML/DMR 결과 로드
dml_file <- file.path(OUT_DIR, "All_DMLs.tsv.gz")
dmr_file <- file.path(OUT_DIR, "Final_DMRs.tsv")

if(!file.exists(dml_file)) stop("DML results not found! Run 07_dss_modeling.R first.")

dml_res <- fread(dml_file)
dmrs <- if(file.exists(dmr_file)) fread(dmr_file) else NULL

logit(paste("Loaded:", nrow(dml_res), "DMLs"))
if(!is.null(dmrs)) logit(paste("Loaded:", nrow(dmrs), "DMRs"))

# Parameters
DSS_DELTA <- params$dss$delta %||% 0.1
DSS_P_VAL <- params$dss$p_threshold %||% 0.001

# ------------------------------------------------------------------------------
# 2. Volcano Plot (DML)
# ------------------------------------------------------------------------------
logit("1. Creating Volcano Plot...")

dml_plot <- copy(dml_res)
dml_plot[, neg_log_p := -log10(pval)]
dml_plot[, neg_log_p := pmin(neg_log_p, 50)]  # Cap at 50 for visualization
dml_plot[, significance := "Not Significant"]
dml_plot[pval < DSS_P_VAL & diff > DSS_DELTA, significance := "Hyper (Case)"]
dml_plot[pval < DSS_P_VAL & diff < -DSS_DELTA, significance := "Hypo (Case)"]

p_volcano <- ggplot(dml_plot, aes(x = diff, y = neg_log_p, color = significance)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("Hyper (Case)" = "red", "Hypo (Case)" = "blue", "Not Significant" = "grey60")) +
  geom_hline(yintercept = -log10(DSS_P_VAL), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-DSS_DELTA, DSS_DELTA), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: Differentially Methylated Loci",
       x = "Methylation Difference (Case - Control)",
       y = "-log10(p-value)",
       color = "Status") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "01_Volcano_Plot.png"), p_volcano, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 3. Manhattan Plot
# ------------------------------------------------------------------------------
logit("2. Creating Manhattan Plot...")

# 염색체 순서 정의
chr_order <- c(paste0("chr", c(1:22, "X", "Y")), as.character(c(1:22, "X", "Y")))
dml_manhattan <- copy(dml_res)
dml_manhattan <- dml_manhattan[chr %in% chr_order]
dml_manhattan[, chr_num := gsub("chr", "", chr)]
dml_manhattan[, chr_num := factor(chr_num, levels = c(1:22, "X", "Y"))]

# 누적 위치 계산
dml_manhattan <- dml_manhattan[order(chr_num, pos)]
chr_lengths <- dml_manhattan[, .(max_pos = max(pos)), by = chr_num]
chr_lengths[, cumsum_pos := cumsum(as.numeric(max_pos)) - max_pos]
dml_manhattan <- merge(dml_manhattan, chr_lengths[, .(chr_num, cumsum_pos)], by = "chr_num")
dml_manhattan[, bp_cum := pos + cumsum_pos]

# 염색체 중앙 위치 (x축 라벨용)
axis_df <- dml_manhattan[, .(center = mean(bp_cum)), by = chr_num]

# 색상
dml_manhattan[, chr_color := as.numeric(chr_num) %% 2]

p_manhattan <- ggplot(dml_manhattan, aes(x = bp_cum, y = -log10(pval), color = factor(chr_color))) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("0" = "#1f78b4", "1" = "#a6cee3"), guide = "none") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_num) +
  geom_hline(yintercept = -log10(DSS_P_VAL), linetype = "dashed", color = "red") +
  labs(title = "Manhattan Plot: Genome-wide DML Distribution",
       x = "Chromosome", y = "-log10(p-value)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "02_Manhattan_Plot.png"), p_manhattan, width = 14, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 4. PCA/MDS Plot
# ------------------------------------------------------------------------------
logit("3. Creating PCA Plot...")

# Methylation 값 추출 (샘플링하여 메모리 절약)
set.seed(42)
n_sites <- min(50000, length(BSobj))
sample_idx <- sort(sample(1:length(BSobj), n_sites))

M_mat <- as.matrix(getMeth(BSobj[sample_idx, ], type = "raw"))
M_mat[is.na(M_mat)] <- 0.5  # NA를 0.5로 대체

# PCA 수행
pca_res <- prcomp(t(M_mat), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Sample = colnames(M_mat),
  Group = targets$Group
)

var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
  labs(title = "PCA: Sample Clustering by Methylation",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "03_PCA_Plot.png"), p_pca, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Sample Heatmap (Top DMRs)
# ------------------------------------------------------------------------------
if(!is.null(dmrs) && nrow(dmrs) > 0) {
  logit("4. Creating Heatmap...")
  
  # 상위 DMR 선택
  top_n_dmr <- min(50, nrow(dmrs))
  dmrs_sorted <- dmrs[order(-abs(diff.Methy))][1:top_n_dmr]
  
  # DMR 영역의 평균 methylation 계산
  dmr_meth <- matrix(NA, nrow = top_n_dmr, ncol = ncol(BSobj))
  colnames(dmr_meth) <- sampleNames(BSobj)
  rownames(dmr_meth) <- paste0(dmrs_sorted$chr, ":", dmrs_sorted$start, "-", dmrs_sorted$end)
  
  gr_bs <- granges(BSobj)
  
  for(i in 1:top_n_dmr) {
    dmr_gr <- GRanges(seqnames = dmrs_sorted$chr[i],
                      ranges = IRanges(start = dmrs_sorted$start[i], end = dmrs_sorted$end[i]))
    overlaps <- findOverlaps(gr_bs, dmr_gr)
    if(length(overlaps) > 0) {
      idx <- queryHits(overlaps)
      meth_vals <- getMeth(BSobj[idx, ], type = "raw")
      dmr_meth[i, ] <- colMeans(as.matrix(meth_vals), na.rm = TRUE)
    }
  }
  
  # NA 행 제거
  dmr_meth <- dmr_meth[complete.cases(dmr_meth), , drop = FALSE]
  
  if(nrow(dmr_meth) > 1) {
    # Annotation
    anno_col <- data.frame(Group = targets$Group, row.names = targets$SampleID)
    anno_colors <- list(Group = c("Control" = "#4DAF4A", "Case" = "#E41A1C"))
    
    png(file.path(plot_dir, "04_Heatmap_TopDMRs.png"), width = 1200, height = 1000, res = 150)
    pheatmap(dmr_meth,
             scale = "row",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             annotation_col = anno_col,
             annotation_colors = anno_colors,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = paste("Top", nrow(dmr_meth), "DMRs - Sample Clustering"),
             fontsize_row = 6,
             fontsize_col = 8)
    dev.off()
  }
}

# ------------------------------------------------------------------------------
# 6. DMR Length Distribution
# ------------------------------------------------------------------------------
if(!is.null(dmrs) && nrow(dmrs) > 0) {
  logit("5. Creating DMR Length Distribution...")
  
  dmrs[, length := end - start]
  
  p_length <- ggplot(dmrs, aes(x = length)) +
    geom_histogram(bins = 50, fill = "#3288BD", color = "white", alpha = 0.8) +
    scale_x_log10() +
    labs(title = "DMR Length Distribution",
         x = "DMR Length (bp, log scale)",
         y = "Count") +
    theme_bw()
  
  ggsave(file.path(plot_dir, "05_DMR_Length_Distribution.png"), p_length, width = 8, height = 6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 7. Chromosome Bar Plot (DMR count per chromosome)
# ------------------------------------------------------------------------------
if(!is.null(dmrs) && nrow(dmrs) > 0) {
  logit("6. Creating Chromosome Distribution...")
  
  chr_counts <- dmrs[, .N, by = chr]
  chr_counts[, chr_num := gsub("chr", "", chr)]
  chr_counts[, chr_num := factor(chr_num, levels = c(1:22, "X", "Y"))]
  chr_counts <- chr_counts[!is.na(chr_num)]
  
  p_chr <- ggplot(chr_counts, aes(x = chr_num, y = N, fill = N)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#FEE08B", high = "#D53E4F") +
    labs(title = "DMR Distribution by Chromosome",
         x = "Chromosome",
         y = "Number of DMRs") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(plot_dir, "06_DMR_Chromosome_Distribution.png"), p_chr, width = 10, height = 6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 8. Hyper/Hypo Methylation Pie Chart
# ------------------------------------------------------------------------------
if(!is.null(dmrs) && nrow(dmrs) > 0) {
  logit("7. Creating Hyper/Hypo Pie Chart...")
  
  dmrs[, direction := ifelse(diff.Methy > 0, "Hypermethylated", "Hypomethylated")]
  direction_counts <- dmrs[, .N, by = direction]
  direction_counts[, pct := round(N / sum(N) * 100, 1)]
  direction_counts[, label := paste0(direction, "\n(n=", N, ", ", pct, "%)")]
  
  p_pie <- ggplot(direction_counts, aes(x = "", y = N, fill = direction)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = c("Hypermethylated" = "#E41A1C", "Hypomethylated" = "#377EB8")) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = "DMR Direction: Hyper vs Hypo Methylation") +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(file.path(plot_dir, "07_DMR_Direction_PieChart.png"), p_pie, width = 8, height = 8, dpi = 300)
}

# ------------------------------------------------------------------------------
# 9. Top DMR Methylation Tracks (PDF)
# ------------------------------------------------------------------------------
if(!is.null(dmrs) && nrow(dmrs) > 0) {
  logit("8. Creating Top DMR Detail Plots...")
  
  # DML 결과를 data.frame으로 변환 (showOneDMR용)
  dml_df <- as.data.frame(dml_res)
  dmrs_df <- as.data.frame(dmrs)
  dmrs_df <- dmrs_df[order(-abs(dmrs_df$diff.Methy)), ]
  
  top_n <- min(10, nrow(dmrs_df))
  
  pdf(file.path(plot_dir, "08_Top_DMRs_Detail.pdf"), width = 12, height = 6)
  
  for(i in 1:top_n) {
    dmr <- dmrs_df[i, ]
    
    tryCatch({
      showOneDMR(dmrs_df, dml_df, i,
                 ext = 500,
                 main = paste0("Rank ", i, ": ", dmr$chr, ":", dmr$start, "-", dmr$end,
                              " (diff=", round(dmr$diff.Methy, 3), ")"))
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error plotting DMR", i))
    })
  }
  
  dev.off()
}

# ------------------------------------------------------------------------------
# 10. Summary Statistics
# ------------------------------------------------------------------------------
logit("9. Creating Summary Report...")

summary_stats <- data.frame(
  Metric = c(
    "Total CpG Sites Tested",
    "Significant DMLs (p < threshold)",
    "Hypermethylated DMLs",
    "Hypomethylated DMLs",
    "Total DMRs",
    "Hypermethylated DMRs",
    "Hypomethylated DMRs",
    "Median DMR Length (bp)",
    "Control Samples",
    "Case Samples"
  ),
  Value = c(
    nrow(dml_res),
    nrow(dml_res[pval < DSS_P_VAL & abs(diff) > DSS_DELTA]),
    nrow(dml_res[pval < DSS_P_VAL & diff > DSS_DELTA]),
    nrow(dml_res[pval < DSS_P_VAL & diff < -DSS_DELTA]),
    if(!is.null(dmrs)) nrow(dmrs) else 0,
    if(!is.null(dmrs)) nrow(dmrs[diff.Methy > 0]) else 0,
    if(!is.null(dmrs)) nrow(dmrs[diff.Methy < 0]) else 0,
    if(!is.null(dmrs) && nrow(dmrs) > 0) median(dmrs$end - dmrs$start) else NA,
    sum(targets$Group == "Control"),
    sum(targets$Group == "Case")
  )
)

fwrite(summary_stats, file.path(plot_dir, "00_Summary_Statistics.csv"))

logit("======================================================")
logit(" SUCCESS! All visualizations saved to:")
logit(paste(" ", plot_dir))
logit("")
logit(" Generated files:")
logit("   00_Summary_Statistics.csv")
logit("   01_Volcano_Plot.png")
logit("   02_Manhattan_Plot.png")
logit("   03_PCA_Plot.png")
logit("   04_Heatmap_TopDMRs.png")
logit("   05_DMR_Length_Distribution.png")
logit("   06_DMR_Chromosome_Distribution.png")
logit("   07_DMR_Direction_PieChart.png")
logit("   08_Top_DMRs_Detail.pdf")
logit("======================================================")