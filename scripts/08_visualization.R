# scripts/08_visualization.R
suppressPackageStartupMessages({
  library(yaml)
  library(DSS)
  library(bsseq)
  library(ggplot2)
})

paths <- read_yaml("config/paths.yaml")
obj_dir <- file.path(paths$output_dir, "r_objects")
plot_dir <- file.path(paths$output_dir, "plots")

if(!dir.exists(plot_dir)) dir.create(plot_dir)

message("[R-DSS] Loading data for visualization...")
BSobj <- readRDS(file.path(obj_dir, "bsseq_raw.rds"))
dmrs <- readRDS(file.path(obj_dir, "dmrs_final.rds"))
dmlTest_res <- readRDS(file.path(obj_dir, "dml_test_result.rds"))

# DMR 정렬 (Area 기준 내림차순 등)
dmrs_sorted <- dmrs[order(dmrs$areaStat, decreasing = TRUE), ]

# 상위 10개 DMR 시각화
top_n <- 10
message(paste("[R-DSS] Generating plots for top", top_n, "DMRs..."))

pdf(file.path(plot_dir, "Top_DMRs.pdf"), width = 10, height = 6)

for(i in 1:min(nrow(dmrs_sorted), top_n)) {
  dmr <- dmrs_sorted[i,]
  
  # showOneDMR은 기본 R plot을 사용
  showOneDMR(dmrs_sorted, dmlTest_res, i, 
             ext = 500, # DMR 앞뒤로 500bp 확장해서 보여줌
             main = paste0("Rank ", i, ": ", dmr$chr, ":", dmr$start, "-", dmr$end))
}

dev.off()

message(paste("[R-DSS] Visualization saved to:", file.path(plot_dir, "Top_DMRs.pdf")))