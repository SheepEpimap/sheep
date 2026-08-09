#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(cluster)
  library(grid)
})

# =========================
# Processing note.
# =========================
setwd("/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/AA_one_count/")

data <- read.csv("all_super_enhancer_Gs_one_count.csv", sep = "\t", header = TRUE, check.names = FALSE)
df2  <- as.matrix(data[, 4:46])
df2  <- apply(df2, 2, as.numeric) 
rownames(df2) <- NULL

# =========================
# Processing note.
# =========================
set.seed(1123)
pa <- pam(df2, k = 10)

cluster_fac <- factor(
  paste0("C", pa$clustering),
  levels = paste0("C", 1:10)
)

# =========================
# Processing note.
# =========================
df2 <- apply(df2, 2, as.numeric) 
rng <- range(df2, na.rm = TRUE)

high_col <- "#CB4B56"

pal <- c(
  "#FFFFFF", 
  "#F2E8AB",  
  "#E7C252",  
  "#E18647",  
  high_col    
)

breaks <- seq(rng[1], rng[2], length.out = length(pal))
mycol  <- circlize::colorRamp2(breaks, pal)


# =========================
# Processing note.
# =========================
pdf("23_all_super_enhancer_Gs_one_count_heatmap.pdf", width = 8, height = 15)

Heatmap(
  df2,
  name = "SE_count",
  col  = mycol,

  show_row_names    = FALSE,
  show_column_names = TRUE,

  split = cluster_fac,  # Processing note.
  cluster_rows = TRUE,  # Processing note.
  cluster_row_slices = FALSE,  # Processing note.
  row_dend_reorder = FALSE,  # Processing note.

  cluster_columns = TRUE,

  # Processing note.
  rect_gp = grid::gpar(col = NA)  # Processing note.
)

dev.off()
