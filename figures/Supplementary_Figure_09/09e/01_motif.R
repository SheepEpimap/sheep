# Plot tissue-specific motif gene expression heatmap
library(ComplexHeatmap)
library(circlize)

# Suppress all startup messages
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

# Set working directory
setwd('/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/Target_gene_cloest/AA_cluster/')
ht_opt$message = FALSE

# Processing note.
# Processing note.
data <- read.table(
  'all_cluster_gene_expression_last.csv', 
  header = TRUE,
  fill = TRUE,
  sep = "",  # Processing note.
  stringsAsFactors = FALSE,
  check.names = FALSE  # Processing note.
) 

# Processing note.
cat("=== Data basic info ===\n")
cat("Total columns:", ncol(data), "\n")
cat("First 10 column names:\n")
print(head(colnames(data), 10))
cat("Last column name (cluster):", colnames(data)[ncol(data)], "\n")

# Processing note.
data$cluster <- data[, ncol(data)]
# Processing note.
data$cluster <- gsub("^\\s+|\\s+$", "", data$cluster)  # Processing note.
data <- data[!is.na(data$cluster) & data$cluster != "", ]  # Processing note.

# Processing note.
cat("\n=== Verify cluster column ===\n")
cat("First 10 cluster values:\n")
print(head(data$cluster, 10))
cat("Unique cluster values:\n")
print(unique(data$cluster))
cat("Cluster distribution:\n")
print(table(data$cluster))

# Processing note.
# Processing note.
tissue_data <- data[, 7:(ncol(data)-1)]
# Processing note.
colnames(tissue_data) <- gsub("[\\.\\-\\s]+", "-", colnames(tissue_data))
colnames(tissue_data) <- tolower(colnames(tissue_data))  # Processing note.

# Processing note.
tissue_matrix <- as.matrix(tissue_data)
tissue_matrix_scaled <- t(scale(t(tissue_matrix)))  # Processing note.

# Processing note.
split <- factor(data$cluster, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))

# Processing note.
row_ha <- rowAnnotation(
  TSR = data$cluster,
  col = list(TSR = c(
    "C1"="#36b5f1", "C2"="#f18264", "C3"="#ffc4f1",
    "C4"="#d84b4b", "C5"="#f4d578", "C6"="#daaa6c",
    "C7"="#efd80b", "C8"="#dcd71a", "C9"="#b6d7a9",
    "C10"="#f2c063"
  )),
  annotation_width = unit(1.2, "cm"),
  annotation_legend_param = list(
    TSR = list(
      at = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"),
      labels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"),
      title = "Cluster",
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    )
  )
)

# Processing note.
pdf(
  "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/27_new_all_cluster_gene_expression_last.pdf",
  width = 18,
  height = 10
)

# Processing note.
Heatmap(
  tissue_matrix_scaled, 
  name = "expression",
  left_annotation = row_ha, 
  row_split = split, 
  row_gap = unit(0.1, "mm"),  # Processing note.
  border = TRUE,
  cluster_column_slices = FALSE, 
  row_title = NULL, 
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 17),  # Processing note.
  show_row_names = FALSE, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  col = colorRamp2(c(0, 3), c("white", "red")),
  use_raster = TRUE,  # Processing note.
  heatmap_legend_param = list(title = "Expression (Z-score)"),
  border_gp = gpar(col = "gray50", lwd = 0.5)  # Processing note.
)

dev.off()

# Processing note.
cat("\n=== Run completed successfully ===\n")
cat("Heatmap saved to:\n/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/27_new_all_cluster_gene_expression_last.pdf\n")