#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

output_dir <- "complexheatmap_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

heatmap_data <- read.delim("heatmap_data.tsv", check.names = FALSE)
bottom_bar_data <- read.delim("bottom_bar_data.tsv")
right_stacked_data <- read.delim("right_stacked_data.tsv")
tissue_colors_df <- read.delim("tissue_colors.tsv")
system_groups <- read.delim("system_groups.tsv", check.names = FALSE)

tissue_colors <- setNames(tissue_colors_df$color, tissue_colors_df$tissue)

heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$cluster
heatmap_matrix_log <- log2(heatmap_matrix + 1)

motif_tissue_count <- right_stacked_data %>%
  group_by(cluster) %>%
  summarise(
    n_tissues = n_distinct(tissue[seqlets > 0]),
    .groups = "drop"
  )

right_annotation_data <- motif_tissue_count$n_tissues
names(right_annotation_data) <- motif_tissue_count$cluster

bottom_annotation_data <- bottom_bar_data$n_clusters
names(bottom_annotation_data) <- bottom_bar_data$tissue

col_fun <- colorRamp2(
  breaks = c(0, quantile(heatmap_matrix_log, 0.95, na.rm = TRUE)),
  colors = c("white", "#A00000")
)

column_ha <- HeatmapAnnotation(
  Tissue = anno_simple(
    colnames(heatmap_matrix_log),
    col = tissue_colors,
    height = unit(0.5, "cm")
  ),
  `Motif Count` = anno_barplot(
    bottom_annotation_data[colnames(heatmap_matrix_log)],
    bar_width = 0.8,
    gp = gpar(fill = "black", col = "white"),
    axis_param = list(side = "left")
  ),
  annotation_name_side = "left"
)

row_ha <- rowAnnotation(
  `Tissue Count` = anno_barplot(
    right_annotation_data[rownames(heatmap_matrix_log)],
    bar_width = 0.8,
    gp = gpar(fill = "black"),
    axis_param = list(side = "bottom")
  ),
  width = unit(2, "cm")
)

ht <- Heatmap(
  heatmap_matrix_log,
  name = "log2(seqlets+1)",
  col = col_fun,
  na_col = "white",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 0.5),
  border = TRUE,
  top_annotation = column_ha,
  right_annotation = row_ha,
  heatmap_legend_param = list(
    title = "log2(seqlets+1)",
    title_position = "leftcenter-rot",
    legend_height = unit(4, "cm"),
    legend_width = unit(0.5, "cm")
  )
)

pdf(file.path(output_dir, "complexheatmap.pdf"), width = 12, height = 10)
draw(ht, heatmap_legend_side = "right")
dev.off()

png(
  file.path(output_dir, "complexheatmap.png"),
  width = 1200,
  height = 1000,
  res = 150
)
draw(ht, heatmap_legend_side = "right")
dev.off()

ht_draw <- draw(ht)
clustered_row_order <- row_order(ht_draw)
clustered_column_order <- column_order(ht_draw)
sorted_rows <- rownames(heatmap_matrix_log)[clustered_row_order]
sorted_columns <- colnames(heatmap_matrix_log)[clustered_column_order]

row_order_df <- data.frame(
  motif = sorted_rows,
  order = seq_along(sorted_rows)
)
write.table(
  row_order_df,
  file.path(output_dir, "clustered_motif_order.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

column_order_df <- data.frame(
  tissue = sorted_columns,
  order = seq_along(sorted_columns)
)
write.table(
  column_order_df,
  file.path(output_dir, "clustered_tissue_order.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

right_bar_data_sorted <- data.frame(
  cluster = factor(sorted_rows, levels = sorted_rows),
  n_tissues = right_annotation_data[sorted_rows]
)

p_right_bar <- ggplot(
  right_bar_data_sorted,
  aes(x = n_tissues, y = cluster)
) +
  geom_col(fill = "black", width = 0.8) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

bottom_bar_data_sorted <- data.frame(
  tissue = factor(sorted_columns, levels = sorted_columns),
  n_clusters = bottom_annotation_data[sorted_columns]
)

p_bottom <- ggplot(
  bottom_bar_data_sorted,
  aes(x = tissue, y = n_clusters)
) +
  geom_col(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = function(x) pretty(x, n = 5)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    axis.text.y = element_text(size = 7),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  labs(y = "# unique motifs")

ggsave(
  file.path(output_dir, "clustered_right_bar.pdf"),
  p_right_bar,
  width = 3,
  height = 10,
  device = "pdf"
)
ggsave(
  file.path(output_dir, "clustered_right_bar.png"),
  p_right_bar,
  width = 3,
  height = 10,
  dpi = 300,
  device = "png"
)
ggsave(
  file.path(output_dir, "clustered_bottom_bar.pdf"),
  p_bottom,
  width = 12,
  height = 3,
  device = "pdf"
)
ggsave(
  file.path(output_dir, "clustered_bottom_bar.png"),
  p_bottom,
  width = 12,
  height = 3,
  dpi = 300,
  device = "png"
)

message("ComplexHeatmap analysis completed.")
