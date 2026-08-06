#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidytext)   # reorder_within / scale_y_reordered

setwd('/vol2/zhangshiwen/sheep_cor/h3k27ac/summary/GO')

# =========================
# =========================
read_go_topn <- function(file, group_name, n_top = 10) {
  df <- read_tsv(file, show_col_types = FALSE)

  need_cols <- c("Description", "p.adjust", "Count")
  miss_cols <- setdiff(need_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop(paste("file :", paste(miss_cols, collapse = ", "), "\nfile:", file))
  }

  df %>%
    mutate(
      p.adjust = as.numeric(p.adjust),
      Count = as.numeric(Count),
      group = group_name
    ) %>%
    filter(!is.na(p.adjust), !is.na(Count)) %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = n_top) %>%
    mutate(log10FDR = -log10(p.adjust))
}

# =========================
# =========================
top_df <- read_go_topn(
  "GO_ENSEMBL_H3K27ac_output_E5_gene_count_top005_56_gene.txt",
  "Top",
  n_top = 10
)

mid_df <- read_go_topn(
  "GO_ENSEMBL_H3K27ac_output_E5_gene_count_medium_13_18_gene.txt",
  "Middle",
  n_top = 10
)

bot_df <- read_go_topn(
  "GO_ENSEMBL_H3K27ac_output_E5_gene_count_down005_2_gene.txt",
  "Bottom",
  n_top = 10
)

plot_df <- bind_rows(top_df, mid_df, bot_df)

# =========================
# =========================
plot_df <- plot_df %>%
  mutate(
    Description_wrap = str_wrap(Description, width = 38),
    group = factor(group, levels = c("Top", "Middle", "Bottom"))
  ) %>%
  group_by(group) %>%
  arrange(log10FDR, .by_group = TRUE) %>%
  ungroup()

plot_df <- plot_df %>%
  mutate(
    Description_wrap = reorder_within(Description_wrap, log10FDR, group)
  )

# =========================
# =========================
group_cols <- c(
  "Top"    = "#355C7D",   #
  "Middle" = "#C06C2B",   #
  "Bottom" = "#6C7A59"    #
)

# =========================
# =========================
p <- ggplot(plot_df, aes(x = log10FDR, y = Description_wrap)) +
  geom_point(aes(size = Count, color = group), alpha = 0.95) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  scale_y_reordered() +
  scale_color_manual(values = group_cols) +
  scale_size_continuous(range = c(2.2, 7)) +
  labs(
    x = expression(-log[10]("FDR")),
    y = NULL,
    size = "Gene count",
    color = NULL
  ) +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold", color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 6.5, color = "black"),
    axis.title.x = element_text(size = 8, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6.5),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )

# =========================
# =========================
ggsave(
  "GO_top10_facet_dotplot_nature_style.pdf",
  p,
  width = 7.1,
  height = 4.2,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)

ggsave(
  "GO_top10_facet_dotplot_nature_style.png",
  p,
  width = 7.1,
  height = 4.2,
  units = "in",
  dpi = 600
)
