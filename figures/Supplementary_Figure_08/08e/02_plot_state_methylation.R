#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

############################################
############################################
INFILE <- "/vol2/zhangshiwen/rrbs/wgbs_bismark/state_methylation_by_sample_tissue/all_samples_sample_state_weighted_methylation.tsv"
OUTPDF <- "/vol2/zhangshiwen/rrbs/wgbs_bismark/state_methylation_by_sample_tissue/plot/all_samples_state_weighted_methylation_final.pdf"
OUTTSV <- "/vol2/zhangshiwen/rrbs/wgbs_bismark/state_methylation_by_sample_tissue/plot/all_samples_state_weighted_methylation_final.summary.tsv"

############################################
############################################
data <- read.table(
  INFILE,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

data$weighted_meth_pct <- suppressWarnings(as.numeric(data$weighted_meth_pct))

############################################
############################################
data <- data %>%
  filter(
    !is.na(weighted_meth_pct),
    state != "E11",
    weighted_meth_pct < 99.5
  )

############################################
############################################
state_map <- c(
  "E1"  = "1 TssA",
  "E2"  = "2 TssFlnk",
  "E3"  = "3 TSSWk",
  "E4"  = "4 TssBiv",
  "E5"  = "5 EnhA",
  "E6"  = "6 EnhAMe",
  "E7"  = "7 EnhAHet",
  "E8"  = "8 EnhPois",
  "E9"  = "9 Repr",
  "E10" = "10 QuiW"
)

state_colors <- c(
  "1 TssA"    = "#FF0000",
  "2 TssFlnk" = "#FF2C79",
  "3 TSSWk"   = "#F08080",
  "4 TssBiv"  = "#fce0e0",
  "5 EnhA"    = "#FFFF00",
  "6 EnhAMe"  = "#d6a32b",
  "7 EnhAHet" = "#f0c36b",
  "8 EnhPois" = "#f5e1a1",
  "9 Repr"    = "#b0b0b0",
  "10 QuiW"   = "#999999"
)

data$state_label <- ifelse(
  data$state %in% names(state_map),
  state_map[data$state],
  data$state
)

state_levels <- names(state_colors)
data$state_label <- factor(data$state_label, levels = rev(state_levels))

############################################
############################################
summary_df <- data %>%
  group_by(state, state_label) %>%
  summarise(
    n_samples   = sum(!is.na(weighted_meth_pct)),
    mean_meth   = mean(weighted_meth_pct, na.rm = TRUE),
    median_meth = median(weighted_meth_pct, na.rm = TRUE),
    sd_meth     = sd(weighted_meth_pct, na.rm = TRUE),
    min_meth    = min(weighted_meth_pct, na.rm = TRUE),
    max_meth    = max(weighted_meth_pct, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  summary_df,
  file = OUTTSV,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################
############################################
p <- ggplot(data, aes(x = state_label, y = weighted_meth_pct, fill = state_label)) +
  geom_violin(
    trim = TRUE,
    alpha = 0.8,
    color = "white",
    linewidth = 0.2,
    scale = "width",
    width = 0.42,
    na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    alpha = 1,
    color = "black",
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = "Sample-level weighted methylation across chromatin states",
    x = NULL,
    y = "Weighted methylation (%)"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),

    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(color = "black", fill = "white", linewidth = 0.8)
  )

############################################
# output PDF
############################################
pdf(OUTPDF, width = 4, height = 8)
print(p)
dev.off()

message("PDF saved to: ", OUTPDF)
message("Summary saved to: ", OUTTSV)
