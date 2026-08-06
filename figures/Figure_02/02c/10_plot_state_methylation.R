#!/usr/bin/env Rscript

# Figure 2c: plot sample-level weighted DNA methylation across E1-E11.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

arguments <- commandArgs(trailingOnly = TRUE)

default_input <- paste0(
  "/vol2/zhangshiwen/rrbs/wgbs_bismark/",
  "state_methylation_by_sample_tissue/",
  "all_samples_sample_state_weighted_methylation.tsv"
)
default_output_dir <- paste0(
  "/vol2/zhangshiwen/rrbs/wgbs_bismark/",
  "state_methylation_by_sample_tissue/plot"
)

input_file <- if (length(arguments) >= 1) arguments[[1]] else default_input
output_dir <- if (length(arguments) >= 2) arguments[[2]] else default_output_dir

if (!file.exists(input_file)) {
  stop("Input file was not found: ", input_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- file.path(
  output_dir,
  "Figure_02c_chromatin_state_methylation_E1_E11.pdf"
)
output_summary <- file.path(
  output_dir,
  "Figure_02c_chromatin_state_methylation_E1_E11.summary.tsv"
)

required_columns <- c(
  "sample",
  "sample_tissue",
  "state_tissue",
  "state",
  "methylated_count",
  "unmethylated_count",
  "total_count",
  "cpg_sites",
  "weighted_meth_pct"
)

methylation <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

missing_columns <- setdiff(required_columns, colnames(methylation))
if (length(missing_columns) > 0) {
  stop(
    "The input table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

states <- paste0("E", seq_len(11))

methylation <- methylation %>%
  mutate(
    state = as.character(state),
    weighted_meth_pct = suppressWarnings(
      as.numeric(weighted_meth_pct)
    )
  ) %>%
  filter(
    state %in% states,
    !is.na(weighted_meth_pct),
    weighted_meth_pct >= 0,
    weighted_meth_pct <= 100
  )

if (nrow(methylation) == 0) {
  stop("No valid E1-E11 methylation values were available for plotting.")
}

duplicated_sample_states <- methylation %>%
  count(sample, state, name = "row_count") %>%
  filter(row_count != 1)

if (nrow(duplicated_sample_states) > 0) {
  stop("At least one sample-state combination has more than one value.")
}

missing_states <- setdiff(states, unique(methylation$state))
if (length(missing_states) > 0) {
  stop(
    "No valid methylation values were available for: ",
    paste(missing_states, collapse = ", ")
  )
}

state_labels <- c(
  "E1" = "E1 TssA",
  "E2" = "E2 TssFlnk",
  "E3" = "E3 TssWk",
  "E4" = "E4 TssWBiv",
  "E5" = "E5 EnhA",
  "E6" = "E6 EnhAMe",
  "E7" = "E7 EnhAHet",
  "E8" = "E8 EnhPois",
  "E9" = "E9 Repr",
  "E10" = "E10 QuiW",
  "E11" = "E11 Qui"
)

state_colors <- c(
  "E1 TssA" = "#FF0000",
  "E2 TssFlnk" = "#FF2C79",
  "E3 TssWk" = "#F08080",
  "E4 TssWBiv" = "#FCE0E0",
  "E5 EnhA" = "#FFFF00",
  "E6 EnhAMe" = "#D6A32B",
  "E7 EnhAHet" = "#F0C36B",
  "E8 EnhPois" = "#F5E1A1",
  "E9 Repr" = "#B0B0B0",
  "E10 QuiW" = "#F2F2F2",
  "E11 Qui" = "#FFFFFF"
)

methylation <- methylation %>%
  mutate(
    state_label = unname(state_labels[state]),
    state_label = factor(
      state_label,
      levels = rev(unname(state_labels))
    )
  )

summary_table <- methylation %>%
  mutate(
    state = factor(state, levels = states)
  ) %>%
  group_by(state) %>%
  summarise(
    n_samples = n(),
    total_cpg_sites = sum(cpg_sites, na.rm = TRUE),
    mean_methylation_pct = mean(weighted_meth_pct),
    median_methylation_pct = median(weighted_meth_pct),
    sd_methylation_pct = sd(weighted_meth_pct),
    minimum_methylation_pct = min(weighted_meth_pct),
    maximum_methylation_pct = max(weighted_meth_pct),
    .groups = "drop"
  ) %>%
  arrange(state)

write.table(
  summary_table,
  file = output_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

figure_02c <- ggplot(
  methylation,
  aes(
    x = state_label,
    y = weighted_meth_pct,
    fill = state_label
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    width = 0.72,
    color = "black",
    linewidth = 0.25,
    na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.10,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.25,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = state_colors,
    breaks = names(state_colors),
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Methylation level (%)"
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 6
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(
      color = "black",
      size = 7,
      margin = margin(t = 2)
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.35
    ),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2)
  )

pdf(
  output_pdf,
  width = 1.45,
  height = 3.10,
  useDingbats = FALSE
)
print(figure_02c)
invisible(dev.off())

message("Figure saved to: ", output_pdf)
message("Summary saved to: ", output_summary)

