#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
tissue_color <- readr::read_tsv(
  "/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv",
  col_types = cols(tissue = col_character(), color = col_character())
)
cmap_organ <- setNames(tissue_color$color, tissue_color$tissue)

total_hits_per_tissue <- hits_all %>%
  group_by(tissue) %>%
  summarise(total_hits = sum(n), .groups = "drop") %>%
  left_join(
    hits_per_peak %>% distinct(tissue, peak_id) %>% count(tissue, name = "n_peaks"),
    by = "tissue"
  )

cor_peaks_total <- lm(total_hits ~ n_peaks, data = total_hits_per_tissue)

p2 <- total_hits_per_tissue %>%
  ggplot(aes(x = n_peaks, y = total_hits)) +
  geom_point(aes(fill = tissue), shape = 21, size = 3, alpha = 0.7, color = "white") +
  scale_fill_manual(values = cmap_organ) +
  scale_y_continuous(labels = scales::label_number(scale = 1/1e6)) +
  scale_x_continuous(labels = scales::comma) +
  annotate("text",
           label = glue("cor = {round(sqrt(summary(cor_peaks_total)$r.squared), 3)}, \nR² = {round(summary(cor_peaks_total)$r.squared, 3)}"),
           x = max(total_hits_per_tissue$n_peaks) * 0.9,
           y = max(total_hits_per_tissue$total_hits) * 0.9,
           size = 5) +
  xlab("Total peaks") +
  ylab("Total hits (×1M)") +
  cowplot::theme_cowplot() +
  ggtitle("Total hits vs.\n# peaks") +
  theme(legend.position = "none")

hits_per_peak_count <- hits_per_peak %>%
  count(tissue, peak_id, name = "n_hits_per_peak")

median_hits_per_tissue <- hits_per_peak_count %>%
  group_by(tissue) %>%
  summarise(median_pos_patterns = median(n_hits_per_peak, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    hits_per_peak_count %>% count(tissue, name = "n_peaks"),
    by = "tissue"
  )

cor_peaks_median <- lm(median_pos_patterns ~ n_peaks, data = median_hits_per_tissue)

p4 <- median_hits_per_tissue %>%
  ggplot(aes(x = n_peaks, y = median_pos_patterns)) +
  geom_point(aes(fill = tissue), shape = 21, size = 3, alpha = 0.7, color = "white") +
  scale_fill_manual(values = cmap_organ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  annotate("text",
           label = glue("cor = {round(sqrt(summary(cor_peaks_median)$r.squared), 3)}, \nR² = {round(summary(cor_peaks_median)$r.squared, 3)}"),
           x = max(median_hits_per_tissue$n_peaks) * 0.9,
           y = max(median_hits_per_tissue$median_pos_patterns) * 0.9,
           size = 5) +
  xlab("Total peaks") +
  ylab("Median # positive hits per peak") +
  cowplot::theme_cowplot() +
  ggtitle("Median hits per peak vs.\n# peaks") +
  theme(legend.position = "none")

median_hits_per_tissue <- median_hits_per_tissue %>%
  mutate(n_peak_quartile = ntile(n_peaks, n = 4))

p6 <- median_hits_per_tissue %>%
  ggplot(aes(x = factor(n_peak_quartile), y = median_pos_patterns, group = n_peak_quartile)) +
  geom_boxplot(fill = "gray90", outliers = FALSE) +
  ggbeeswarm::geom_quasirandom(aes(fill = tissue), shape = 21, color = "white",
                               size = 4, alpha = 0.6, width = 0.2) +
  scale_fill_manual(values = cmap_organ) +
  scale_x_discrete(labels = c("Q1", "Q2", "Q3", "Q4")) +
  xlab("Quantile (total peaks)") +
  ylab("Median # positive hits per peak") +
  cowplot::theme_cowplot() +
  ggtitle("Median hits per peak vs.\npeaks quartile") +
  theme(legend.position = "none")

p_peaks_only <- plot_grid(p2, p4, p6,
                          nrow = 1,
                          align = "h",
                          axis = "tb",
                          rel_widths = c(1, 1.2, 1))

ggsave(file.path(figout, "hits_vs_peaks.pdf"), p_peaks_only,
       width = 15, height = 4, units = "in", device = cairo_pdf)
ggsave(file.path(figout, "hits_vs_peaks.png"), p_peaks_only,
       width = 15, height = 4, units = "in", dpi = 300)
