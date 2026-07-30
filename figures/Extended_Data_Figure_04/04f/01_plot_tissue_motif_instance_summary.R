#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
tissue_color <- readr::read_tsv("/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv",
                                col_types = cols(tissue = col_character(), color = col_character()))
p0 <- peaks_w_hits %>%
      mutate(pattern_class = "pos_patterns",
             tissue        = factor(tissue, levels = tissue_color$tissue)) %>%
      ggplot(aes(x = tissue, y = Prop_peaks_w_hits)) +
      geom_jitter(stat = "identity", aes(fill = tissue), shape = 21,
                  color = "white", size = 2, alpha = 0.6, width = 0.25) +
      scale_fill_manual(values = setNames(tissue_color$color, tissue_color$tissue)) +
      rotate_x() +
      facet_wrap(~ pattern_class, ncol = 1) +
      labs(y = "Proportion of peaks with hits", title = "Peaks with hits") +
      theme(legend.position = "none",
            aspect.ratio = 0.2,
            strip.background = element_blank(),
            panel.grid.major.y = element_line(color = "gray90"),
            panel.grid.minor.y = element_line(color = "gray90"))

tissue_color <- readr::read_tsv("/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv",
                                col_types = cols(tissue = col_character(), color = col_character()))
p1 <- hits_all %>%
      mutate(pattern_class = "pos_patterns") %>%
      group_by(tissue, pattern_class) %>%
      summarise(total_hits = sum(n), .groups = "drop") %>%
      mutate(tissue = factor(tissue, levels = tissue_color$tissue)) %>%
      ggplot(aes(x = tissue, y = total_hits)) +
      geom_jitter(stat = "identity", aes(fill = tissue), shape = 21,
                  color = "white", size = 2, alpha = 0.6, width = 0.25) +
      scale_fill_manual(values = setNames(tissue_color$color, tissue_color$tissue)) +
      rotate_x() +
      facet_wrap(~ pattern_class, ncol = 1) +
      labs(y = "Total number of hits", title = "Total hits") +
      theme(legend.position = "none",
            aspect.ratio = 0.2,
            strip.background = element_blank(),
            panel.grid.major.y = element_line(color = "gray90"),
            panel.grid.minor.y = element_line(color = "gray90"))

tissue_color <- readr::read_tsv("/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv",
                                col_types = cols(tissue = col_character(), color = col_character()))

p2 <- hits_per_peak %>%
      count(tissue, peak_id) %>%
      group_by(tissue) %>%
      summarise(median_n_hits = median(n, na.rm = TRUE), .groups = "drop") %>%
      mutate(pattern_class = "median_pos_patterns",
             tissue = factor(tissue, levels = tissue_color$tissue)) %>%
      ggplot(aes(x = tissue, y = median_n_hits)) +
      ggbeeswarm::geom_quasirandom(aes(fill = tissue), shape = 21,
                                   color = "white", size = 2, alpha = 0.6) +
      scale_fill_manual(values = setNames(tissue_color$color, tissue_color$tissue)) +
      facet_wrap(~ pattern_class, ncol = 1) +
      rotate_x() +
      labs(y = "Median number of hits per peak", title = "Median hits per peak") +
      theme(legend.position = "none",
            aspect.ratio = 0.2,
            strip.background = element_blank(),
            panel.grid.major.y = element_line(color = "gray90"),
            panel.grid.minor.y = element_line(color = "gray90")) +
      scale_y_continuous(breaks = seq(0, 15, 1), limits = c(0, 15))

p0 <- p0 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

p1 <- p1 +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

p2 <- p2 +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
library(patchwork)
p_final <- (p0 / p1 / p2) +  #
           plot_layout(heights = c(1, 1.2, 1))  #   rel_heights

ggsave(file.path(figout, "hits_per_peak.pdf"), p_final,
       width = 6, height = 10, units = "in", device = cairo_pdf)
ggsave(file.path(figout, "hits_per_peak.png"), p_final,
       width = 6, height = 10, units = "in", dpi = 300)
