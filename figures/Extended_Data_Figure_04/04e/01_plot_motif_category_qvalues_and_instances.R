#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
p1 <- motifs_compiled_unique %>%
  mutate(similarity = qval0) %>%
  ggplot(aes(x = category, y = -log10(qval0))) +
  geom_boxplot(aes(fill = category), outliers = FALSE) +
  geom_jitter(shape = 21, width = 0.2, size = 0.5) +
  scale_fill_manual(values = cmap_category) +
  geom_hline(yintercept = -log10(0.05)) +
  no_legend() +
  xlab(NULL) +
  rotate_x()

p2 <- motifs_compiled_unique %>%
  ggplot(aes(x = category, y = log10(total_hits))) +
  geom_boxplot(aes(fill = category), outlier.colour = NA) +
  geom_jitter(color = "black", width = 0.2, size = 0.5, shape = 21) +
  scale_fill_manual(values = cmap_category) +
  rotate_x() +
  no_legend() +
  # ggtitle("Total # of hits per unique motif \neach point is a motif (n=118)") +
  coord_cartesian(ylim = c(1, 7))

p=plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb")
 ggsave(file.path(figout, "qvalues_hits.pdf"), p, w = 5, h = 5, u = "in",device = cairo_pdf); ggsave(file.path(figout, "qvalues_hits.png"), p, w = 5, h = 5, u = "in", dpi = 300)
