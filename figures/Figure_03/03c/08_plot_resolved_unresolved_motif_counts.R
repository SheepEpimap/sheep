#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
cmap_category <- c(
  resolved   = "#B3D9FF",   #
  unresolved = "#FFF2B3"    #
)
library(ggplot2)

cmap_category <- c(resolved = "#B3D9FF", unresolved = "#FFF2B3")

p1 <- motifs_compiled_all %>%
  group_by(category) %>%
  count() %>%
  mutate(category = factor(category, levels = rev(names(cmap_category)))) %>%
  ggplot(aes(x = category, y = n)) +
  geom_col(aes(fill = category)) +
  geom_text(aes(label = n), hjust = -0.5, fontface = "bold", size = 4) +
  coord_flip(ylim = c(0, 280), expand = FALSE) +   #   ylim
  scale_fill_manual(values = cmap_category) +
  ggtitle("De novo deep learning derived motifs \n(after MoDISco agggregation step)") +
  theme(legend.position = "none",
        panel.grid.major.x = element_line(color = "gray90"),
        panel.grid.minor.x = element_line(color = "gray90")) +
  theme(aspect.ratio = .25)

p2 <- motifs_compiled_all %>%
  ggplot(aes(x = "breakdown")) +
  geom_bar(aes(fill = category), position = "fill") +
  scale_fill_manual(values = cmap_category) +
  theme(aspect.ratio = .9)

p=plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", rel_widths = c(0.65, 0.35))


 ggsave(file.path(figout, "pattern_motif_n.pdf"), p, w = 10, h = 5, u = "in"); ggsave(file.path(figout, "pattern_motif_n.png"), p, w = 10, h = 5, u = "in", dpi = 300)
