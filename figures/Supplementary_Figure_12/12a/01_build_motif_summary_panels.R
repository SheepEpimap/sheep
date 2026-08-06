#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
resolved_broad <- motifs_compiled_unique %>%
  filter(category %in% c("resolved"))

length(resolved_broad$motif_name %>% unique())


unresolved_broad <- motifs_compiled_unique %>%
  filter(category %in% c("unresolved"))

length(unresolved_broad$motif_name %>% unique())
motif_eqtl <- read_tsv(
  "/data/home/sczd644/run/zsw_chrombpnet/03-syntax/01/motif_eqtl_split/motif_eqtl_count.tsv",
  skip = 1,                                    #   1  ( )
  col_names = c("merged_pattern", "n_variants"), #
  col_types = cols(merged_pattern = col_character(),
                   n_variants = col_integer())
)

motifs_to_plot <- unresolved_broad %>%
  left_join(motif_eqtl, by = "merged_pattern")


motif_order <- hits_all_anno %>%
  filter(motif_name %in% motifs_to_plot$motif_name) %>%
  group_by(motif_name) %>%
  summarize(mean_dist = mean(median_distToTSS)) %>%
  arrange(mean_dist) %>%  #
  pull(motif_name)

motifs_to_plot <- motifs_to_plot %>%
  mutate(motif_name = factor(motif_name, levels = motif_order)) %>%
  arrange(motif_name)


cwm_list_subset <- cwm_list[motifs_to_plot$merged_pattern_clean] %>%
  map(trim_cwm, flank = 2)


names(cwm_list_subset) <- plyr::mapvalues(
  names(cwm_list_subset),
  from = motifs_to_plot$merged_pattern_clean,
  to = as.character(motifs_to_plot$motif_name)
)

to_revcomp <- c("unresolved#1")
cwm_list_subset[to_revcomp] <- map(cwm_list_subset[to_revcomp], revcomp)


p1 <- stack_logo(cwm_list_subset, method = "custom", title = "De novo CWM")



known_motifs_subset <- jaspar_motifs[motifs_to_plot$match0]
names(known_motifs_subset) <- motifs_to_plot$TF

if (length(unique(names(known_motifs_subset))) < length(names(known_motifs_subset))) {
  names(known_motifs_subset) <- paste0(seq_along(known_motifs_subset), " ", names(known_motifs_subset))
}

known_motifs_subset <- map(known_motifs_subset, function(motif) {
  if (is.null(rownames(motif))) {
    rownames(motif) <- c("A", "C", "G", "T")
  }
  return(motif)
})

p2 <- stack_logo(known_motifs_subset, method = "bits", title = "Known PPM")

#library(viridis)
p3 <- motifs_to_plot %>%
  distinct(motif_name, total_hits) %>%
  mutate(motif_name = factor(motif_name, levels = rev(motif_order))) %>%
  ggplot(aes(y = motif_name, x = "1")) +
  geom_tile(aes(fill = log10(total_hits)), alpha = 1) +
  geom_text(aes(label = scales::comma(total_hits)), color = "white", size = 3) +
  scale_fill_gradientn(colors = viridis::plasma(100)[0:95], limits = c(0, 8)) +
  ylab(NULL) + xlab(NULL) + ggtitle("total # \nhits") +
  rotate_x() + hide_ticks() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), legend.position = "bottom")

p4 <- motifs_to_plot %>%
  distinct(motif_name, n_variants) %>%
  mutate(n_variants = as.numeric(n_variants),
         motif_name = factor(motif_name, levels = rev(motif_order))) %>%
  ggplot(aes(y = motif_name, x = n_variants)) +
  geom_col(fill = "gray70") +
  geom_text(aes(label = n_variants), hjust = -0.5) +
  ylab(NULL) + xlab(NULL) + ggtitle("# motifs") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlim(c(0, 20000))


p6 <- hits_all_anno %>%
  filter(motif_name %in% motif_order) %>%
  dplyr::select(motif_name, tissue, n) %>%
  group_by(motif_name, tissue) %>%
  summarise(n_hits = sum(n), .groups = "drop") %>%
  mutate(motif_name = factor(motif_name, levels = rev(motif_order))) %>%
  ggplot(aes(y = motif_name, x = n_hits)) +
  geom_col(aes(fill = tissue), position = "fill") +
  scale_fill_manual(values = cmap_organ) +
  ylab(NULL) + xlab(NULL) + ggtitle("% hits \nby tissue") +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    legend.position = "none"     # ←
  )


p7 <- hits_all_anno %>%
  filter(motif_name %in% motif_order) %>%
  dplyr::select(motif_name, tissue, Distal, Enhancer, Promoter) %>%
  pivot_longer(cols = c(Distal, Enhancer, Promoter), names_to = "region_type", values_to = "n_hits") %>%
  group_by(motif_name, region_type) %>%
  summarize(n_hits_total = sum(n_hits, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    region_type = factor(region_type, levels = names(cmap_region_type)),
    motif_name = factor(motif_name, levels = rev(motif_order))
  ) %>%   # ←   )
  ggplot(aes(y = motif_name, x = n_hits_total)) +
  geom_col(aes(fill = region_type), position = "fill", alpha = 0.7) +
  scale_fill_manual(values = cmap_region_type) +
  ylab(NULL) + xlab(NULL) + ggtitle("% hits \nby type") +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    legend.position = "none"     # ←
  )

p8 <- hits_all_anno %>%
  filter(motif_name %in% motif_order) %>%
  dplyr::select(motif_name, tissue, median_distToTSS) %>%
  mutate(motif_name = factor(motif_name, levels = rev(motif_order))) %>%
  ggplot(aes(y = motif_name, x = median_distToTSS)) +
  geom_boxplot(fill = "gray70", outliers = FALSE) +
  ylab(NULL) + xlab(NULL) + ggtitle("median distance \nto TSS") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_cartesian(xlim = c(0, 50000)) +
  scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale()))


p9 <- hits_all_anno %>%
  filter(motif_name %in% motif_order) %>%
  dplyr::select(motif_name, tissue, median_distToPeakSummit) %>%
  mutate(motif_name = factor(motif_name, levels = rev(motif_order))) %>%
  ggplot(aes(y = motif_name, x = median_distToPeakSummit)) +
  geom_boxplot(fill = "gray70", outliers = FALSE) +
  ylab(NULL) + xlab(NULL) + ggtitle("median distance \nto summit") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_cartesian(xlim = c(0, 300))

hit_dyad_dist_agg <- hit_dyad_dist %>%
  filter(motif_name %in% motif_order) %>%
  group_by(motif_name, bin_dist) %>%
  summarize(sum_counts = sum(n)) %>%
  mutate(zscore = scale(sum_counts)) %>%
  mutate(motif_name = factor(motif_name, levels = rev(motif_order)),
         bin_dist = as.numeric(bin_dist)) %>%
  ungroup()

p10 <- hit_dyad_dist_agg %>%
  ggplot(aes(x = bin_dist, y = motif_name)) +
  geom_tile(aes(fill = zscore)) +
  scale_fill_gradientn(colours = rdbu2,
                       rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  scale_x_continuous(breaks = seq(10, 250, by = 10)) +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    legend.position = "none"     # ←
  ) +
  xlab("distance bin (bp)") +
  ylab(NULL) +
  ggtitle("binned distances to \nnucleosome dyad")





final_plot <- cowplot::plot_grid(
  p1, p2,  p3, p7, p8, p9, p10, p6, p4,
  nrow = 1, align = "h", axis = "tb",
  rel_widths = c(2.5, 2.5, 1, 1, 2, 2, 1, 1, 1)
)

print(final_plot)


ggsave(file.path(figout, "unresolved_fig", "summary.pdf"), final_plot,
       width = 25, height = 20, units = "in", device = cairo_pdf)
ggsave(file.path(figout, "unresolved_fig", "summary.png"), final_plot,
       width = 25, height = 20, units = "in", dpi = 300)
ggsave(file.path(figout, "p8.pdf"), p8,
       width = 10, height = 15, units = "in", device = cairo_pdf)
ggsave(file.path(figout, "p8.png"), p8,
       width = 10, height = 15, units = "in", dpi = 300)
plot_list <- list(
  p1 = p1,
  p2 = p2,
  p3 = p3,
  p7 = p7,
  p8 = p8,
  p9 = p9,
  p10 = p10,
  p6 = p6,
  p4 = p4
)

out_one <- function(p, name) {
  pdf_path <- file.path(figout, "unresolved_fig", paste0(name, ".pdf"))
  png_path <- file.path(figout, "unresolved_fig", paste0(name, ".png"))
  # PDF
  ggsave(pdf_path, p, width = 10, height = 15, units = "in", device = cairo_pdf)
  # PNG
  ggsave(png_path, p, width = 10, height = 15, units = "in", dpi = 300)
  invisible(NULL)
}

walk(names(plot_list), ~ out_one(plot_list[[.x]], .x))

save.image(file = "/data/home/sczd644/run/zsw_chrombpnet/03-syntax/01/output/MySession.RData")
load("/data/home/sczd644/run/zsw_chrombpnet/03-syntax/01/output/MySession.RData")
