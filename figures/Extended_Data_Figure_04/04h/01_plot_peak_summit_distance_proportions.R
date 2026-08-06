#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
hit_summit_dist_coarse <- hit_summit_dist %>%
  dplyr::select(tissue, motif_name, category, bin_dist, n_per_tissue_per_bin = n) %>%
  filter(category %in% c("resolved", "unresolved") &
           motif_name %in% motifs_keep) %>%
  mutate(bin_dist_coarse = case_when(
    bin_dist <= 50 ~ "0-50 bp",
    bin_dist <= 100 ~ "51-100 bp",
    bin_dist <= 150 ~ "101-150 bp",
    bin_dist <= 200 ~ "151-200 bp",
    bin_dist <= 250 ~ "201-250 bp",
    TRUE ~ "> 251 bp"
  )) %>%
  group_by(motif_name, , bin_dist_coarse) %>%
  summarize(n_per_bin = sum(n_per_tissue_per_bin)) %>%
  mutate(n_total = sum(n_per_bin),
         prop_ber_bin = n_per_bin / n_total)

cmap_bins <- RColorBrewer::brewer.pal(6, "BuGn")
names(cmap_bins) <- rev(c("0-50 bp", "51-100 bp", "101-150 bp", "151-200 bp", "201-250 bp", "> 251 bp"))

motif_order_dist <- hit_summit_dist_coarse %>% filter(bin_dist_coarse == "0-50 bp") %>% arrange(desc(prop_ber_bin)) %>% pull(motif_name) %>% unique()

p=hit_summit_dist_coarse %>%
  mutate(bin_dist_coarse = factor(bin_dist_coarse, levels = names(cmap_bins)),
         motif_name = factor(motif_name, levels = motif_order_dist)) %>%
  ggplot(aes(x = motif_name, y = prop_ber_bin)) +
  geom_bar(stat = "identity", aes(fill = bin_dist_coarse)) +
  scale_fill_manual(values = cmap_bins) +
  rotate_x()


ggsave(file.path(figout, "proportions_of_motifs_per_bin.pdf"), p,
       width = 13, height = 5, units = "in", device = cairo_pdf)
ggsave(file.path(figout, "proportions_of_motifs_per_bin.png"), p,
       width = 13, height = 5, units = "in", dpi = 300)
