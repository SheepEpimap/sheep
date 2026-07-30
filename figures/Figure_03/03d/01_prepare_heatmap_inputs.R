#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
library(dplyr)
library(tidyr)

# ---- inputfile ----
infile <- "../final_merged.tsv"

df <- read.delim(infile, stringsAsFactors = FALSE, check.names = FALSE)

df_clean <- df %>%
  dplyr::select(
    tissue = input_tissue,
    seqlets = n_seqlets,
    cluster = Transcription_Factor
  ) %>%
  dplyr::filter(
    !is.na(cluster),
    !is.na(seqlets),
    cluster != "Transcription_Factor"
  )

df_clean$seqlets <- as.numeric(df_clean$seqlets)

custom_sort <- function(cluster_names) {
  unresolved <- cluster_names[grepl("unresolved", cluster_names)]
  resolved <- cluster_names[!grepl("unresolved", cluster_names)]

  resolved_sorted <- sort(resolved)

  if(length(unresolved) > 0) {
    unresolved_nums <- as.numeric(gsub("unresolved#", "", unresolved))
    unresolved_sorted <- unresolved[order(unresolved_nums)]
  } else {
    unresolved_sorted <- character(0)
  }

  c(resolved_sorted, unresolved_sorted)
}

tissue_colors <- c(
   "cerebral-cortex" = "#dcd71a",
   "midbrain" = "#f9ed19",
   "cerebellum" = "#efd80b",
   "brainstem" = "#f4d578",
   "hippocampus" = "#f1cb05",
   "hypothalamus" = "#825e19",
   "medulla-oblongata" = "#d0b35b",
   "optic-chiasm" = "#fece01",
   "pineal" = "#807120",
   "pituitary" = "#f1d95d",
   "pons" = "#fcd222",
   "splenium" = "#7c6919",
   "rumen" = "#fc9891",
   "reticulum" = "#d97c68",
   "omasum" = "#f8ae81",
   "abomasum" = "#f18264",
   "duodenum" = "#eb9d63",
   "jejunum" = "#eb951c",
   "ileum" = "#ce9639",
   "cecum" = "#daaa6c",
   "colon" = "#f2c063",
   "rectum" = "#efe0a8",
   "cervix" = "#b6d7a9",
   "cornua-uteri" = "#69d683",
   "corpus-uteri" = "#80d897",
   "ovary" = "#69d28c",
   "oviduct" = "#79ffaa",
   "epididymis" = "#70d24b",
   "testis" = "#7ef351",
   "mammary-gland" = "#fed9d0",
   "bone-marrow" = "#d84b4b",
   "lymph-node" = "#c33a11",
   "thymus" = "#ff3a32",
   "thyroid" = "#f359d1",
   "spleen" = "#962932",
   "liver" = "#ad8c8b",
   "kidney" = "#4f3136",
   "lung" = "#36b5f1",
   "heart" = "#bc58e3",
   "muscle" = "#a180ca",
   "adipose" = "#ffc4f1",
   "skin" = "#d09dc5",
   "soft-horn" = "#a25d73"
)

df_sum <- df_clean %>%
  group_by(tissue, cluster) %>%
  summarise(seqlets = sum(seqlets, na.rm = TRUE), .groups = "drop")

all_clusters <- unique(df_sum$cluster)
cluster_order <- custom_sort(all_clusters)

tissue_order <- names(tissue_colors)
tissue_order <- tissue_order[tissue_order %in% unique(df_sum$tissue)]

df_sum$tissue <- factor(df_sum$tissue, levels = tissue_order)
df_sum$cluster <- factor(df_sum$cluster, levels = cluster_order)

heatmap_data <- df_sum %>%
  pivot_wider(
    names_from = tissue,
    values_from = seqlets,
    values_fill = 0
  ) %>%
  arrange(match(cluster, cluster_order))

write.table(heatmap_data, "heatmap_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message("✅  save: heatmap_data.tsv")

bottom_bar_data <- df_sum %>%
  group_by(tissue) %>%
  summarise(n_clusters = n_distinct(cluster)) %>%
  arrange(match(tissue, tissue_order))

write.table(bottom_bar_data, "bottom_bar_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message("✅  save: bottom_bar_data.tsv")

right_stacked_data <- df_sum %>%
  group_by(cluster) %>%
  mutate(
    total_seqlets = sum(seqlets),
    prop = seqlets / total_seqlets
  ) %>%
  ungroup()

write.table(right_stacked_data, "right_stacked_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message("✅  save: right_stacked_data.tsv")

tissue_color_df <- data.frame(
  tissue = names(tissue_colors),
  color = tissue_colors
)
write.table(tissue_color_df, "tissue_colors.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message("✅  save: tissue_colors.tsv")

cluster_order_df <- data.frame(
  cluster = cluster_order,
  order = 1:length(cluster_order)
)
write.table(cluster_order_df, "cluster_order.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message("✅  save: cluster_order.tsv")

message("📊  file ！")
message("📁  file :")
message("  - heatmap_data.tsv")
message("  - bottom_bar_data.tsv")
message("  - right_stacked_data.tsv")
message("  - tissue_colors.tsv")
message("  - cluster_order.tsv")
