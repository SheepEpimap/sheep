
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(stringr)
library(dplyr)



##### seq-specific #####

library(UpSetR)

df = read.csv("251204_4upsetR", sep = "\t", header = F)


head(df)
unique(df$V3)
# paste0(df[df$V1=="ATAC","V2"], ":",df[df$V1=="ATAC","V3"])
dev.off()



listInput <- list(
  ATAC = paste0(df[df$V3=="ATAC","V2"], ":",df[df$V3=="ATAC","V1"]), 
  H3K27ac = paste0(df[df$V3=="H3K27ac","V2"], ":",df[df$V3=="H3K27ac","V1"]), 
  H3K27me3 = paste0(df[df$V3=="H3K27me3","V2"], ":",df[df$V3=="H3K27me3","V1"]),
  H3K4me1 = paste0(df[df$V3=="H3K4me1","V2"], ":",df[df$V3=="H3K4me1","V1"]),
  H3K4me3 = paste0(df[df$V3=="H3K4me3","V2"], ":",df[df$V3=="H3K4me3","V1"]),
  RNASeq = paste0(df[df$V3=="RNASeq","V2"], ":",df[df$V3=="RNASeq","V1"])
)

# Processing note.
# png("251204_Upset_Assay.png", width = 3600, height = 2400, res = 300)
pdf("251204_Upset_Assay.pdf", width = 36, height = 24)
#
upset(fromList(listInput), nsets = 6, order.by = "freq", point.size = 3.5,
      text.scale = c(1.3, 1.3, 1, 1, 2, 2),
      expression = "intersections_size > 100"
)
dev.off()


# ggsave("251204_Upset_Assay.png", width=20, height = 17.5, units = "cm")
# ?upset


##### tissue-specific #####
df = read.csv("251125_AS_TisNum", sep = "\t", header = F)

head(df)
unique(df$V1)
df$TisNum = str_split_fixed(df$V2,":",2)[,2]
df$TisNum = as.numeric(df$TisNum)
head(df)

# Processing note.
df_group = df %>%
  group_by(V1) %>%
  mutate(
    group_total = sum(V3),  # Processing note.
    percentage = V3 / group_total * 100,  # Processing note.
    proportion = V3 / group_total  # Processing note.
  ) %>%
  ungroup() %>%
  select(V1, V2, V3, TisNum, group_total, percentage, proportion) %>%
  arrange(V1, V2, V3)

color_map = c(
  "ATAC" = "#f768a1",  # Processing note.
  "H3K27ac" = "#DC0000FF",  # Processing note.
  "H3K4me3" = "#41ab5d",  # Processing note.
  "H3K27me3" = "#969696",  # Processing note.
  "H3K4me1" ="#fed976",  # Processing note.
  "RNASeq"="#41b6c4"
)


ggplot(df_group) +
  geom_line(aes(x=TisNum, y=percentage, color=V1), linewidth=1.25, alpha=0.7) +
  geom_point(aes(x=TisNum, y=percentage, color=V1), size=2, alpha=0.7) +
  theme_classic() +
  labs(y="Proportion (%)") +
  # facet_wrap(~V1, ncol = 1) +
  theme(
    axis.text = element_text(size = 14, colour = "black", vjust = 0.5,hjust = 1),
    axis.title = element_text(size = 14, colour = "black"),
    legend.position = "none",
    legend.text = element_text(size = 14, colour = "black"),
    legend.title = element_blank()
  ) +
  annotate('text', x=8, y=89.3, label="ATAC (89.3%)    ", size=5) +
  annotate('text', x=8, y=98.2-3, label="H3K27ac (98.2%) ", size=5) +
  annotate('text', x=8, y=76.7, label="H3K27me3 (76.7%)", size=5) +
  annotate('text', x=8, y=99.4, label="H3K4me1 (99.4%) ", size=5) +
  
  annotate('text', x=8, y=82.6, label="H3K4me3 (82.6%) ", size=5) +
  annotate('text', x=8, y=69.0, label="RNASeq (69.0%)  ", size=5) +
  geom_segment(mapping = aes(x=1, xend=4,y=89.3), linetype=3) +
  geom_segment(mapping = aes(x=1, xend=4,y=98.2, yend=98.2-3), linetype=3) +
  geom_segment(mapping = aes(x=1, xend=4,y=76.7), linetype=3) +
  geom_segment(mapping = aes(x=1, xend=4,y=99.4), linetype=3) +
  
  geom_segment(mapping = aes(x=1, xend=4,y=82.6), linetype=3) +
  geom_segment(mapping = aes(x=1, xend=4,y=69.0), linetype=3) +
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = color_map  # Processing note.
    # Processing note.
  )

ggsave("260408_Tis-specific.png", width=17.5, height = 12.5, units = "cm")
ggsave("260408_Tis-specific.pdf")

