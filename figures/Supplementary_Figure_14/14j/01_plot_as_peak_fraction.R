setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(stringr)
library(ggbeeswarm)



df = read.csv("251222_ASPeakNum", sep = "\t", header = F)
head(df)

df$omics = str_split_fixed(df$V1,"_",2)[,1]
df$tis = str_split_fixed(df$V1,"_",2)[,2]
head(df)

unique(df$tis)

color_map = c(
  "ATAC" = "#f768a1",  # Processing note.
  "H3K27ac" = "#DC0000FF",  # Processing note.
  "H3K4me3" = "#41ab5d",  # Processing note.
  "H3K27me3" = "#969696",  # Processing note.
  "H3K4me1" ="#fed976",  # Processing note.
  "RNASeq"="#41b6c4"
)
ggplot(df) +
  geom_violin(aes(x=omics, y=(V2/V3)*100), color=NA, fill="lightgrey", alpha=.5) +
  geom_beeswarm(aes(x=omics, y=(V2/V3)*100, color=omics), alpha=0.75, shape=16) +
  theme_classic() +
  labs(x="", y="Proportion of the AS peaks (%)") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = "right",
    axis.text = element_text(size=14, color="black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size=16, color="black"),
    plot.title = element_text(size=16, color="black"),
  ) +
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = color_map  # Processing note.
    # Processing note.
  )

ggsave("fraction_of_peaks_containing_AS_SNPs.png")  
ggsave("fraction_of_peaks_containing_AS_SNPs.pdf")  
