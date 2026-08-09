setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(stringr)
library(dplyr)


df = read.csv("251204_TSS_Dis", sep = "\t", header = F)
head(df)


# Processing note.
min_val <- min(df$V2, na.rm = TRUE)
max_val <- max(df$V2, na.rm = TRUE)
bin_width <- 10000  # Processing note.
breaks <- seq(floor(min_val/bin_width)*bin_width, 
              ceiling(max_val/bin_width)*bin_width, 
              by = bin_width)

# Processing note.
df_bin <- df %>%
  group_by(V1) %>%
  mutate(
    value_bin = cut(V2, 
                    breaks = breaks,
                    include.lowest = TRUE)
  ) %>%
  count(V1, value_bin, name = "count") %>%
  arrange(V1, value_bin)

# prop
head(df_bin)
df_bin = df_bin %>%
  group_by(V1) %>%
  mutate(
    group_total = sum(count),  # Processing note.
    percentage = count / group_total * 100,  # Processing note.
    # Processing note.
  ) %>%
  ungroup() 

head(df_bin)
df_bin$bin_end = as.numeric(sub("]","", str_split_fixed(df_bin$value_bin,",",2)[,2]))

ggplot(df_bin) +
  geom_line(aes(x=bin_end/1000, y=percentage, color=V1), linewidth=1) +
  theme_classic() +
  scale_x_continuous(limits = c(-100000/1000,100000/1000)) +
  labs(x="Distance to TSS (Kb)", y="Proportion of AS SNPs (%)") +
  theme(
    legend.title = element_blank(),
    legend.text =  element_text(size = 14, color="black"),
    legend.position = c(0.8,0.7),
    axis.title = element_text(size = 16, color="black"),
    axis.text = element_text(size = 14, color="black"),
  ) +
  scale_color_manual(
    limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = c(
      "ATAC" = "#f768a1",  # Processing note.
      "H3K27ac" = "#DC0000FF",  # Processing note.
      "H3K4me3" = "#41ab5d",  # Processing note.
      "H3K27me3" = "#969696",  # Processing note.
      "H3K4me1" ="#fed976",  # Processing note.
      "RNASeq"="#41b6c4")  # Processing note.
    # Processing note.
  )

ggsave("TSS.png", width = 11, height = 8.5, units = "cm")
ggsave("TSS.pdf", width = 11, height = 8.5, units = "cm")

unique(df_bin$V1)
