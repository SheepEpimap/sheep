setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(circlize)



#
library(circlize)


library(ggplot2)
library(stringr)
dir()
library(CMplot)
library(dplyr)


df = read.csv("251125_AS_Dis", sep = "\t", header = F)
head(df)

# Processing note.
df_count = df %>% group_by(V1) %>%
  mutate(count=n())  %>%
  ungroup() %>%
  select(V1, count) %>%
  unique()

df_count[order(df_count$count,decreasing = T),]


# custom_chr
custom_chr = df[,c("V2","V3")] %>%
  group_by(V2) %>%  # Processing note.
  slice_max(V3, n = 1) %>%  # Processing note.
  ungroup()  
custom_chr$start = rep(0, length(custom_chr$V2))
custom_chr = custom_chr[,c(1,3,2)]
head(custom_chr)
# head(custom_chr)
# custom_chr = rbind(c("Z",0,min(custom_chr$V3)),custom_chr)
# custom_chr$start = as.numeric(custom_chr$start)
# custom_chr$V3 = as.numeric(custom_chr$V3)

# 
head(df)
# ATAC_ASSnp_density <- df[df$V1=="ATAC", c("V2", "V3")]
ls_snp_density = list()
i = 0

for (reg in as.data.frame(df_count)[order(df_count$count,decreasing = T),"V1"]) {
  i = i + 1
  print(reg)
  ls_snp_density[[i]] = df[df$V1==reg, c("V2", "V3")]
}

#
# unique(df$V1)
circos.clear()
as.data.frame(df_count)[order(df_count$count,decreasing = T),"V1"]
custom_color = c("#3498db", "#4d34db", "#cb34db", "#db346e", "#db7734", "#c2db34", "#45db34", "#34dba1")
custom_color = c("#f768a1", "#fed976", "#41b6c4", "#41ab5d", "#DC0000FF", "#969696")
bin_size = 10e6

dev.off()
pdf("1.pdf")
# circos.initializeWithIdeogram(custom_chr)
circos.initializeWithIdeogram(custom_chr, 
                              axis.labels.cex = 1e-100,
                              labels.cex = 1.5
                              )

circos.trackHist(ls_snp_density[[1]]$V2, ls_snp_density[[1]]$V3,
                 bin.size = bin_size,
                 col = custom_color[[1]],
                 border = NA,
                 track.height = 0.12,
                 bg.border = NA)
circos.trackHist(ls_snp_density[[2]]$V2, ls_snp_density[[2]]$V3,
                 bin.size = bin_size,
                 border = NA,
                 col = custom_color[[2]],
                 track.height = 0.12,
                 bg.border = NA)


circos.trackHist(ls_snp_density[[3]]$V2, ls_snp_density[[3]]$V3,
                 bin.size = bin_size,
                 border = NA,
                 col = custom_color[[3]],
                 track.height = 0.12,
                 bg.border = NA)

circos.trackHist(ls_snp_density[[4]]$V2, ls_snp_density[[4]]$V3,
                 bin.size = bin_size,
                 border = NA,
                 col = custom_color[[4]],
                 track.height = 0.12,
                 bg.border = NA)

circos.trackHist(ls_snp_density[[5]]$V2, ls_snp_density[[5]]$V3,
                 bin.size = bin_size,
                 border = NA,
                 col = custom_color[[5]],
                 track.height = 0.12,
                 bg.border = NA)


# dev.off()
# circos.initializeWithIdeogram(custom_chr)
# custom_chr[custom_chr$V2=="chr25",]
# ls_snp_density[[6]]$V2
# ls_snp_density[[6]]$V3
circos.trackHist(
  # ls_snp_density[[6]]$V2, ls_snp_density[[6]]$V3,
  c("chr14","chr25"), c(63861494, 44795092),
                 bin.size = bin_size,
                 border = NA,
                 col = custom_color[[6]],
                 track.height = 0.12,
                 bg.border = NA)

# png("circos_plot.png", width = 800, height = 800, res = 150)
circos.clear()
dev.off()

# legend 
tmp4legend = as.data.frame(df_count)[order(df_count$count,decreasing = T), ]
tmp4legend$V1
tmp4legend$V1 = factor(tmp4legend$V1, levels=tmp4legend$V1)

ggplot(tmp4legend) +
  geom_col(aes(x=V1, y=count, fill=V1)) +
  scale_fill_manual(values = custom_color) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    axis.text.x = element_text(size=14, color="black", angle = 90, vjust = 0.5)
  ) +
  theme(
    plot.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "cm")  # Processing note.
  ) +
  guides(fill = guide_legend(nrow = 2))  # Processing note.


ggsave("1-legend.pdf", width = 22, height = 10, units = "cm")





# ##### demo #####
# # https://jokergoo.github.io/circlize_examples/
# 
# x = rnorm(2600)
# factors = sample(letters, 2600, replace = TRUE)
# circos.initialize(factors = factors, x = x)
# circos.trackHist(factors = factors, x = x, track.height = 0.1, col = "#999999", border = "#999999")
# circos.trackHist(factors = factors, x = x, force.ylim = FALSE, bin.size = 0.1, track.height = 0.1, col = "#999999", border = "#999999")
# # circos.trackHist(factors = factors, x = x, draw.density = TRUE, track.height = 0.1, col = "#999999", border = "#999999")
# # circos.trackHist(factors = factors, x = x, draw.density = TRUE, force.ylim = FALSE, track.height = 0.1, col = "#999999", border = "#999999")
# 
# circos.clear()
# 
# # 
# # circos.initializeWithIdeogram()
# Processing note.
# custom_chr <- data.frame(
#   chr = c("chr1", "chr2", "chr3", "chr4", "chr5"),
#   start = rep(0, 5),
# Processing note.
# )
# 
# circos.initializeWithIdeogram(custom_chr)
# 
# Processing note.
# Processing note.
# snp_density <- data.frame(
#   chr = sample(paste0("chr", 1:5), 2000, replace = TRUE),
#   pos = runif(2000, 1, 200e6)
# )
# 
# Processing note.
# gene_density <- data.frame(
#   chr = sample(paste0("chr", 1:5), 500, replace = TRUE),
#   pos = runif(500, 1, 200e6)
# )
# 
# circos.trackHist(snp_density$chr, snp_density$pos,
#                  bin.size = 10e6,
#                  col = "red",
#                  track.height = 0.12,
#                  bg.border = NA)
# 
# Processing note.
# circos.trackHist(gene_density$chr, gene_density$pos,
#                  bin.size = 10e6,
#                  col = "blue",
#                  track.height = 0.12,
#                  bg.border = NA)
# 
# circos.clear()