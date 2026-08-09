setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(ggplot2)
library(patchwork)
library(ggrepel)
library(ggbeeswarm)
library(stringr)
library(tidyr)
library(tibble)
library(scales)

color_map = c(
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

color_map_assay = c(
  ATAC = "#f768a1",  # Processing note.
  H3K27ac = "#DC0000FF",  # Processing note.
  H3K4me3 = "#41ab5d",  # Processing note.
  H3K27me3 = "#969696",  # Processing note.
  H3K4me1 ="#fed976",  # Processing note.
  RNASeq ="#41b6c4"
)
##### 260406_AS-SNP_ChrState_fold_enrichment #####

df = read.csv("260406_AS-SNP_ChrState_fold_enrichment", sep = "\t", header = F)

df$V3 = factor(df$V3, levels = paste0("E",c(1:10)))
df$V2 = str_split_fixed(df$V2,"_",2)[,2]
head(df)


# colored by assay
ggplot(df) +
  geom_quasirandom(aes(x=V3, y=(V4), color=V2), alpha=0.5) +
  theme_classic() +
  theme_classic()+
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = color_map_assay  # Processing note.
    # Processing note.
  ) +
  theme(
    axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
    axis.title.y = element_text(size=16, color="black"),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.text = element_text(size=14, color="black"),
    legend.title = element_blank()
  ) +
  labs(y='Fold enrichment', title = "AS-SNP (colored by assay)")
ggsave("260408_AS-SNP_ChrState_FoldEnrichment(byAssay).png")
ggsave("260408_AS-SNP_ChrState_FoldEnrichment(byAssay).pdf")

# colored by tissue
ggplot(df) +
  geom_quasirandom(aes(x=V3, y=(V4), color=V1), alpha=0.5) +
  theme_classic() +
  theme_classic()+
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = color_map  # Processing note.
    # Processing note.
  ) +
  theme(
    axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
    axis.title.y = element_text(size=16, color="black"),
    axis.title.x = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=14, color="black"),
    legend.title = element_blank()
  ) +
  labs(y='Fold enrichment', title = "AS-SNP (colored by tissue)")
ggsave("260408_AS-SNP_ChrState_FoldEnrichment(byTissue).png")
ggsave("260408_AS-SNP_ChrState_FoldEnrichment(byTissue).pdf")



##### QTL #####
dir(pattern = "NON")
df = read.csv("260406_QTL_state_fold_enrichment", sep = "\t", header = F)
head(df)

df_nonAS = read.csv("260406_QTL_NON-AS-state_fold_enrichment", sep = "\t", header = F)
head(df_nonAS)


df$V3 = factor(df$V3, levels = paste0("E",c(1:10)))
df_nonAS$V3 = factor(df_nonAS$V3, levels = paste0("E",c(1:10)))

unique(df$V2)



# ggplot(df[(df$V2=="eQTL") & (df$V3=="E5"), ]) +
#   geom_violin(aes(x=V3, y=V4), fill="lightgrey", color=NA, alpha=0.5) +
#   geom_quasirandom(aes(x=V3, y=V4, color=V1)) +
#   theme_classic() +
#   scale_color_manual(
#     # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
# Processing note.
# Processing note.
#   ) +
#   theme(
#     axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
#     axis.title.y = element_text(size=16, color="black"),
#     axis.title.x = element_blank(),
#     legend.position = "none"
#   ) +
#   labs(y='Fold enrichment', title = "eQTL")


# E5 only
tsign = t.test(df[(df$V2=="eQTL") & (df$V3=="E5"), "V4"], df_nonAS[(df$V2=="eQTL") & (df_nonAS$V3=="E5"), "V4"], paired = TRUE)
tsign$p.value
tmp_max = max(c(df[(df$V2=="eQTL") & (df$V3=="E5"), "V4"], df_nonAS[(df$V2=="eQTL") & (df_nonAS$V3=="E5"), "V4"]))
tmp_p = ifelse(tsign$p.value < 0.05, ifelse(tsign$p.value<0.01, ifelse(tsign$p.value<0.001, "***", "**"),"*"), "")

ggplot() +
  geom_violin(df[(df$V2=="eQTL") & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
  geom_quasirandom(df[(df$V2=="eQTL") & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4, color=V1)) +
  geom_violin(df_nonAS[(df$V2=="eQTL") & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
  geom_quasirandom(df_nonAS[(df$V2=="eQTL") & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4, color=V1)) +
  geom_text(mapping = aes(x=1.5, y=tmp_max, label = tmp_p), size=6) +
  geom_segment(mapping = aes(x=1, xend=2, y=tmp_max*0.98)) +
  theme_classic()+
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = color_map  # Processing note.
    # Processing note.
  ) +
  theme(
    axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
    axis.title.y = element_text(size=16, color="black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y='Fold enrichment', title = "eQTL in AS-E5")
ggsave("260408_AS-ChrState_eQTL_AS-E5.png", width = 5, height = 7.5, units = "cm")
ggsave("260408_AS-ChrState_eQTL_AS-E5.pdf", width = 5, height = 7.5, units = "cm")


# E5 enQTL
tsign = t.test(df[(df$V2=="enQTL") & (df$V3=="E5"), "V4"], df_nonAS[(df$V2=="enQTL") & (df_nonAS$V3=="E5"), "V4"], paired = T)
tsign$p.value
tmp_max = max(c(df[(df$V2=="enQTL") & (df$V3=="E5"), "V4"], df_nonAS[(df$V2=="enQTL") & (df_nonAS$V3=="E5"), "V4"]))
tmp_p = ifelse(tsign$p.value < 0.05, ifelse(tsign$p.value<0.01, ifelse(tsign$p.value<0.001, "***", "**"),"*"), "")

ggplot() +
  geom_violin(df[(df$V2=="enQTL") & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
  geom_quasirandom(df[(df$V2=="enQTL") & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4, color=V1)) +
  geom_violin(df_nonAS[(df$V2=="enQTL") & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
  geom_quasirandom(df_nonAS[(df$V2=="enQTL") & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4, color=V1)) +
  geom_text(mapping = aes(x=1.5, y=tmp_max, label = tmp_p), size=6) +
  geom_segment(mapping = aes(x=1, xend=2, y=tmp_max*0.98)) +
  theme_classic()+
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = color_map  # Processing note.
    # Processing note.
  ) +
  theme(
    axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
    axis.title.y = element_text(size=16, color="black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y='Fold enrichment', title = "enQTL in AS-E5")

ggsave("260408_AS-ChrState_enQTL_AS-E5.png", width = 5, height = 7.5, units = "cm")
ggsave("260408_AS-ChrState_enQTL_AS-E5.pdf", width = 5, height = 7.5, units = "cm")


# # for (qtl in unique(df$V2))
#     
#   {
#   
#   ggplot(df[df$V2==qtl, ]) +
#     geom_violin(aes(x=V3, y=V4), fill="lightgrey", color=NA, alpha=0.5) +
#     geom_quasirandom(aes(x=V3, y=V4, color=V1)) +
#     theme_classic() +
#     scale_color_manual(
#       # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
# Processing note.
# Processing note.
#     ) +
#     theme(
#       axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
#       axis.title.y = element_text(size=16, color="black"),
#       axis.title.x = element_blank(),
#       legend.position = "none"
#     ) +
#     labs(y='Fold enrichment', title = qtl)
#   ggsave(paste0("260407_",qtl,".png"), width = 12, height = 5.5, units = "cm")
# }

##### gwas #####

df = read.csv("260406_GWAS_state_fold_enrichment", sep = "\t", header = F)
head(df)
df$V3 = factor(df$V3, levels = paste0("E",c(1:10)))

df_nonAS = read.csv("260406_GWAS_NON-AS-state_fold_enrichment", sep = "\t", header = F)
df_nonAS$V3 = factor(df_nonAS$V3, levels = paste0("E",c(1:10)))
head(df_nonAS)


unique(df$V2)
unique(df_nonAS$V2)

# demo
tmp_trait = "CH"
tsign = t.test(df[(df$V2==tmp_trait) & (df$V3=="E5"), "V4"], df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), "V4"], paired = TRUE)
tsign$p.value
tmp_max = max(c(df[(df$V2==tmp_trait) & (df$V3=="E5"), "V4"], df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), "V4"]))
tmp_p = ifelse(tsign$p.value < 0.05, ifelse(tsign$p.value<0.01, ifelse(tsign$p.value<0.001, "***", "**"),"*"), ".")

ggplot() +
  geom_violin(df[(df$V2==tmp_trait) & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
  geom_quasirandom(df[(df$V2==tmp_trait) & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4, color=V1)) +
  geom_violin(df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
  geom_quasirandom(df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4, color=V1)) +
  geom_text(mapping = aes(x=1.5, y=tmp_max, label = tmp_p), size=6, alpha=ifelse(tsign$p.value>0.05, 0 ,1)) +
  geom_segment(mapping = aes(x=1, xend=2, y=tmp_max*0.98), alpha=ifelse(tsign$p.value>0.05, 0 ,1)) +
  theme_classic()+
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = color_map  # Processing note.
    # Processing note.
  ) +
  theme(
    axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
    axis.title.y = element_text(size=16, color="black"),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y='Fold enrichment', title = tmp_trait)



for (tmp_trait in unique(df_nonAS$V2)) {
  tsign = t.test(df[(df$V2==tmp_trait) & (df$V3=="E5"), "V4"], df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), "V4"], paired = TRUE)
  print(tsign$p.value)
  tmp_max = max(c(df[(df$V2==tmp_trait) & (df$V3=="E5"), "V4"], df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), "V4"]))
  tmp_p = ifelse(tsign$p.value < 0.05, ifelse(tsign$p.value<0.01, ifelse(tsign$p.value<0.001, "***", "**"),"*"), ".")

  ggplot() +
    geom_violin(df[(df$V2==tmp_trait) & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
    geom_quasirandom(df[(df$V2==tmp_trait) & (df$V3=="E5"), ], mapping = aes(x="AS", y=V4, color=V1)) +
    geom_violin(df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4), fill="lightgrey", color=NA, alpha=0.5) +
    geom_quasirandom(df_nonAS[(df$V2==tmp_trait) & (df_nonAS$V3=="E5"), ], mapping = aes(x="non-AS", y=V4, color=V1)) +
    geom_text(mapping = aes(x=1.5, y=tmp_max, label = tmp_p), size=6, alpha=ifelse(tsign$p.value>0.05, 0 ,1)) +
    geom_segment(mapping = aes(x=1, xend=2, y=tmp_max*0.98), alpha=ifelse(tsign$p.value>0.05, 0 ,1)) +
    theme_classic()+
    scale_color_manual(
      # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
      values = color_map  # Processing note.
      # Processing note.
    ) +
    theme(
      axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
      axis.title.y = element_text(size=16, color="black"),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    labs(y='Fold enrichment', title = tmp_trait)
  # ggsave(paste0("260408_AS-ChrState_",tmp_trait,"_AS-E5.png"), width = 5, height = 7.5, units = "cm")
  ggsave(paste0("260408_AS-ChrState_",tmp_trait,"_AS-E5.pdf"), width = 5, height = 7.5, units = "cm")
  

}

# for (qtl in unique(df$V2)
# ) {
#   ggplot(df[df$V2==qtl, ]) +
#     geom_violin(aes(x=V3, y=V4), fill="lightgrey", color=NA, alpha=0.5) +
#     geom_quasirandom(aes(x=V3, y=V4, color=V1)) +
#     theme_classic() +
#     scale_color_manual(
#       # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
# Processing note.
# Processing note.
#     ) +
#     theme(
#       axis.text = element_text(size=14, color="black", angle=0, vjust=0.5, hjust=0.5),
#       axis.title.y = element_text(size=16, color="black"),
#       axis.title.x = element_blank(),
#       legend.position = "none"
#     ) +
#     labs(y='Fold enrichment', title = qtl)
#   ggsave(paste0("260407_",qtl,".png"), width = 12, height = 5.5, units = "cm")
# }







# ##### DOCK5 x AnimalQTL #####
# 
# dir(pattern = "DOCK5")
# 
# reg = read.csv("REG_DOCK5_all", sep = "\t", header = F)
# reg_tsr = read.csv("REG_DOCK5", sep = "\t", header = F)
# head(reg)
# head(reg_tsr)
# # ggplot(reg) +
# #   geom_rect(aes(xmin = V7, xmax = V8, fill=V2, ymin=0,ymax=1)) +
# #   facet_wrap(~V1, ncol = 1) +
# #   theme_void()
# # head(reg)
# unique(reg$V1)
# 
# color_map <- c(
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# )
# 
# re_p = list()
# i = 0
# for (tis in unique(reg$V1)) {
#   i = i+1
#   re_p[[i]] = 
#     ggplot(reg[reg$V1==tis,]) +
#       geom_rect(aes(xmin = V7, xmax = V8, fill=V2, ymin=0,ymax=1), alpha=1) +
#       theme_void() +
#       xlim(40517375, 40519376) +
#       labs(y=tis)+
#       scale_fill_manual(values = color_map) +
#       theme(axis.title.y = element_text(angle = 0, size=14),
#             legend.position = 'none')
#       # geom_vline(xintercept = 40518375, linetype="dotted")
# }
# 
# # re_p[[1]] / re_p[[2]]  / re_p[[3]]  / re_p[[4]]  / re_p[[5]]
# 
# re_p[[3]]
# # paste(paste0("re_p[[", 1:43, "]]"), collapse = " / ")
# cat(paste(paste0("re_p[[", 1:43, "]]"), collapse = " / "))
# 
# unique(reg$V1)[13]
# 
# re_p[[1]] / re_p[[2]] / re_p[[3]] / re_p[[4]] / re_p[[5]] / re_p[[6]] / re_p[[7]] / re_p[[8]] / re_p[[9]] / re_p[[10]] /
#   re_p[[11]] / re_p[[12]] / (re_p[[13]] + 
#                                geom_rect(mapping = aes(xmin=40518200, xmax=40518600, ymin=0, ymax=1, fill="E8")) +
#                                geom_point(mapping = aes(x=40518375, y=0.5), size=3.25, shape=18, color="#70d24b")) / re_p[[14]] / re_p[[15]] / re_p[[16]] / re_p[[17]] / re_p[[18]] / re_p[[19]] / 
#   re_p[[20]] / re_p[[21]] / re_p[[22]] / re_p[[23]] / re_p[[24]] / re_p[[25]] / re_p[[26]] / re_p[[27]] / re_p[[28]] / 
#   re_p[[29]] / re_p[[30]] / re_p[[31]] / re_p[[32]] / re_p[[33]] / re_p[[34]] / re_p[[35]] / re_p[[36]] / re_p[[37]] / 
#   re_p[[38]] / re_p[[39]] / re_p[[40]] / re_p[[41]] / re_p[[42]] / 
#   re_p[[43]] + theme(
#     axis.ticks.x = element_line(),
#     axis.line = element_line(linetype = "dotted"), 
#     # axis.text.x = element_text(size=14, color="black"),
#     axis.ticks.length.x = unit(0.05, "cm"),
#     axis.title.x = element_text(size=16, color="black")) +
#   labs(x="TSR")
# 
# ggsave("260407_REG_DOCK5.png", width = 10.5, height = 20, units = "cm")
# 
# dev.off()
# ggplot()+
#   geom_col(mapping = aes(x=rownames(data.frame(color_map)), y=1, fill=rownames(data.frame(color_map)))) +
#   theme_void() +
#   scale_fill_manual(values = color_map) +
#   theme(legend.position = "top",
#         legend.title = element_blank(),
#         legend.text = element_text(size=14, color="black")
#   )
# ggsave("260407_REG_DOCK5_legend.png", width = 10, height = 2.5, units = "cm")
# 
# 
# 
# ##### DOCK5 GTF #####
# gtf = read.csv("GTF_DOCK5", sep = "\t", header = F)
# head(gtf[,c(1:6)])
# 
# ggplot(gtf) +
#   geom_rect(aes(xmin=V4/1000000, xmax=V5/1000000, ymin=0, ymax=1), fill="darkblue") +
#   geom_vline(xintercept = 40518375/1000000, color="lightgrey", linewidth=2, alpha=0.9) +
#   theme_void() +
#   theme(
#     axis.ticks.x = element_line(),
#     axis.line.x = element_line(), 
#     axis.ticks.length.x = unit(0.05, "cm"),
#     axis.text.x = element_text(size=14, color="black"),
#     axis.title.x = element_text(size=16, color="black")) +
#   labs(x="Chromosome2 (Mb)", title="DOCK5")
# 
# ggplot(gtf[gtf$V3=="exon",]) +
#   geom_rect(aes(xmin=V4/1000000, xmax=V5/1000000, ymin=0, ymax=1), fill="darkblue") +
#   geom_segment(mapping = aes(x=40232384/1000000, xend=40508685/1000000,y=0.5), color="darkblue" ) +
#   geom_vline(xintercept = 40518375/1000000, color="grey", linetype="solid", linewidth=2, alpha=0.9) +
#   
#   theme_void()+
#   theme(
#     axis.ticks.x = element_line(),
#     axis.line.x = element_line(), 
#     axis.ticks.length.x = unit(0.05, "cm"),
#     axis.text.x = element_text(size=14, color="black"),
#     axis.title.x = element_text(size=16, color="black")) +
#   labs(x="Chromosome2 (Mb)", title="DOCK5")
# 
# ggsave("260407_GTF_DOCK5.png", width = 12.5, height = 2.5, units = "cm")
# 
# 
# ##### DOCK5 TPM #####
# 
# dir(pattern = "DOCK5")
# tpm = read.csv("TPM_DOCK5", sep = "\t", header = F)
# head(tpm)
# color_map = c(
#   "cerebral-cortex" = "#dcd71a",
#   "midbrain" = "#f9ed19",
#   "cerebellum" = "#efd80b",
#   "brainstem" = "#f4d578",
#   "hippocampus" = "#f1cb05",
#   "hypothalamus" = "#825e19",
#   "medulla-oblongata" = "#d0b35b",
#   "optic-chiasm" = "#fece01",
#   "pineal" = "#807120",
#   "pituitary" = "#f1d95d",
#   "pons" = "#fcd222",
#   "splenium" = "#7c6919",
#   "rumen" = "#fc9891",
#   "reticulum" = "#d97c68",
#   "omasum" = "#f8ae81",
#   "abomasum" = "#f18264",
#   "duodenum" = "#eb9d63",
#   "jejunum" = "#eb951c",
#   "ileum" = "#ce9639",
#   "cecum" = "#daaa6c",
#   "colon" = "#f2c063",
#   "rectum" = "#efe0a8",
#   "cervix" = "#b6d7a9",
#   "cornua-uteri" = "#69d683",
#   "corpus-uteri" = "#80d897",
#   "ovary" = "#69d28c",
#   "oviduct" = "#79ffaa",
#   "epididymis" = "#70d24b",
#   "testis" = "#7ef351",
#   "mammary-gland" = "#fed9d0",
#   "bone-marrow" = "#d84b4b",
#   "lymph-node" = "#c33a11",
#   "thymus" = "#ff3a32",
#   "thyroid" = "#f359d1",
#   "spleen" = "#962932",
#   "liver" = "#ad8c8b",
#   "kidney" = "#4f3136",
#   "lung" = "#36b5f1",
#   "heart" = "#bc58e3",
#   "muscle" = "#a180ca",
#   "adipose" = "#ffc4f1",
#   "skin" = "#d09dc5",
#   "soft-horn" = "#a25d73"
# ) 
# 
# ggplot(tpm) +
#   geom_col(aes(x=reorder(V1,-V5), y=V5, fill=V1)) +
#   # geom_col(aes(x=reorder(V1,-V5), y=V5)) +
#   
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.ticks.x = element_line(),
#     axis.line.x = element_line(), 
#     axis.ticks.length.x = unit(0.05, "cm"),
#     axis.text = element_text(size=14, color="black",),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     axis.title = element_text(size=16, color="black")) +
#   labs(x="Tissue", y="TPM", title="DOCK5") +
#   scale_fill_manual(
#     limits = c("epididymis"),
# Processing note.
#   ) 
# 
# ggplot(tpm) +
#   geom_col(aes(x=V5, y=V1, fill=V1)) +
#   theme_void() +
#   scale_fill_manual(
#     limits = c("epididymis"),
# Processing note.
#   ) +
#   theme(
#     legend.position = "none",
#     axis.ticks.x = element_line(),
#     axis.line.x = element_line(), 
#     axis.ticks.length.x = unit(0.05, "cm"),
#     axis.text = element_text(size=14, color="black",),
#     # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     axis.title = element_text(size=16, color="black")) +
#   labs(x="TPM", y="", title="DOCK5") +
#   scale_y_discrete(limits = rev) +
#   scale_x_continuous(trans = "reverse")
#   # scale_y_reverse()
# 
# color_map
# 
# ggsave("260407_TPM_DOCK5.png", width = 10, height = 20, units = "cm")
