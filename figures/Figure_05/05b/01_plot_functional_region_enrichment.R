
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(stringr)
library(dplyr)
library(patchwork)



df = read.csv("251211_func_rate", sep = "\t", header = F)

head(df)


# ggplot(df) +
#   geom_point(aes(x=V2, y=V3/V4, color=V1)) +
#   # facet_wrap(~V2, scales = "free_x", nrow=1) +
#   theme_classic()
# 
# ggplot(df[df$V2=="intergenic",]) +
#   geom_col(aes(x=paste0(V2,"_",V1), y=V3/V4, fill=V1), width=0.025) +
#   geom_point(aes(x=paste0(V2,"_",V1), y=V3/V4, color=V1), size=3) +
#   theme_classic() +
#   theme(
#     legend.position = "none"
#   )
# ggplot(df) +
#   geom_point(aes(x=paste0(V2,"_",V1), y=V3/V4, color=V1)) +
#   facet_wrap(~V2, scales = "free", nrow=1) +
#   theme_classic()


head(df)
df$rate = df$V3 / df$V4
head(df)
str(df)

library(tidyr)
library(dplyr)

# Processing note.
df_complete <- df %>%
  complete(V1, V2,
           fill = list(V3 = 0, V4 = 0, rate = 0)) %>%
  arrange(V1, V2)


# head(df_complete)
# p_ls = list()
# i = 0
# for (func in unique(df_complete$V2)) {
#   i = i + 1
#   print(func)
#   p_ls[[i]] = 
#     ggplot(df_complete[df_complete$V2==func,]) +
#       geom_col(aes(x=paste0(V1,"_",V2), y=V3/V4, fill=V1), width=0.025) +
#       geom_point(aes(x=paste0(V1,"_",V2), y=V3/V4, color=V1), size=3) +
#       theme_classic() +
#       labs(title=func) +
#       theme(
#         legend.position = "none"
#       )
#       labs(title=func)
# }
# 
# p_ls[[1]] | p_ls[[2]]
# p_ls[[8]]



# fold_enrichment
unique(df_complete$V1)
df_AllSNP = df_complete[df_complete$V1=="AllSNP",]
df_complete = merge(df_complete, df_AllSNP[,c("V2","rate")], by="V2" , all = T)
df_complete$fold_enrichment = df_complete$rate.x/df_complete$rate.y

head(df_complete)

# ggplot(df_complete[df_complete$V1!="AllSNP",]) +
#   geom_col(aes(x=paste0(V2,"_",V1), y=fold_enrichment, fill=V1), width=0.025) +
#   geom_point(aes(x=paste0(V2,"_",V1), y=fold_enrichment, color=V1), size=3) +
#   geom_hline(yintercept = 1, linetype=2) + 
#   # facet_wrap(~V2, nrow = 1, scales = "free_x") + 
#   theme_classic() 
#   # coord_flip()
# dev.off()
p_ls = list()
i = 0
# c("")
# for (func in unique(df_complete$V2)) {
head(df_complete)
unique(df_complete$V1)
for (func in c("intergenic", "upstream", "downstream", "exonic", "intronic", "UTR5", "UTR3")) {
  i = i + 1
  print(func)
  print(i%%2==0)
  p_ls[[i]] =
    ggplot(df_complete[(df_complete$V2==func) & (df_complete$V1!="AllSNP") & (df_complete$V1!="RNASeq"),]) +
      geom_rect(mapping = aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), 
                fill=ifelse(i%%2!=0, "lightgrey", "white"), alpha=0.05) + 
      geom_col(aes(y=V1, x=log(fold_enrichment), fill=V1), width=0.05) +
      geom_point(aes(y=V1, x=log(fold_enrichment), color=V1), size=3) +
      geom_vline(xintercept = 0, linetype=2) + 
      theme_void() +
        labs(y=func) +
      scale_x_continuous(limits = c(-2.5,5)) +
      # coord_flip() +
      theme(
        axis.title.y = element_text(size = 14, color="black", angle = 0, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        legend.position = "none"
      ) +
    scale_fill_manual(
      limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
      values = c(
        "ATAC" = "#f768a1",  # Processing note.
        "H3K27ac" = "#DC0000FF",  # Processing note.
        "H3K4me3" = "#41ab5d",  # Processing note.
        "H3K27me3" = "#969696",  # Processing note.
        "H3K4me1" ="#fed976",  # Processing note.
        "RNASeq"="#41b6c4")  # Processing note.
      # Processing note.
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
}
# p_ls[[14]] + theme(axis.text = element_text(color="black", size = 10))
# unique(df_complete$V2)

# p_ls[[1]]  / p_ls[[2]] / p_ls[[4]] /  p_ls[[5]] /  p_ls[[11]]  /  p_ls[[13]] /  p_ls[[14]] 
dev.off()
p_ls[[1]]  / p_ls[[2]] / p_ls[[3]] /  p_ls[[4]] /  p_ls[[5]]  /  p_ls[[6]] /  p_ls[[7]] + 
  theme(axis.text.x = element_text(size=14, color="black"),
        axis.title.x =  element_text(size=16, color="black")) +
  labs(x="Fold enrichment (log)")
# ggsave("functional_region_enrichment.png", unit="cm", width=10, height=15)

# ggsave("functional_region_enrichment.png", width = 11, height = 17, units = "cm")
ggsave("functional_region_enrichment.pdf")

p_ls[[6]] + 
  theme(
    legend.position = "right",
    legend.text = element_text(size=14, color="black"),
    legend.title = element_blank()) 

# ggplot(df_complete[(df_complete$V2=="exonic") & (df_complete$V1!="AllSNP"),]) +
#   geom_rect(mapping = aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), fill="lightgrey", alpha=0.05) + 
#   geom_col(aes(y=V1, x=log(fold_enrichment), fill=V1), width=0.025) +
#   geom_point(aes(y=V1, x=log(fold_enrichment), color=V1), size=3) +
#   geom_vline(xintercept = 0, linetype=2) + 
#   theme_void() +
#   labs(y="exonic") +
#   scale_x_continuous(limits = c(-2.5,5)) +
#   # coord_flip() +
#   theme(
#     axis.title.y = element_text(size = 14, color="black", angle = 0, vjust = 0.5, hjust = 1),
#     axis.title.x = element_blank(),
#     axis.text = element_blank(),
#     legend.position = "none"
#   ) 

(max(df_complete$fold_enrichment))
(min(df_complete$fold_enrichment))

log(1)
log(0.5)
log(2)

