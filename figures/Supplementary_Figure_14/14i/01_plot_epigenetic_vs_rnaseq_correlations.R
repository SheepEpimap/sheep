R.version
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(ggplot2)
library(stringr)
dir()
library(CMplot)
library(dplyr)

library(ggchicklet)

library(patchwork)

library(reshape2)
library(tibble)   # column_to_rownames

# Processing note.

df = read.csv("251129_SeqCon_EffectSize", sep = "\t", header = F)


head(df)
unique(df$V1)

head(df[grepl("RNASeq",df$V1),])
unique(df[grepl("RNASeq",df$V1), "V1"])

cor.test(log(df[grepl("RNASeq",df$V1), "V3"]), log(df[grepl("RNASeq",df$V1), "V4"]))

ggplot(df[grepl("RNASeq",df$V1),]) +
  geom_point(aes(x=log(V3), y=log(V4), color=str_split_fixed(V1,"_",2)[,1]), size=3) +
  theme_classic() +
  geom_smooth(aes(x=log(V3), y=log(V4)), method = "lm") +
  labs(x="Effect size (log): Epigenetics" , y="Effect size (log): RNA-seq") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = c(0.7,0.3),
    axis.text = element_text(size=14, color="black"),
    axis.title = element_text(size=16, color="black"),
    plot.title = element_text(size=16, color="black"),
  ) +
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    limits = c("ATAC","H3K4me3", "H3K27ac"),
    values = c(
      "ATAC" = "#f768a1",  # Processing note.
      "H3K27ac" = "#DC0000FF",  # Processing note.
      "H3K4me3" = "#41ab5d",  # Processing note.
      "H3K27me3" = "#969696",  # Processing note.
      "H3K4me1" ="#fed976",  # Processing note.
      "RNASeq"="#41b6c4")  # Processing note.
    # Processing note.
  ) 
  # annotate("text", x = -0.5, y = 2, label = "Some text")
  # ggtitle("cor = 0.543; p-value = 0.001295")
  # Processing note.

ggsave("EffectSize_RNA-Epigenetics.png", width = 11, height = 8.5, units = "cm")
ggsave("EffectSize_RNA-Epigenetics.pdf")

# Processing note.
re_cor = list()
df_cor = data.frame()
i = 0
for (pair in unique(df$V1)) {
  i = i + 1
  re_cor[[i]] = cor.test(df[df$V1==pair,"V3"], df[df$V1==pair,"V4"])
  # cat(pair, re_cor[[i]]$p.value, re_cor[[i]]$estimate, "\n",sep = "\t")
  tmp_ = data.frame(tis=pair, p=re_cor[[i]]$p.value, r2=re_cor[[i]]$estimate)
  df_cor = rbind(df_cor, tmp_)
}
re_cor[[1]]$p.value
re_cor[[1]]$estimate
rownames(df_cor) = NULL
head(df_cor)

# Processing note.
i = 0
ls_cor_geompoint = list()
for (pair in unique(df$V1)) {
  print(pair)
  i = i + 1
  tmp_x_title = as.vector(str_split_fixed(pair,"_",2)[,1])
  tmp_y_title = as.vector(str_split_fixed(pair,"_",2)[,2])

  if(re_cor[[i]]$p.value < 0.01){
    tmp_sign = "**"
  } else if (re_cor[[i]]$p.value < 0.05){
    tmp_sign = "*"
  } else {
    tmp_sign = ""
  }

  ls_cor_geompoint[[i]] =
    ggplot(df[df$V1==pair,]) +
    geom_point(aes(x=log(V3), y=log(V4)), alpha=0.5) +
    geom_smooth(aes(x=log(V3), y=log(V4)), method="lm") +
    # annotate("text", x=log(max(df[df$V1==pair,"V3"]))*0.8, y=log(max(df[df$V1==pair,"V4"]))*0.8,
    # Processing note.
    theme_classic()+
    labs(x=paste0("Effect size (",tmp_x_title,")"), y=paste0("Effect size (",tmp_y_title,")"),
         title = paste0("cor ≈ ", round(re_cor[[i]]$estimate,2), tmp_sign)
    )+
    theme(
      axis.text = element_text(size=14, color="black"),
      axis.title = element_text(size=16, color="black"),
    )

}

unique(df$V1)
complex_design <- "AB##C"
(ls_cor_geompoint[[5]] | ls_cor_geompoint[[8]] | ls_cor_geompoint[[10]]) +
  plot_layout(design = complex_design)


df_cor$tis_a = str_split_fixed(df_cor$tis, "_", 2)[,1]
df_cor$tis_b = str_split_fixed(df_cor$tis, "_", 2)[,2]

# Processing note.
df_complete <- df_cor[,-1] %>%
  bind_rows(df_cor %>% 
              rename(tis_a = tis_b, tis_b = tis_a) %>%
              select(tis_a, tis_b, p, r2)) %>%
  distinct(tis_a, tis_b, .keep_all = TRUE)
# Processing note.

p_matrix <- dcast(df_complete, tis_a ~ tis_b, value.var = "p") %>%
  column_to_rownames("tis_a") %>%
  as.matrix()
r2_matrix <- dcast(df_complete, tis_a ~ tis_b, value.var = "r2") %>%
  column_to_rownames("tis_a") %>%
  as.matrix()
# Processing note.
diag(p_matrix) <- 0  # Processing note.
diag(r2_matrix) <- 1  # Processing note.

head(p_matrix)
head(r2_matrix)

p_matrix[is.na(p_matrix)] = 1
r2_matrix[is.na(r2_matrix)] = 0

# library(ggcorrplot)
# ggcorrplot(r2_matrix,
#            p.mat = p_matrix, hc.order = TRUE,
#            type = "lower", 
#            insig = "blank"
# )

library(corrplot)
col1=colorRampPalette(colors =c("blue","white","red"),space="Lab")  # Processing note.
corrplot(corr =r2_matrix, p.mat = p_matrix,method = "circle",type = "lower",
         insig="label_sig",sig.level = c(.01, .05),pch.cex = 2, pch.col = "black",
         tl.col="black", tl.cex = 1.5,
         col = col1(100),
         # order = "AOE",
         diag = F,
)

