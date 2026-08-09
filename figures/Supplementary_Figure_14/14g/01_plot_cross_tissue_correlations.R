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


##### X tissues #####
df = read.csv("251129_TisCon_EffectSize", sep = "\t", header = F)
head(df)
unique(df$V1)

# all pairs' effect sizes
cor.test(df$V3, df$V4)
head(df)
ggplot(df) +
  geom_point(aes(x=log(V3), y=log(V4), color=V2), alpha=0.25)+
  geom_smooth(aes(x=log(V3), y=log(V4)), method="lm")+
  theme_classic() +
  theme(
    axis.text = element_text(size=14, color="black"),
    axis.title = element_text(size=16, color="black"),
    plot.title = element_text(size=16, color="black"),
  ) + 
  labs(x="Effect size (log)", y="Effect size (log)", title="All tissue-shared AS SNPs") +
  annotate("text", x=0.1, y=0.4, label="cor ≈ 0.73**", size=5)

# cor.test by Sequencing
head(df)

re_cor = list()
re_graphics = list()
df_cor = data.frame()
i = 0
for (reg in unique(df$V2)) {
  i = i + 1
  re_cor[[i]] = cor.test(df[df$V2==reg,"V3"], df[df$V2==reg,"V4"])
  # cat(pair, re_cor[[i]]$p.value, re_cor[[i]]$estimate, "\n",sep = "\t")
  tmp_ = data.frame(reg=reg, p=re_cor[[i]]$p.value, r2=re_cor[[i]]$estimate)
  df_cor = rbind(df_cor, tmp_)
  
  re_graphics[[i]] = ggplot(df[df$V2==reg,]) +
    geom_point(aes(x=log(V3), y=log(V4)), alpha=0.25, color="grey")+
    geom_smooth(aes(x=log(V3), y=log(V4)), method="lm")+
    theme_classic() +
    theme(
      axis.text = element_text(size=14, color="black"),
      axis.title = element_text(size=16, color="black"),
      plot.title = element_text(size=16, color="black"),
    ) + 
    labs(x="Effect size (log)", y="Effect size (log)", 
         title= paste0("r2=", round(re_cor[[i]]$estimate, 2), " ", "p-val=", re_cor[[i]]$p.value)
    ) 
  # Processing note.
}
rownames(df_cor) = NULL
df_cor

# re_graphics[[1]] / re_graphics[[2]] / re_graphics[[3]] / re_graphics[[4]] / re_graphics[[5]] / re_graphics[[6]] +
# plot_layout(axis_titles = "collect")
re_graphics[[1]] | re_graphics[[2]] 

# # facet_wrap
# ggplot(df) +
#   geom_point(aes(x=log(V3), y=log(V4), color=V2), alpha=0.25)+
#   geom_smooth(aes(x=log(V3), y=log(V4)), method="lm")+
#   theme_classic() +
#   facet_wrap(~V2, ncol = 1, scales = "free") +
#   theme(
#     axis.text = element_text(size=14, color="black"),
#     axis.title = element_text(size=16, color="black"),
#     plot.title = element_text(size=16, color="black"),
#   )

# NUM count
head(df)
df_count = df[,c("V1","V2")] %>%
  count(V1, V2, name = "count")

# head(df_count)
# ggplot(df_count[df_count$count>350,]) +
#   geom_col(aes(x=reorder(V2, count),fill=V1, y=(count)))+
#   # facet_wrap(~V2) +
#   theme_classic() +
#   theme(
#     axis.text = element_text(size=14, color="black"),
#     # axis.text.x = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     legend.title = element_blank(),
#     legend.text = element_text(size=14, color="black"),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size=16, color="black"),
#   ) +
# Processing note.
# Processing note.
# Processing note.
head(df_count)
length(df_count[df_count$count<=3, "count"])
# ggplot(df_count[df_count$count,]) +
#   geom_col(aes(x=reorder(V2, count), y=log10(count), fill=V2)) +
#   theme_classic() +
#   labs(x="", y="Number of the tissue-shared AS SNP (log10)") +
#   theme(
#     legend.position = "none",
#     axis.text = element_text(size=14, color="black"),
#     axis.title = element_text(size=16, color="black"),
#   ) +
#   coord_flip()

# cor
re_cor = list()
df_cor = data.frame()
i = 0
head(df)
for (pair in unique(df$V1)) {
  i = i + 1
  if(length(df[df$V1==pair,"V3"])>=3){
    re_cor[[i]] = cor.test(df[df$V1==pair,"V3"], df[df$V1==pair,"V4"])
    tmp_ = data.frame(tis=pair, p=re_cor[[i]]$p.value, r2=re_cor[[i]]$estimate)
  }
  # else{
  #   print(pair)
  #   tmp_ = data.frame(tis=pair, p=1, r2=0.01)
  # }
  # cat(pair, re_cor[[i]]$p.value, re_cor[[i]]$estimate, "\n",sep = "\t")
  df_cor = rbind(df_cor, tmp_)
}

df_cor$tis_a = str_split_fixed(df_cor$tis, "_", 2)[,1]
df_cor$tis_b = str_split_fixed(df_cor$tis, "_", 2)[,2]
rownames(df_cor) = NULL
head(df_cor[,-1])

unique(df_cor$tis_a)
# Processing note.
# p_matrix <- dcast(df_cor[,-1], tis_a ~ tis_b, value.var = "p")
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

r2_matrix[is.na(r2_matrix)] = 0
p_matrix[is.na(p_matrix)] = 1

# is.infinite(r2_matrix)
# # corrplot
# library(corrplot)
# Processing note.
# corrplot(corr =r2_matrix, p.mat = p_matrix,method = "circle",type = "lower",
#          insig="label_sig",sig.level = c(.05, 0.01),pch.cex = 0.5, pch.col = "black",
#          tl.col="black", tl.cex = 0.7,
#          col = col1(100),
#          order = "AOE",
#          diag = F,
# )

# igraph
library(igraph)
library(ggplot2)
library(reshape2)
set.seed(123)

head(r2_matrix)
head(p_matrix)
# Processing note.
R2_long <- melt(r2_matrix)
P_long <- melt(p_matrix)
# Processing note.
edge_data <- data.frame(
  from = R2_long$Var1,
  to = R2_long$Var2,
  R2 = R2_long$value,
  p_value = P_long$value
)
# Processing note.
edge_data <- edge_data[edge_data$p_value < 0.01 & edge_data$from != edge_data$to, ]
# edge_data <- edge_data[edge_data$p_value < 0.01 & edge_data$from != edge_data$to, ]

# Processing note.
g <- graph_from_data_frame(edge_data, directed = FALSE)
# Processing note.
E(g)$width <- edge_data$R2 * 0.05  # Processing note.
E(g)$color <- ifelse(edge_data$R2 > 0.7, "#DC0000FF", 
                     ifelse(edge_data$R2 > 0.4, "#fed976", "#969696"))


# Processing note.
plot(g,
     vertex.label = V(g)$name,
     vertex.size = 15,
     vertex.color = "lightgrey",
     vertex.frame.color = "gray",
     vertex.label.color = "black",
     edge.curved = 0.25,
     main = "Correlation Network (R² with p < 0.05)",
     layout = layout_with_fr(g))



color_mapping  <- list(
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
V(g)$name
color_mapping
# Processing note.
vertex_colors <- sapply(V(g)$name, function(node_name) {
  if (node_name %in% names(color_mapping)) {
    return(color_mapping[[node_name]])
  } else {
    return("lightgray")  # Processing note.
  }
})
set.seed(75369)
plot(g,
     vertex.label = V(g)$name,
     vertex.shape = "circle", 
     vertex.size = 15,
     vertex.label.cex = 1e-100,  # Processing note.
     # Processing note.
     
     vertex.color = vertex_colors,
     vertex.frame.color = NA,
     
     vertex.label.color = "black",
     edge.curved = 0.2,
     main = NA,
     layout = layout_with_mds(g))



# Processing note.
layouts <- list(
  layout_with_fr(g),  # Processing note.
  layout_with_kk(g),  # Processing note.
  layout_with_dh(g),  # Processing note.
  layout_in_circle(g),  # Processing note.
  layout_as_star(g),  # Processing note.
  layout_as_tree(g),  # Processing note.
  layout_with_lgl(g),  # Processing note.
  layout_with_mds(g)  # Processing note.
)


data.frame(vertex_colors)$vertex_colors
rownames(data.frame(vertex_colors))

ggplot() +
  geom_point(mapping = aes(x=rownames(data.frame(vertex_colors)), y=1, color=rownames(data.frame(vertex_colors))),
             size=4) +
  theme_void() +
  scale_color_manual(
    # limits = c("ATAC", "H3K4me1","RNASeq", "H3K4me3", "H3K27ac","H3K27me3"),
    values = c(
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
    )  # Processing note.
    # Processing note.
  ) +
  guides(color = guide_legend(ncol = 5)) +
  theme(
    legend.position = "left",
    legend.text = element_text(size=14 ,color='black'),
    legend.title = element_blank()
  )+
  coord_fixed()

ggsave("tissue_color_legend.png", width = 45, height = 25, units = 'cm')

# # pheatmap
# library(pheatmap)
# pheatmap(r2_matrix)
# ggplot(df[df$V1=="pons_testis",]) +
#   geom_point(aes(x=log(V3), y=log(V4))) +
#   # facet_wrap(~V1, nrow = 1)+
#   theme_classic()
