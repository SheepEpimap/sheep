setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(UpSetR)
library(tools)  # Processing note.
# Processing note.
setwd("./CrossGroup")

files <- list.files(pattern = "", full.names = TRUE)

# Processing note.
data_list <- lapply(files, function(x) {
  read.csv(x, sep = "\t", header = F)[,c(1,2)]
})

# Processing note.
# names(data_list) <- file_path_sans_ext(basename(files))
names(data_list) <- gsub("by|Omics|Tis","",file_path_sans_ext(basename(files)))

# Processing note.
# Processing note.
head(data_list$Ind39)
head(paste0(data_list$Ind39$V1, ":", data_list$Ind39$V2))

gsub("by|Omics|Tis","",file_path_sans_ext(basename(files)))

head(data_list[[1]],2)
for (i in c(1:length(data_list))) {
  print(i)
  print(head(data_list[[i]],5))
  data_list[[i]]$ID = paste0(data_list[[i]]$V1,":",data_list[[i]]$V2)
}
head(data_list[[1]],2)


# Processing note.

# Processing note.
group_names <- names(data_list)

# Processing note.
intersection_matrix <- matrix(0, 
                              nrow = length(group_names), 
                              ncol = length(group_names),
                              dimnames = list(group_names, group_names))


# Processing note.
# Processing note.

# Processing note.
for(i in seq_along(group_names)) {
  for(j in seq_along(group_names)) {
    if(i == j) {
      # Processing note.
      intersection_matrix[i, j] <- length(data_list[[i]]$ID)
    } else if(i < j) {
      # Processing note.
      inter <- intersect(data_list[[i]]$ID, data_list[[j]]$ID)
      intersection_matrix[i, j] <- length(inter)
      intersection_matrix[j, i] <- length(inter)  # Processing note.
    }
  }
}

# Processing note.
print(intersection_matrix)

# Processing note.
intersection_df <- as.data.frame(intersection_matrix)
intersection_df

# ?pheatmap
library(pheatmap)

pheatmap(intersection_df,
         cluster_rows=F, cluster_cols=F,
         cellwidth=10, cellheight = 10,
         )

# proportion?
intersection_matrix_prop =
  matrix(0, nrow = length(group_names), 
         ncol = length(group_names),
         dimnames = list(group_names, group_names))
# Processing note.
for(i in seq_along(group_names)) {
  for(j in seq_along(group_names)) {
    if(i == j) {
      # Processing note.
      intersection_matrix_prop[i, j] = 1
    } else if(i < j) {
      # Processing note.
      inter <- intersect(data_list[[i]]$ID, data_list[[j]]$ID)
      # Processing note.
      total_num = length(unique(c(data_list[[i]]$ID, data_list[[j]]$ID)))
      intersection_matrix_prop[i, j] <- length(inter)/total_num
      intersection_matrix_prop[j, i] <- length(inter)/total_num
    }
  }
}


# Processing note.
files
substring(files, 5, 7)
annotation_col = data.frame(
  group = substring(files, 5, 7)
)
row.names(annotation_col) <- group_names

intersection_matrix_prop

# Processing note.
for(i in 1:nrow(intersection_matrix_prop)) {
  intersection_matrix_prop[i, i] <- 0
}
# print(df)
(annotation_col)

annotation_colors <- list(
  group = c(
    Ind39="white",
    Ind40 = "white",
    ATAC = "#f768a1",  # Processing note.
    H3K27ac = "#DC0000FF",  # Processing note.
    H3K4me3 = "#41ab5d",  # Processing note.
    H3K27me3 = "#969696",  # Processing note.
    H3K4me1 ="#fed976",  # Processing note.
    RNASeq ="#41b6c4",
    `cerebral-cortex` = "#dcd71a",
    midbrain = "#f9ed19",
    cerebellum = "#efd80b",
    brainstem = "#f4d578",
    hippocampus = "#f1cb05",
    hypothalamus = "#825e19",
    `medulla-oblongata` = "#d0b35b",
    `optic-chiasm` = "#fece01",
    pineal = "#807120",
    pituitary = "#f1d95d",
    pons = "#fcd222",
    splenium = "#7c6919",
    rumen = "#fc9891",
    reticulum = "#d97c68",
    omasum = "#f8ae81",
    abomasum = "#f18264",
    duodenum = "#eb9d63",
    jejunum = "#eb951c",
    ileum = "#ce9639",
    cecum = "#daaa6c",
    colon = "#f2c063",
    rectum = "#efe0a8",
    cervix = "#b6d7a9",
    `cornua-uteri` = "#69d683",
    `corpus-uteri` = "#80d897",
    ovary = "#69d28c",
    oviduct = "#79ffaa",
    epididymis = "#70d24b",
    testis = "#7ef351",
    `mammary-gland` = "#fed9d0",
    `bone-marrow` = "#d84b4b",
    `lymph-node `= "#c33a11",
    thymus = "#ff3a32",
    thyroid = "#f359d1",
    spleen = "#962932",
    liver = "#ad8c8b",
    kidney = "#4f3136",
    lung = "#36b5f1",
    heart = "#bc58e3",
    muscle = "#a180ca",
    adipose = "#ffc4f1",
    skin = "#d09dc5",
    `soft-horn` = "#a25d73"
    )  
)

# annotation_colors <- list(
#   group = c(
#     "Ind39"="white",
#     "Ind40" = "white",
# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
#     "RNASeq" ="#41b6c4",
#     "cerebral-cortex" = "#dcd71a",
#     "midbrain" = "#f9ed19",
#     "cerebellum" = "#efd80b",
#     "brainstem" = "#f4d578",
#     "hippocampus" = "#f1cb05",
#     "hypothalamus" = "#825e19",
#     "medulla-oblongata" = "#d0b35b",
#     "optic-chiasm" = "#fece01",
#     "pineal" = "#807120",
#     "pituitary" = "#f1d95d",
#     "pons" = "#fcd222",
#     "splenium" = "#7c6919",
#     "rumen" = "#fc9891",
#     "reticulum" = "#d97c68",
#     "omasum" = "#f8ae81",
#     "abomasum" = "#f18264",
#     "duodenum" = "#eb9d63",
#     "jejunum" = "#eb951c",
#     "ileum" = "#ce9639",
#     "cecum" = "#daaa6c",
#     "colon" = "#f2c063",
#     "rectum" = "#efe0a8",
#     "cervix" = "#b6d7a9",
#     "cornua-uteri" = "#69d683",
#     "corpus-uteri" = "#80d897",
#     "ovary" = "#69d28c",
#     "oviduct" = "#79ffaa",
#     "epididymis" = "#70d24b",
#     "testis" = "#7ef351",
#     "mammary-gland" = "#fed9d0",
#     "bone-marrow" = "#d84b4b",
#     "lymph-node" = "#c33a11",
#     "thymus" = "#ff3a32",
#     "thyroid" = "#f359d1",
#     "spleen" = "#962932",
#     "liver" = "#ad8c8b",
#     "kidney" = "#4f3136",
#     "lung" = "#36b5f1",
#     "heart" = "#bc58e3",
#     "muscle" = "#a180ca",
#     "adipose" = "#ffc4f1",
#     "skin" = "#d09dc5",
#     "soft-horn" = "#a25d73"
#   )  
# )


pheatmap(intersection_matrix_prop,
         cluster_rows=F, cluster_cols=F,
         cellwidth=10, cellheight = 10,
         color = colorRampPalette(colors = c("white","#DC0000FF"))(100),
         annotation_row = annotation_col, annotation_col = annotation_col,
         gaps_row = c(2,8), gaps_col = c(2,8),
         annotation_colors = list(group=c(Ind="#f768a1", Omi="#fed976", Tis="#41b6c4")),
         filename = "../AS_SNP_overlap_between_groups.pdf"
)
dev.off()


exit(0)



##### ######
intersection_matrix_prop[c(1,2),c(1,2)]

df_matrix <- as.matrix(intersection_matrix_prop[c(3:8),c(3:8)])
diag(df_matrix) <- NA  # Processing note.
mean(df_matrix, na.rm = TRUE)

df_matrix <- as.matrix(intersection_matrix_prop[c(9:51),c(9:51)])
diag(df_matrix) <- NA  # Processing note.
mean(df_matrix, na.rm = TRUE)

# 
df_matrix <- as.matrix(intersection_matrix_prop[c(3:8),c(1:2)])
diag(df_matrix) <- NA  # Processing note.
mean(df_matrix, na.rm = TRUE)

df_matrix <- as.matrix(intersection_matrix_prop[c(9:51),c(1:2)])
diag(df_matrix) <- NA  # Processing note.
mean(df_matrix, na.rm = TRUE)

df_matrix <- as.matrix(intersection_matrix_prop[c(9:51),c(3:8)])
diag(df_matrix) <- NA  # Processing note.
mean(df_matrix, na.rm = TRUE)


mean_cor = 
  data.frame(ind=c( 0.2073699, 0.136281, 0.04042218),
           omics=c(0.136281,  0.01600266, 0.02310529), 
           tis=c(0.04042218,0.02310529,0.06448509))
rownames(mean_cor) = c("ind", "omics", "tis")
mean_cor
pheatmap(mean_cor,
         cluster_rows=F, cluster_cols=F,
         cellwidth=20, cellheight = 20,
         color = colorRampPalette(colors = c("white","pink"))(10),
         gaps_row = c(1,2), gaps_col = c(1,2)
         )



# ?corrplot
library(corrplot)
corrplot(as.matrix(mean_cor), 
         type = "lower",
         tl.col = 'black',
         col.lim = c(0, 1), 
         addCoef.col = "black",
         number.cex = 1.5,  # Processing note.
         col = colorRampPalette(c( "white", "pink"))(100),
         # Processing note.
         )

# data("mtcars")
# Processing note.
# M <- cor(mtcars)
# class(M)
# M
# dev.off()
# Processing note.
# corrplot(M, method = "circle",
#          number.cex = 1,
#          addCoef.col="black"
# )
# str(M)


library(reshape2)  # Processing note.
# Processing note.
df_melt <- melt(as.matrix(mean_cor))
head(df_melt)

ggplot(df_melt, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +



##### UpsetR #####
listInput <- list(
  ATAC = paste0(data_list$ATAC$V1, ":", data_list$ATAC$V2),
  H3K4me3 = paste0(data_list$H3K4me3$V1, ":", data_list$H3K4me3$V2),
  cervix = paste0(data_list$cervix$V1, ":", data_list$cervix$V2),
  
  Ind40 = paste0(data_list$Ind40$V1, ":", data_list$Ind40$V2),
  Ind39 = paste0(data_list$Ind39$V1, ":", data_list$Ind39$V2)
)
upset(fromList(listInput), order.by = "freq")


# Processing note.
group_names <- names(data_list)
