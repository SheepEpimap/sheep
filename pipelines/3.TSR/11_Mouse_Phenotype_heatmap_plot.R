library(ComplexHeatmap) 
library(circlize)
setwd('/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/GO_GREAT/Gorich/') #mac
#setwd('/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/GO_GREAT/HPrich/') #mac
#setwd('/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/GO_GREAT/MPrich/') #mac
data <- read.csv('TSR_go_enhancer_enrichment.csv', sep = "\t", header = T) 
#data <- read.csv('TSR_go_HP_enrichment.csv', sep = "\t", header = T) 
#data <- read.csv('TSR_go_MP_enrichment.csv', sep = "\t", header = T) 
colnames(data) <- gsub("\\.", "-", colnames(data))  
data1 = data[,c(4:46)]
rownames(data1)=data[,2] 
data1 <- data1[rowSums(data1)>0,]
split = factor(colnames(data1), levels= c("abomasum", "adipose", "bone-marrow", "brainstem", "cecum", "cerebellum", "cerebral-cortex", "cervix", "colon", "cornua-uteri", "corpus-uteri", "duodenum", "epididymis", "heart", "hippocampus", "hypothalamus", "ileum", "jejunum", "kidney", "liver", "lung", "lymph-node", "mammary-gland", "medulla-oblongata", "midbrain", "muscle", "omasum", "optic-chiasm", "ovary", "oviduct", "pineal", "pituitary", "pons", "rectum", "reticulum", "rumen", "skin", "soft-horn", "spleen", "splenium", "testis", "thymus", "thyroid"))
colnames(data1)= factor(colnames(data1), levels=c("abomasum", "adipose", "bone-marrow", "brainstem", "cecum", "cerebellum", "cerebral-cortex", "cervix", "colon", "cornua-uteri", "corpus-uteri", "duodenum", "epididymis", "heart", "hippocampus", "hypothalamus", "ileum", "jejunum", "kidney", "liver", "lung", "lymph-node", "mammary-gland", "medulla-oblongata", "midbrain", "muscle", "omasum", "optic-chiasm", "ovary", "oviduct", "pineal", "pituitary", "pons", "rectum", "reticulum", "rumen", "skin", "soft-horn", "spleen", "splenium", "testis", "thymus", "thyroid"))

type = colnames(data1)
ha = HeatmapAnnotation(tissue = type,annotation_name_side = "right",annotation_legend_param = list(at = colnames(data1), labels= c("abomasum", "adipose", "bone-marrow", "brainstem", "cecum", "cerebellum", "cerebral-cortex", "cervix", "colon", "cornua-uteri", "corpus-uteri", "duodenum", "epididymis", "heart", "hippocampus", "hypothalamus", "ileum", "jejunum", "kidney", "liver", "lung", "lymph-node", "mammary-gland", "medulla-oblongata", "midbrain", "muscle", "omasum", "optic-chiasm", "ovary", "oviduct", "pineal", "pituitary", "pons", "rectum", "reticulum", "rumen", "skin", "soft-horn", "spleen", "splenium", "testis", "thymus", "thyroid")),
                       col = list(tissue =c("abomasum" = "#f18264",
      "adipose" = "#ffc4f1",
      "bone-marrow" = "#d84b4b",
      "brainstem" = "#f4d578",
      "cecum" = "#daaa6c",
      "cerebellum" = "#efd80b",
      "cerebral-cortex" = "#dcd71a",
      "cervix" = "#b6d7a9",
      "colon" = "#f2c063",
      "cornua-uteri" = "#69d683",
      "corpus-uteri" = "#80d897",
      "duodenum" = "#eb9d63",
      "epididymis" = "#70d24b",
      "heart" = "#bc58e3",
      "hippocampus" = "#f1cb05",
      "hypothalamus" = "#825e19",
      "ileum" = "#ce9639",
      "jejunum" = "#eb951c",
      "kidney" = "#4f3136",
      "liver" = "#ad8c8b",
      "lung" = "#36b5f1",
      "lymph-node" = "#c33a11",
      "mammary-gland" = "#fed9d0",
      "medulla-oblongata" = "#d0b35b",
      "midbrain" = "#f9ed19",
      "muscle" = "#a180ca",
      "omasum" = "#f8ae81",
      "optic-chiasm" = "#fece01",
      "ovary" = "#69d28c",
      "oviduct" = "#79ffaa",
      "pineal" = "#807120",
      "pituitary" = "#f1d95d",
      "pons" = "#fcd222",
      "rectum" = "#efe0a8",
      "reticulum" = "#d97c68",
      "rumen" = "#fc9891",
      "skin" = "#d09dc5",
      "soft-horn" = "#a25d73",
      "spleen" = "#962932",
      "splenium" = "#7c6919",
      "testis" = "#7ef351",
      "thymus" = "#ff3a32",
      "thyroid" = "#f359d1")))
pdf("/vol2/mengzhu/jupyter/figure/TSR_go_GO_enrichment.pdf", width = 12, height = 8)
#pdf("/vol2/mengzhu/jupyter/figure/TSR_go_HP_enrichment.pdf", width = 12, height = 8)
#pdf("/vol2/mengzhu/jupyter/figure/TSR_go_MP_enrichment.pdf", width = 12, height = 8)
Heatmap(data1, border = TRUE, show_column_names = F,show_row_names = F, column_gap = unit(0, "mm"), cluster_column_slices = FALSE, column_title =NULL, column_split =split, bottom_annotation = ha, row_names_gp = gpar(fontsize = 7), cluster_rows = FALSE, cluster_columns = FALSE,col = colorRamp2(c(0, 0, 20),c("white", "white", "red"))) 
dev.off()
