library(ggplot2)
setwd('/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/')
data <- read.csv("all_super_enhancer_summary_1.txt", sep = "\t", header = T) 

data$Tissue = factor(data$Tissue, levels=c("abomasum", "adipose", "bone-marrow", "brainstem", "cecum", "cerebellum", "cerebral-cortex", "cervix", "colon", "cornua-uteri", "corpus-uteri", "duodenum", "epididymis", "heart", "hippocampus", "hypothalamus", "ileum", "jejunum", "kidney", "liver", "lung", "lymph-node", "mammary-gland", "medulla-oblongata", "midbrain", "muscle", "omasum", "optic-chiasm", "ovary", "oviduct", "pineal", "pituitary", "pons", "rectum", "reticulum", "rumen", "skin", "soft-horn", "spleen", "splenium", "testis", "thymus", "thyroid"))

data$gap=data$Super_average_size-data$Origin_size
colnames(data)
pdf("all_super_enhancer_summary_1.pdf", width=2.7, height=5)
makers =c("Number","Origin_size", "Super_size","Super_average_size","gap")
for (maker in makers) {
  print(ggplot(data=data, aes(x = Tissue, y = !!sym(maker), fill=Tissue)) +
          geom_col(position=position_dodge(1), width=0.8 )+
          ylab(paste(maker, " of super enhancer                               ", sep=""))+xlab("")+
          theme_bw()+theme(legend.position = "none")+
          scale_x_discrete(limits = rev(levels(data$Tissue)))+
          coord_flip()+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          theme(axis.text.x=element_text(angle=30, hjust = 1,colour="black",family="Times",size=10),
                axis.text.y=element_text(family="Times",size=10,face="plain"),
                axis.title.y=element_text(family="Times",size = 15,face="plain"),
                axis.title.x=element_text(family="Times",size = 10,face="plain"))+
          scale_fill_manual(values=c("abomasum" = "#f18264",
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
      "thyroid" = "#f359d1"))
  )      
}

dev.off()