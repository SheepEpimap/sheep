#PCA
library(ggplot2)
setwd('/vol2/mengzhu/snakemake_sheep/Figures/') 
makers =c("RNASeq", "ATAC", "H3K27ac","H3K27me3","H3K4me3","H3K4me1")
pdf("PCA_all_Assay.pdf", width=7.4, height=5)
for(maker in makers)
{
data <- read.table(paste(maker, "_PCA_last.txt", sep=""), header = T) 
  # Convert the Rep column to a factor
  data$Rep <- as.factor(data$Rep)

print(ggplot(data)+geom_point(aes(x=X1,y=X2,color= group, shape = Rep),size=3)+
        theme(legend.title =element_blank())+labs(x="PCA1",y="PCA2")+labs(title = maker)+
        scale_shape_manual(values=c(16, 17, 12, 10))+
  scale_color_manual(values=c("abomasum" = "#f18264",
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
      "thyroid" = "#f359d1"))+
        theme(axis.text.x=element_text(colour="black",family="Times",size=15), #Set the x-axis tick-label font, rotate it by 15 degrees, shift it down by 1 (hjust = 1), and use Times at size 20
              axis.text.y=element_text(family="Times",size=15,face="plain"), #Set the y-axis tick-label font family, size, and plain style
              axis.title.y=element_text(family="Times",size = 15,face="plain"),
              axis.title.x=element_text(family="Times",size = 15,face="plain"))+#Set the y-axis title font properties
        theme(legend.text=element_text(family="Times", colour="black",  #Set the legend-label font properties
                                       size=16))+
        theme(legend.title=element_text(family="Times", colour="black", #Set the legend-title font properties
                                        size=16))+theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_hline(yintercept = 0, colour="blue", linetype="dashed" )+
    geom_vline(xintercept = 0, colour="blue", linetype="dashed" )
  
)
}
dev.off() 
