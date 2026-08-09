
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(stringr)
dir()
library(CMplot)
library(dplyr)

library(ggchicklet)

library(patchwork)

##### AS SNP num X chr length (each chr) #####
df = read.csv("251126_SNPNumOnEachChr", sep = "\t", header = F)
head(df)
#
tmp_group_by_seq = df %>%
  group_by(across(c(V1))) %>%
  summarise(sum = sum(V3))
head(tmp_group_by_seq)

df = merge(df, tmp_group_by_seq, by="V1")
df$prop = df$V3 / df$sum
# #
# ggplot(df) +
#   geom_point(aes(x=log(V4), y=(prop), color=V1, size=prop)) +
#   facet_wrap(~V1, nrow=1, scales = "free") +
#   theme_classic() +
#   labs(x="Chr length", y="")


head(df)
re_cor = list()
i=0
for (reg in unique(df$V1)) {
  i = i+1
  re_cor[[i]] = cor.test(df[df$V1==reg,"V4"], df[df$V1==reg,"prop"], method="pearson")
  # print(reg)
  # print(re_cor[[i]]$p.value)
  # print(re_cor[[i]]$estimate)
  print(paste(reg, round(as.numeric(re_cor[[i]]$estimate),2), (re_cor[[i]]$p.value), sep = "    "))
}

i=0
re_cor_graphics = list()
for (reg in unique(df$V1)) {
  print(reg)
  i = i+1

  if(as.numeric(re_cor[[i]]$p.value)<0.01){
    tmp_p_val = "P<0.01**"
  }else if(as.numeric(re_cor[[i]]$p.value)<0.05){
    tmp_p_val = "P<0.05*"
  } else {
    tmp_p_val = "P>0.05"
  }
  re_cor_graphics[[i]] =
      ggplot(df[df$V1==reg,]) +
        geom_point(aes(x=(V4)/1000000, y=(prop)*100, color=V1)) +
        geom_smooth(aes(x=(V4)/1000000, y=(prop)*100), method = "lm") +
        theme_classic() +
        labs(x="Chr length (Mb)", y="Proportion (%)", title = reg) +
        annotate('text', x=100, y=max(df[df$V1==reg,"prop"])*85, 
                 label=paste0("R2=",round(as.numeric(re_cor[[i]]$estimate),2),"  ", tmp_p_val), size=5
                 # label=expression(-log[10]*'(p value)'),size=8,color='red'
                 ) +  # Processing note.
         theme(
           legend.position = "none",
            axis.text = element_text(size=14, color="black"),
            axis.title = element_text(size=16, color="black"),
            plot.title = element_text(size=16, color="black"),
         )

}

re_cor_graphics[[1]] / re_cor_graphics[[4]] / re_cor_graphics[[6]] / 
  re_cor_graphics[[5]]/ re_cor_graphics[[2]]/ re_cor_graphics[[3]] + plot_layout(axis_titles = "collect")


# Processing note.
ggsave("260408_AS-SNPinformation.pdf")

# ###### ##### 
# tmp = read.csv("251125_AS_Dis", sep = "\t", header = F)
# head(tmp)
# 
# Processing note.
# each_omic_count = tmp %>% group_by(V1) %>%
#   mutate(count=n())  %>%
#   ungroup() %>%
#   select(V1, count) %>%
#   unique()
# 
# each_omic_count[each_omic_count(each_omic_count$count,decreasing = T),]


# Processing note.
df = read.csv("251205_SNPNumEachSeqAndTis", sep = "\t", header = F)
head(df)

library(ggalluvial)

color_map = c(
  "ATAC" = "#f768a1",  # Processing note.
  "H3K27ac" = "#DC0000FF",  # Processing note.
  "H3K4me3" = "#41ab5d",  # Processing note.
  "H3K27me3" = "#969696",  # Processing note.
  "H3K4me1" ="#fed976",  # Processing note.
  "RNASeq"="#41b6c4"
)

ggplot(df,
       aes(x = reorder(V2,-V3), y = (V3), fill = reorder(V1,V3),
           stratum = reorder(V1,V3),
           alluvium = reorder(V1,V3)
       )) +
  geom_alluvium(alpha=.35)+
  geom_stratum(width=0.6, size=0.1) +
  geom_bar(position='stack',stat='identity',width=0.6, color="white", linewidth=0.05) +
  theme_classic() +
  # coord_flip() + 
  theme(
    # legend.position = c(0.8,0.6),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    axis.text = element_text(size=14, color="black"),
    axis.title = element_text(size=14, color="black"),
    axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
  ) + 
  labs(y="Count", x="Tissue") +
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
  )

ggsave("AS-SNPinformation.png")
ggsave("AS-SNPinformation.pdf")


head(df)
df_prop = df %>%
  group_by(V1) %>%
  mutate(
    group_total = sum(V3),  # Processing note.
    percentage = V3 / group_total * 100,  # Processing note.
    proportion = V3 / group_total  # Processing note.
  ) %>%
  ungroup() %>%
  # select(V1, V2, V3, TisNum, group_total, percentage, proportion) %>%
  arrange(V1, V2, V3)

head(df_prop)
ggplot(df_prop,
       aes(x = reorder(V2,-percentage), y = (percentage), fill = reorder(V1,percentage),
           stratum = reorder(V1,percentage),
           alluvium = reorder(V1,percentage)
       )) +
  geom_alluvium(alpha=.35)+
  geom_stratum(width=0.6, size=0.1) +
  geom_bar(position='stack',stat='identity',width=0.6, color="white", linewidth=0.05) +
  theme_classic() +
  # coord_flip() + 
  theme(
    # legend.position = c(0.8,0.6),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    axis.text = element_text(size=14, color="black"),
    axis.title = element_text(size=14, color="black"),
    axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
  ) + 
  labs(y="Percentage (%)", x="Tissue")+
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
  )



ggsave("AS_SNP_counts_by_tissue_and_assay_pct.png")
ggsave("AS_SNP_counts_by_tissue_and_assay_pct.pdf")


# Processing note.
library(CMplot)
df = read.csv("251125_AS_Dis", sep = "\t", header = F)
head(df)

df$SNP = rownames(df)
head(df)

CMplot(df[,c(4,2,3)],plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="png",file.name=NULL,dpi=300,
       main="",file.output=T,verbose=TRUE,width=9,height=6)




##### REF-ALT freq #####
dir()
df = read.csv("251125_AS_Tis_REF-ALT", sep = "\t", header = F)

head(df)

# Processing note.
df_group <- df %>%
  group_by(across(c(V1,V3))) %>%
  summarise(sum = sum(V4))

head(df_group)

# Processing note.
df_group = df_group %>%
  group_by(V1) %>%
  mutate(
    group_total = sum(sum),  # Processing note.
    percentage = sum / group_total * 100,  # Processing note.
    proportion = sum / group_total  # Processing note.
  ) %>%
  ungroup() %>%
  select(V1, V3, sum, group_total, percentage, proportion) %>%
  arrange(V1, V3)

ggplot(df_group) +
  geom_col(aes(x = (percentage), y = reorder(V3,sum), fill = V1)) +
  labs(title = "",
       x = "REF-ALT",
       y = "Count",
       fill = "") +
  theme_classic()


library(ggalluvial)
head(df_group)
ggplot(df_group,
       aes(x = reorder(V3,percentage), y = percentage, fill = reorder(V1,percentage),
           stratum = reorder(V1,percentage),
           alluvium = reorder(V1,percentage)
       )) +
  geom_alluvium(alpha=.35)+
  geom_stratum(width=0.6, size=0.1) +
  geom_bar(position='stack',stat='identity',width=0.6, color="white", linewidth=0.05) +
  theme_classic() +
  coord_flip() + 
  theme(
    legend.position = c(0.7,0.3),
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    axis.text = element_text(size=14, color="black"),
    axis.title = element_text(size=14, color="black"),
    
  ) + 
  labs(y="Proportion (%)", x="Rel-Alt")+
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
  )

ggsave("250408_AS_SNP_mutation_types.png")
ggsave("250408_AS_SNP_mutation_types.pdf")
