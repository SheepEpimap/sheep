#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
setwd('/vol2/zhangshiwen/sheep_cor/h3k27ac/')
data <- read.table('H3K27ac_output_E5_confident.tsv',header = F)


names(data)=c("enhancer","cor","pvalue","id","Distance","qvalue")
data$logdistance=log10(data$Distance+1)
data$logqvalue=-log10(data$qvalue)
mean(data$Distance)
median(data$Distance)


cor.test(data$Distance,data$qvalue)
cor.test(data$logdistance,data$logqvalue)

mean(data$Distance)
max(data$Distance)

library(ggplot2)
library(aplot)
p1<-ggplot(data,aes(log10(Distance),-log10(qvalue)))+
  geom_hex()+
  scale_fill_gradientn(colours = c("#FAE9BA","#EA5C60", "#9C307F","#5C157C","black"))+
  stat_smooth(method = "lm",color = "#2bb1ff", fill = "lightgray")+
  theme_classic(base_size = 15)+
  labs(y = "-log10(qvalue)")+
  xlab("Distance to gene TSS(log10)")+
  annotate("text", x=1, y=10, label= paste("r = ", round(cor(data$logdistance,data$logqvalue),digits = 2),sep = ""), colour = "red",size=4)+
  annotate("text", x=2, y=9, label="p-value < 2.2e-16", colour = "red",size=4)



p2<-ggplot(data,aes(logdistance))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


p3<-ggplot(data,aes(logqvalue))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  coord_flip()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("output_E5_confident_p_distance.pdf", width=5.38, height=4.23)

p1%>%
  insert_top(p2,height = 0.5)%>%
  insert_right(p3,0.5)

dev.off()


p1=ggplot(data,aes(cor,-log10(qvalue)))+
  geom_hex()+
  scale_fill_gradientn(colours = c("#FAE9BA","#EA5C60", "#9C307F","#5C157C","black"))+
  theme_classic(base_size = 15)+
  labs(y = "-log10(qvalue)")+
  xlab("cor")+
  stat_smooth()


p2<-ggplot(data,aes(cor))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


p3<-ggplot(data,aes(logqvalue))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  coord_flip()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))


pdf("output_E5_confident_cor_qvalue.pdf", width=5.38, height=4.23)

p1%>%
  insert_top(p2,height = 0.5)%>%
  insert_right(p3,0.5)

dev.off()
