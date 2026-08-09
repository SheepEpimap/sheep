library(ggplot2)

setwd("/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/AA_one_count/")
data <- read.csv("all_super_enhancer_Gs_one_count.csv", sep = "\t", header = TRUE)
data1 <- as.data.frame(table(data$X0))

ymax <- max(data1$Freq)

pdf("22_super_enhancer_Gs_one_summary.pdf", width = 6.5, height = 4.23)

ggplot(data = data1, aes(x = Var1, y = Freq, fill = Freq)) +
  geom_col(position = position_dodge(1), width = 0.8) +
  ylab("Number of super enhancer") + xlab("Number of tissues") +
  scale_x_discrete(limits = levels(data1$Var1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, ymax * 1.10), clip = "off") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 12),
    axis.text.x = element_text(angle = 360, hjust = 0.5, colour = "black", family = "Times", size = 8),
    axis.text.y = element_text(family = "Times", size = 12, face = "plain"),
    axis.title.y = element_text(family = "Times", size = 15, face = "plain"),
    axis.title.x = element_text(family = "Times", size = 15, face = "plain")
  ) +
  geom_text(
    aes(label = Freq),
    position = position_dodge(width = 0.9),
    angle = 60, vjust = -1, hjust = -0.3, size = 2
  )

dev.off()