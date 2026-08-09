# Processing note.
setwd('/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/Target_gene_cloest/')

# Processing note.
data <- read.table('AA_with_without_super_enhancer_compare_to_human.txt')
View(data)

# Processing note.
names(data) = c("ID","gene","chr","strand","start","end","total",
               "average","type","tau","humanID","humantochicken","conservation")

# Processing note.
library(ggpubr)
library(ggplot2) 

# Processing note.
data$type

# Processing note.
# Processing note.
pdf("04_AA_with_without_super_enhancer_compare_to_human_1.pdf", width=4.5, height=3)
ggplot(data, aes(x=log10(average), fill = type)) + 
  theme_classic() +
  geom_density(alpha = 0.3) +
  labs(x="Gene expression (log10)", y = "gene density")
dev.off()  # Processing note.

# Processing note.
# Processing note.
pdf("04_AA_with_without_super_enhancer_compare_to_human_2.pdf", width=4.5, height=3)
ggplot(data, aes(x=tau, fill = type)) + 
  theme_classic() +
  geom_density(alpha = 0.3) +
  labs(x="tau", y = "gene density")
dev.off()  # Processing note.

# Processing note.
# Processing note.
pdf("04_AA_with_without_super_enhancer_compare_to_human_3.pdf", width=4.5, height=3)
ggplot(data, aes(x=conservation, fill = type)) + 
  theme_classic() +
  geom_density(alpha = 0.3) +
  labs(x="conservation", y = "gene density")
dev.off()  # Processing note.