#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
library(ggpubr)
library(ggplot2)
library(rstatix)
library(dplyr)

setwd("/vol2/zhangshiwen/sheep_cor/h3k27ac/summary")

############################################
# inputfile
############################################
expr_file <- "H3K27ac_output_E5_gene_count_top005_expression_all_2.txt"
cons_file <- "H3K27ac_output_E5_gene_count_top005_expression_conservation_2.txt"

############################################
# read expression + tau file
############################################
data1 <- read.table(
  expr_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)

colnames(data1) <- c(
  "ID", "number", "gene", "chr", "strand", "start", "end",
  "total", "average", "type", "tau"
)[1:ncol(data1)]

data1$number  <- as.numeric(data1$number)
data1$start   <- as.numeric(data1$start)
data1$end     <- as.numeric(data1$end)
data1$total   <- as.numeric(data1$total)
data1$average <- as.numeric(data1$average)
data1$tau     <- as.numeric(data1$tau)

data1$type <- factor(data1$type, levels = c("Top", "Middle", "Bottom"))

############################################
# read conservation file
############################################
data2 <- read.table(
  cons_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)

colnames(data2) <- c(
  "ID", "number", "gene", "chr", "strand", "start", "end",
  "total", "average", "type", "tau", "human_id",
  "humantochicken", "conservation", "q"
)[1:ncol(data2)]

data2$number          <- as.numeric(data2$number)
data2$start           <- as.numeric(data2$start)
data2$end             <- as.numeric(data2$end)
data2$total           <- as.numeric(data2$total)
data2$average         <- as.numeric(data2$average)
data2$tau             <- as.numeric(data2$tau)
data2$humantochicken  <- as.numeric(data2$humantochicken)
data2$conservation    <- as.numeric(data2$conservation)
data2$q               <- as.numeric(data2$q)

data2$type <- factor(data2$type, levels = c("Top", "Middle", "Bottom"))

############################################
############################################
my_comparisons <- list(
  c("Top", "Middle"),
  c("Top", "Bottom"),
  c("Middle", "Bottom")
)

############################################
############################################
make_ypos <- function(ymax, step_min, step_max, n_comp) {
  if (n_comp <= 0) return(numeric(0))
  ymax + seq(step_min, step_max, length.out = n_comp)
}

############################################
# 1. average
############################################
data_clean_average <- do.call(
  rbind,
  lapply(split(data1, data1$type), function(df) {
    bx <- boxplot.stats(df$average)
    df[df$average >= bx$stats[1] & df$average <= bx$stats[5], ]
  })
)

cat("\n========== Average summary ==========\n")
print(
  data_clean_average %>%
    group_by(type) %>%
    summarise(
      n = n(),
      mean_average = mean(average, na.rm = TRUE),
      median_average = median(average, na.rm = TRUE),
      sd_average = sd(average, na.rm = TRUE),
      IQR_average = IQR(average, na.rm = TRUE)
    )
)

# ANOVA
anova_average <- anova_test(
  data = data_clean_average,
  dv = average,
  between = type
)
cat("\nANOVA results for average:\n")
print(anova_average)

pairwise_average <- data_clean_average %>%
  pairwise_t_test(average ~ type, p.adjust.method = "bonferroni")

cat("\nPairwise comparisons for average:\n")
print(pairwise_average)

ymax_avg <- max(data_clean_average$average, na.rm = TRUE)
pval_average <- pairwise_average %>%
  mutate(y.position = make_ypos(ymax_avg, 5, 15, n()))

p1 <- ggboxplot(
  data_clean_average,
  x = "type",
  y = "average",
  color = "type",
  outlier.shape = NA
) +
  labs(y = "Average gene expression", x = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, colour = "black", family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 20)
  ) +
  stat_pvalue_manual(pval_average, label = "p.adj.signif", tip.length = 0) +
  expand_limits(y = ymax_avg + 20)

############################################
# 2. tau
############################################
cat("\n========== Tau summary ==========\n")
print(
  data1 %>%
    group_by(type) %>%
    summarise(
      n = n(),
      mean_tau = mean(tau, na.rm = TRUE),
      median_tau = median(tau, na.rm = TRUE),
      sd_tau = sd(tau, na.rm = TRUE),
      IQR_tau = IQR(tau, na.rm = TRUE)
    )
)

# Kruskal-Wallis
kruskal_tau <- kruskal_test(data1, tau ~ type)
cat("\nKruskal results for tau:\n")
print(kruskal_tau)

pairwise_tau <- data1 %>%
  pairwise_wilcox_test(tau ~ type, p.adjust.method = "bonferroni")

cat("\nPairwise comparisons for tau:\n")
print(pairwise_tau)

tau_cor <- cor.test(data1$number, data1$tau, method = "spearman", exact = FALSE)
cat("\nSpearman correlation between enhancer number and tau:\n")
print(tau_cor)

ymax_tau <- max(data1$tau, na.rm = TRUE)
pval_tau <- pairwise_tau %>%
  mutate(y.position = make_ypos(ymax_tau, 0.03, 0.12, n()))

p2 <- ggboxplot(
  data1,
  x = "type",
  y = "tau",
  color = "type",
  outlier.shape = NA
) +
  labs(y = "Tau value", x = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, colour = "black", family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 20)
  ) +
  stat_pvalue_manual(pval_tau, label = "p.adj.signif", tip.length = 0) +
  expand_limits(y = ymax_tau + 0.15)

############################################
# 3. conservation
############################################
cat("\n========== Conservation summary ==========\n")
print(
  data2 %>%
    group_by(type) %>%
    summarise(
      n = n(),
      mean_conservation = mean(conservation, na.rm = TRUE),
      median_conservation = median(conservation, na.rm = TRUE),
      sd_conservation = sd(conservation, na.rm = TRUE),
      IQR_conservation = IQR(conservation, na.rm = TRUE)
    )
)

# Kruskal-Wallis
kruskal_conservation <- kruskal_test(data2, conservation ~ type)
cat("\nKruskal results for conservation:\n")
print(kruskal_conservation)

pairwise_conservation <- data2 %>%
  pairwise_wilcox_test(conservation ~ type, p.adjust.method = "bonferroni")

cat("\nPairwise comparisons for conservation:\n")
print(pairwise_conservation)

cons_cor <- cor.test(data2$number, data2$conservation, method = "spearman", exact = FALSE)
cat("\nSpearman correlation between enhancer number and conservation:\n")
print(cons_cor)

ymax_cons <- max(data2$conservation, na.rm = TRUE)
pval_conservation <- pairwise_conservation %>%
  mutate(y.position = make_ypos(ymax_cons, 3, 12, n()))

p3 <- ggboxplot(
  data2,
  x = "type",
  y = "conservation",
  color = "type",
  outlier.shape = NA
) +
  labs(y = "Identical to human gene (%)", x = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, colour = "black", family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 20)
  ) +
  stat_pvalue_manual(pval_conservation, label = "p.adj.signif", tip.length = 0) +
  expand_limits(y = ymax_cons + 15)

############################################
############################################
pdf("all_three_boxplots_new.pdf", width = 9, height = 4)
ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

############################################
############################################
pdf("tau_conservation_scatter.pdf", width = 8, height = 4)

p_tau_scatter <- ggplot(data1, aes(x = number, y = tau)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(x = "Enhancer number", y = "Tau")

p_cons_scatter <- ggplot(data2, aes(x = number, y = conservation)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(x = "Enhancer number", y = "Conservation")

ggarrange(p_tau_scatter, p_cons_scatter, ncol = 2, nrow = 1)
dev.off()
