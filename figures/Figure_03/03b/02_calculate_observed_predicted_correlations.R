#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
library(tidyverse)
library(data.table)
library(ggpointdensity)
library(viridis)
library(rstatix)

tissue_file <- "/data/home/sczd644/run/zsw_chrombpnet/tissue.txt"
tissues <- read_lines(tissue_file) %>% str_trim()

output_main_dir <- "/data/home/sczd644/run/zsw_chrombpnet/pred_bw/correlation_analysis"
dir.create(output_main_dir, recursive = TRUE, showWarnings = FALSE)
cat(" Output directory:", output_main_dir, "\n")

all_cor_results <- data.frame()

for (tissue in tissues) {
  cat(" :", tissue, "\n")

  real_file <- paste0("/data/home/sczd644/run/zsw_chrombpnet/ATAC_bams/", tissue, "/data/peaks_no_blacklist.observed.out")
  pred_file <- paste0("/data/home/sczd644/run/zsw_chrombpnet/pred_bw/", tissue, "/", tissue, "_peaks_no_blacklist.predicted.out")

  cat("   file :", real_file, "\n")
  cat("   file :", pred_file, "\n")

  if (!file.exists(real_file)) {
    cat("   :  file , ", tissue, "\n")
    next
  }
  if (!file.exists(pred_file)) {
    cat("   :  file , ", tissue, "\n")
    next
  }

  tryCatch({
    cat("   read ...\n")
    real <- fread(real_file)
    names(real) <- c("name", "size", "covered", "observed", "mean0", "mean")  # bigWigAverageOverBed output
    real <- real %>% select(name, observed)

    pred <- fread(pred_file)
    names(pred) <- c("name", "size", "covered", "predicted", "mean0", "mean")
    pred <- pred %>% select(name, predicted)

    cat("   read , :", nrow(real), " :", nrow(pred), "\n")

    df <- merge(real, pred, by = "name")
    cat("   :", nrow(df), "\n")

    # =========================
    # =========================

    if (nrow(df) < 10) {
      cat("   :  , ", tissue, "\n")
      next
    }

    get_density <- function(x, y, ...) {
      dens <- MASS::kde2d(x, y, ...)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
    }

    df$density <- get_density(df$observed, df$predicted, n = 100)

    # =========================
    # =========================
    eps <- 1e-6
    df$log_observed  <- log10(df$observed  + eps)
    df$log_predicted <- log10(df$predicted + eps)

    cor_df <- df %>% cor_test(log_observed, log_predicted) %>% mutate(Tissue = tissue)  # cor_test

    p <- ggplot(data = df, aes(x = log_observed, y = log_predicted)) +
      geom_pointdensity(show.legend = TRUE) +
      geom_smooth(se = FALSE, color = "red", method = "lm") +
      scale_color_viridis(option = "A", name = "Density") +
      xlab("Observed accessibility (log10)") +
      ylab("Predicted accessibility (log10)") +
      ggtitle(paste("Observed vs Predicted -", tissue)) +
      annotate("text",
               x = min(df$log_observed, na.rm = TRUE),
               y = max(df$log_predicted, na.rm = TRUE),
               hjust = 0, vjust = 1,
               label = paste0("r = ", round(cor_df$cor, 3),
                              "\np = ", format.pval(cor_df$p, digits = 3)),
               size = 5, color = "black") +
      theme_classic(base_size = 20, base_line_size = 0.5) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10)
      )

    output_dir <- paste0(output_main_dir, "/", tissue)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    ggsave(paste0(output_dir, "/peaks_no_blacklist.obs_vs_pred.pdf"), p, width = 4.5, height = 4)
    cat("   save\n")

    fwrite(cor_df, paste0(output_dir, "/peaks_no_blacklist.obs_vs_pred.txt"),
           sep = "\t", row.names = FALSE, quote = FALSE)
    cat("   results save, :", cor_df$cor, "\n")

    all_cor_results <- bind_rows(all_cor_results, cor_df)

  }, error = function(e) {
    cat("   ", tissue, " :", e$message, "\n")
  })

  cat("   ", tissue, " \n\n")
}

if (nrow(all_cor_results) > 0) {
  all_cor_results <- all_cor_results %>% arrange(cor)

  fwrite(all_cor_results, paste0(output_main_dir, "/all_tissues_correlation_results.txt"),
         sep = "\t", row.names = FALSE, quote = FALSE)
  cat(" results save :", paste0(output_main_dir, "/all_tissues_correlation_results.txt"), "\n")

  p_bar <- ggplot(all_cor_results, aes(x = cor, y = reorder(Tissue, cor))) +
    geom_col(fill = "lightblue", color = "black", width = 0.35) +
    geom_text(aes(label = round(cor, 3)),
              hjust = -0.1, size = 3) +
    labs(x = "Pearson Correlation Coefficient (log10 scale)",
         y = "Tissue",
         title = "Correlation between Observed and Predicted Accessibility (log10)") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

  ggsave(paste0(output_main_dir, "/tissue_correlation_barplot.pdf"), p_bar, width = 8, height = 12)
  cat(" save :", paste0(output_main_dir, "/tissue_correlation_barplot.pdf"), "\n")

  cat("\n statistics:\n")
  cat(" :", nrow(all_cor_results), "\n")
  cat(" :", round(mean(all_cor_results$cor, na.rm = TRUE), 3), "\n")
  cat(" :", round(min(all_cor_results$cor, na.rm = TRUE), 3), "-",
      round(max(all_cor_results$cor, na.rm = TRUE), 3), "\n")
} else {
  cat(" \n")
}
