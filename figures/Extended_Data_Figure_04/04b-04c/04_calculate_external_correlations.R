#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.

library(tidyverse)
library(data.table)
library(ggpointdensity)
library(viridis)
library(rstatix)

# ------------------------------------------------------------------
# ------------------------------------------------------------------
tissue_file <- "/data/home/sczd644/run/zsw_chrombpnet/public/tissue.txt"
tissues     <- read_lines(tissue_file) %>% str_trim()

# ------------------------------------------------------------------
# ------------------------------------------------------------------
output_main_dir <- "/data/home/sczd644/run/zsw_chrombpnet/public/pred_bw/correlation_analysis"
dir.create(output_main_dir, recursive = TRUE, showWarnings = FALSE)
cat(" Output directory:", output_main_dir, "\n")

# ------------------------------------------------------------------
# ------------------------------------------------------------------
all_cor_results <- data.frame()

# ------------------------------------------------------------------
# ------------------------------------------------------------------
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix   <- findInterval(x, dens$x)
  iy   <- findInterval(y, dens$y)
  dens$z[cbind(ix, iy)]
}

# ------------------------------------------------------------------
# ------------------------------------------------------------------
for (tissue in tissues) {
  cat(" :", tissue, "\n")

  tissue_dir <- paste0("/data/home/sczd644/run/zsw_chrombpnet/public/ATAC_bams/", tissue)

  if (!dir.exists(tissue_dir)) {
    tissue2 <- tolower(tissue)
    tissue_dir2 <- paste0("/data/home/sczd644/run/zsw_chrombpnet/public/ATAC_bams/", tissue2)
    if (dir.exists(tissue_dir2)) {
      tissue <- tissue2
      tissue_dir <- tissue_dir2
    } else {
      cat("   :  directory , ", tissue, "\n")
      next
    }
  }

  pair_dirs <- list.dirs(tissue_dir, recursive = FALSE, full.names = FALSE)
  pair_dirs <- pair_dirs[grepl("^Rep[0-9]{2}-[0-9]{2}$", pair_dirs)]

  if (length(pair_dirs) == 0) {
    cat("   :   RepXX-YY  directory, ", tissue, "\n")
    next
  }

  for (pair_tag in pair_dirs) {
    cat("   :", pair_tag, "\n")

    real_file <- paste0(tissue_dir, "/", pair_tag,
                        "/data/peaks_no_blacklist.observed.out")

    pred_dir  <- paste0("/data/home/sczd644/run/zsw_chrombpnet/public/pred_bw/", tissue, "/", pair_tag)
    pred_file <- paste0(pred_dir, "/", tissue, "_", pair_tag, "_peaks_no_blacklist.predicted.out")

    if (!file.exists(real_file)) {
      cat("     :  file , ", tissue, pair_tag, "\n")
      next
    }

    if (!file.exists(pred_file)) {
      if (dir.exists(pred_dir)) {
        hits <- list.files(pred_dir, pattern = "predicted\\.out$", full.names = TRUE)
        if (length(hits) > 0) pred_file <- hits[1]
      }
    }

    if (!file.exists(pred_file)) {
      cat("     :  file , ", tissue, pair_tag, "\n")
      next
    }

    tryCatch({
      real <- fread(real_file) %>% select(1, 4) %>% set_names(c("name", "observed"))
      pred <- fread(pred_file) %>% select(1, 4) %>% set_names(c("name", "predicted"))

      df <- merge(real, pred, by = "name")

      # df <- df[grepl("chr1_|chr3_|chr6_", df$name), ]

      if (nrow(df) < 10) {
        cat("     :  , ", tissue, pair_tag, "\n")
        next
      }

      df$density <- get_density(df$observed, df$predicted, n = 100)

      # ------------------------------------------------------------
      # ------------------------------------------------------------
      eps <- 1e-6
      df <- df %>%
        mutate(
          log_obs  = log10(observed + eps),
          log_pred = log10(predicted + eps)
        )

      cor_df <- cor_test(df, log_obs, log_pred) %>%
        mutate(Tissue = tissue, Pair = pair_tag)

      p <- ggplot(df, aes(log_obs, log_pred)) +
        geom_pointdensity(show.legend = TRUE) +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        scale_color_viridis(option = "A", name = "Density") +
        labs(x = "Observed accessibility (log10)",
             y = "Predicted accessibility (log10)",
             title = paste("Observed vs Predicted -", tissue, pair_tag)) +
        annotate("text",
                 x = min(df$log_obs, na.rm = TRUE),
                 y = max(df$log_pred, na.rm = TRUE),
                 hjust = 0, vjust = 1,
                 label = paste0("r = ", round(cor_df$cor, 3),
                                "\np = ", format.pval(cor_df$p, digits = 3)),
                 size = 5, color = "black") +
        theme_classic(base_size = 20, base_line_size = 0.5) +
        theme(legend.position = "none",
              axis.text.x = element_text(color = "black", size = 15),
              axis.text.y = element_text(color = "black", size = 15),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title      = element_text(size = 10))

      out_dir <- file.path(output_main_dir, tissue, pair_tag)
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      ggsave(file.path(out_dir, "peaks_no_blacklist.obs_vs_pred.pdf"),
             p, width = 4.5, height = 4)
      fwrite(cor_df,
             file.path(out_dir, "peaks_no_blacklist.obs_vs_pred.txt"),
             sep = "\t", row.names = FALSE, quote = FALSE)

      all_cor_results <- bind_rows(all_cor_results, cor_df)

    }, error = function(e) {
      cat("     ", tissue, pair_tag, " :", e$message, "\n")
    })

    cat("     ", tissue, pair_tag, "\n\n")
  }
}

# ------------------------------------------------------------------
# ------------------------------------------------------------------
if (nrow(all_cor_results) > 0) {

  all_cor_results <- all_cor_results %>%
    mutate(TissuePair = paste(Tissue, Pair, sep = "_")) %>%
    arrange(cor)

  fwrite(all_cor_results,
         file.path(output_main_dir, "all_tissues_correlation_results.txt"),
         sep = "\t", row.names = FALSE, quote = FALSE)

  p_bar <- ggplot(all_cor_results,
                  aes(x = cor, y = reorder(TissuePair, cor))) +
    geom_col(fill = "lightblue", colour = "black", width = 1) +
    geom_text(aes(label = sprintf("%.3f", cor)),
              hjust = -0.1, size = 3) +
    labs(x = "Pearson Correlation Coefficient (log10 scale)",
         y = "Tissue_Pair",
         title = "Correlation between Observed and Predicted Accessibility (log10)") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 8),
      plot.title       = element_text(hjust = 0.5),
      panel.spacing    = element_blank(),
      plot.margin      = margin(0, 0, 0, 0, "cm")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_y_discrete(expand = c(0, 0))

  bar_height_cm <- 1 * nrow(all_cor_results)
  ggsave(file.path(output_main_dir, "tissue_correlation_barplot.pdf"),
         p_bar,
         width = 8,
         height = bar_height_cm / 2.54)

  cat("\n statistics:\n")
  cat(" (tissue×pair):", nrow(all_cor_results), "\n")
  cat(" :", round(mean(all_cor_results$cor, na.rm = TRUE), 3), "\n")
  cat(" :",
      round(min(all_cor_results$cor, na.rm = TRUE), 3), "-",
      round(max(all_cor_results$cor, na.rm = TRUE), 3), "\n")
} else {
  cat(" \n")
}
