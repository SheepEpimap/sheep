# -*- coding: utf-8 -*-
options(encoding = "UTF-8")

# Install and load required packages
if (!require("VennDiagram")) {
  install.packages(
    "VennDiagram",
    dependencies = TRUE,
    repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
  )
  library(VennDiagram)
}
if (!require("grid")) {
  install.packages("grid", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(grid)
}

# =========================
# Processing note.
# =========================
input_dir <- "/vol2/wulingyun/AS/2_binomial_results/fdr_less_0.1_batch_results/venn_preprocess_5col"
output_dir <- "/vol2/wulingyun/AS/AS_VISUAL/venn"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created venn subdirectory: ", output_dir, "\n")
} else {
  cat("Venn subdirectory already exists: ", output_dir, "\n")
}

# Processing note.
input_files <- c(
  "39"       = file.path(input_dir, "39_gene_table.tsv"),
  "40"       = file.path(input_dir, "40_gene_table.tsv"),
  "combined" = file.path(input_dir, "combined_39_40_gene_table.tsv")
)

# Processing note.
plot_titles <- c(
  "39"       = "Individual 39",
  "40"       = "Individual 40",
  "combined" = "Combined 39 + 40"
)

# Processing note.
output_prefix <- c(
  "39"       = "5venn_39data",
  "40"       = "5venn_40data",
  "combined" = "5venn_combined_39_40"
)

# =========================
# Processing note.
# =========================

# Processing note.
open_pdf_device <- function(filename, width = 10, height = 10) {
  if (capabilities("cairo")) {
    cairo_pdf(filename = filename, width = width, height = height, onefile = TRUE)
  } else {
    pdf(file = filename, width = width, height = height, onefile = TRUE)
  }
}

# Processing note.
read_venn_sets <- function(file) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  venn_data <- read.delim(
    file = file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    na.strings = "",
    check.names = FALSE
  )

  required_cols <- c("ATAC", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3")
  missing_cols <- setdiff(required_cols, colnames(venn_data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in file: ", file, "\n",
      "Missing columns: ", paste(missing_cols, collapse = ", ")
    )
  }

  set_ATAC     <- unique(venn_data$ATAC[!is.na(venn_data$ATAC) & venn_data$ATAC != ""])
  set_H3K27ac  <- unique(venn_data$H3K27ac[!is.na(venn_data$H3K27ac) & venn_data$H3K27ac != ""])
  set_H3K27me3 <- unique(venn_data$H3K27me3[!is.na(venn_data$H3K27me3) & venn_data$H3K27me3 != ""])
  set_H3K4me1  <- unique(venn_data$H3K4me1[!is.na(venn_data$H3K4me1) & venn_data$H3K4me1 != ""])
  set_H3K4me3  <- unique(venn_data$H3K4me3[!is.na(venn_data$H3K4me3) & venn_data$H3K4me3 != ""])

  cat("\nNumber of genes in each set for:", basename(file), "\n")
  cat("ATAC:       ", length(set_ATAC), "\n")
  cat("H3K27ac:    ", length(set_H3K27ac), "\n")
  cat("H3K27me3:   ", length(set_H3K27me3), "\n")
  cat("H3K4me1:    ", length(set_H3K4me1), "\n")
  cat("H3K4me3:    ", length(set_H3K4me3), "\n")

  venn_sets <- list(
    "ATAC"     = set_ATAC,
    "H3K27ac"  = set_H3K27ac,
    "H3K27me3" = set_H3K27me3,
    "H3K4me1"  = set_H3K4me1,
    "H3K4me3"  = set_H3K4me3
  )

  return(venn_sets)
}

# Processing note.
make_venn_grob <- function(venn_sets) {
  venn_plot <- venn.diagram(
    x = venn_sets,
    category.names = c("ATAC", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3"),
    filename = NULL,
    imagetype = "png",

    # Processing note.
    height = 2400,
    width = 2400,
    resolution = 300,

    # Processing note.
    col = "white",
    lty = 1,
    lwd = 1,
    fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2", "#eadff0", "#fff2cd"),
    alpha = 0.9,

    # Processing note.
    label.col = "black",
    label.cex = 0.10,
    fontfamily = "serif",
    fontface = "bold",

    # Processing note.
    cat.col = c("#cb6274", "#7ba498", "#687d94", "#81668b", "#ffcf5c"),
    cat.cex = 1.5,
    cat.fontfamily = "serif",
    cat.fontface = "bold",
    cat.pos = c(0, -30, -130, 130, 40),
    cat.dist = 0.25
  )

  return(venn_plot)
}

# Processing note.
draw_venn_page <- function(venn_plot, title_text) {
  grid.newpage()

  # Processing note.
  vp_main <- viewport(
    x = 0.5,
    y = 0.6,
    width = 0.8,
    height = 0.8
  )
  pushViewport(vp_main)
  grid.draw(venn_plot)
  popViewport()

  # Processing note.
  vp_text <- viewport(
    x = 0.5,
    y = 0.1,
    width = 1,
    height = 0.2
  )
  pushViewport(vp_text)
  grid.text(
    label = title_text,
    x = 0.5,
    y = 0.5,
    gp = gpar(
      fontsize = 40,
      fontfamily = "serif",
      fontface = "bold",
      col = "black"
    )
  )
  popViewport()
}

# Processing note.
save_venn_plot <- function(venn_plot, title_text, out_base) {
  png_file <- paste0(out_base, ".png")
  pdf_file <- paste0(out_base, ".pdf")

  # PNG
  png(
    filename = png_file,
    height = 3000,
    width = 3000,
    res = 300,
    type = "cairo"
  )
  draw_venn_page(venn_plot, title_text)
  dev.off()

  # PDF
  open_pdf_device(
    filename = pdf_file,
    width = 10,
    height = 10
  )
  draw_venn_page(venn_plot, title_text)
  dev.off()

  cat("\nSaved files:\n")
  cat("PNG: ", png_file, "\n")
  cat("PDF: ", pdf_file, "\n")
}

# =========================
# Processing note.
# =========================
for (nm in names(input_files)) {
  cat("\n============================\n")
  cat("Processing dataset:", nm, "\n")
  cat("Input file:", input_files[[nm]], "\n")

  venn_sets <- read_venn_sets(input_files[[nm]])
  venn_plot <- make_venn_grob(venn_sets)

  out_base <- file.path(output_dir, output_prefix[[nm]])
  save_venn_plot(venn_plot, plot_titles[[nm]], out_base)
}

cat("\n============================\n")
cat("All done!\n")
cat("Output directory: ", output_dir, "\n")
cat("============================\n")