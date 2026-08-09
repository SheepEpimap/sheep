#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 3_GO_enrichemnt.R <out_dir> <gene_file>")
}

out_dir   <- args[1]
gene_file <- args[2]

# gene_file 
gene_path <- gene_file
if (!file.exists(gene_path)) {
  gene_path <- file.path(out_dir, gene_file)
}
if (!file.exists(gene_path)) {
  stop("Gene file not found: ", gene_file, " (also tried: ", gene_path, ")")
}

prefix <- sub("\\.txt$", "", basename(gene_path))
prefix <- sub("_gene_human$", "", prefix)
prefix <- sub("^human/", "", prefix)

genes <- readLines(gene_path, warn = FALSE)
genes <- trimws(genes)
genes <- genes[genes != "" & genes != "NA"]
genes <- unique(genes)

if (length(genes) == 0) {
  message("No genes found in: ", gene_path)
  quit(status = 0)
}

is_ensembl <- all(grepl("^ENSG", genes))

genes <- sub("\\..*$", "", genes)

from_keytype <- if (is_ensembl) "ENSEMBL" else "SYMBOL"

map_df <- suppressWarnings(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = genes,
    keytype  = from_keytype,
    columns  = c("ENTREZID", "SYMBOL")
  )
)

map_df <- map_df[!is.na(map_df$ENTREZID) & map_df$ENTREZID != "", , drop = FALSE]
entrez <- unique(map_df$ENTREZID)

if (length(entrez) == 0) {
  message("No ENTREZID mapped from input IDs. keytype=", from_keytype,
          "  file=", gene_path)
  quit(status = 0)
}

# Processing note.
ego <- enrichGO(
  gene          = entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

out_txt <- file.path(out_dir, paste0("GO_ENSEMBL_", prefix, "_BP.txt"))
out_pdf <- file.path(out_dir, paste0("GO_ENSEMBL_", prefix, "_BP.pdf"))

write.table(as.data.frame(ego), out_txt, sep = "\t", quote = FALSE, row.names = FALSE)

pdf(out_pdf, width = 9, height = 7)
df_ego <- as.data.frame(ego)
if (nrow(df_ego) == 0) {
  grid.newpage()
  grid.text(paste0(prefix, " (BP): no significant terms"), x = 0.5, y = 0.5)
} else {
  print(dotplot(ego, showCategory = min(20, nrow(df_ego))))
}
dev.off()

message("Done: ", prefix,
        " | input=", basename(gene_path),
        " | mapped_ENTREZ=", length(entrez),
        " | result_terms=", nrow(as.data.frame(ego)))
