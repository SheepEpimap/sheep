#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(grid)
})

# ================================
# Processing note.
# ================================
OUT_W <- 8
OUT_H <- 8
N_PER_TISSUE <- 1500
USE_RASTER_LINE <- FALSE
RASTER_DPI <- 300
RASTER_DEV <- "ragg_png"
HAS_GGRASTR <- requireNamespace("ggrastr", quietly = TRUE)

if (USE_RASTER_LINE && !HAS_GGRASTR) {
  stop("USE_RASTER_LINE=TRUE, but ggrastr is not installed: install.packages('ggrastr')")
}

# ================================
# Processing note.
# ================================
setwd(paste0(
  "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/",
  "Merge_chromatin_state/state_variability/AA_super_enhancer/",
  "AA_super_result_1/"
))

plotmap_file <- "all_tissues_super_enhancer.txt"
summary_file <- "all_super_enhancer_summary.txt"

if (!file.exists(plotmap_file)) stop("File not found: ", plotmap_file)
if (!file.exists(summary_file)) stop("File not found: ", summary_file)

bar_w_mm <- round(58 * 0.60)
shrink_right_mm <- 10

# ================================
# Processing note.
# ================================
x1 <- 2500
x2 <- 50000
k0 <- 1
k1 <- 0.15
k2 <- 1.0

compress_x <- function(x, x1 = 2500, x2 = 50000,
                       k0 = 0.12, k1 = 0.015, k2 = 1.0) {
  x <- as.numeric(x)
  y1 <- x1 * k0
  y2 <- y1 + (x2 - x1) * k1
  ifelse(
    x <= x1, x * k0,
    ifelse(x <= x2,
           y1 + (x - x1) * k1,
           y2 + (x - x2) * k2)
  )
}

# ================================
# Processing note.
# ================================
system_tissue_order <- list(
  Nervous_System = c(
    "cerebral-cortex", "midbrain", "cerebellum", "brainstem",
    "hippocampus", "hypothalamus", "medulla-oblongata",
    "optic-chiasm", "pineal", "pituitary", "pons", "splenium"
  ),
  Digestive_System = c(
    "rumen", "reticulum", "omasum", "abomasum", "duodenum",
    "jejunum", "ileum", "cecum", "colon", "rectum"
  ),
  Reproductive_System = c(
    "cervix", "cornua-uteri", "corpus-uteri", "ovary", "oviduct",
    "epididymis", "testis", "mammary-gland"
  ),
  Immune_System = c(
    "bone-marrow", "lymph-node", "thymus", "thyroid", "spleen"
  ),
  Endocrine_System = c("liver"),
  Renal_System = c("kidney"),
  Respiratory_System = c("lung"),
  Cardiovascular_System = c("heart"),
  Muscular_System = c("muscle"),
  Adipose = c("adipose"),
  Skin = c("skin", "soft-horn")
)

system_order <- names(system_tissue_order)
tissue_order <- unlist(system_tissue_order, use.names = FALSE)

# Processing note.
if (length(tissue_order) != 43L) {
  stop("The mapping table must contain 43 tissues; found: ", length(tissue_order))
}
if (length(system_order) != 11L) {
  stop("The mapping table must contain 11 systems; found: ", length(system_order))
}
if (anyDuplicated(tissue_order)) {
  stop("The mapping table contains duplicated tissues: ",
       paste(unique(tissue_order[duplicated(tissue_order)]), collapse = ", "))
}

tissue_cols <- c(
  "abomasum"="#f18264", "adipose"="#ffc4f1", "bone-marrow"="#d84b4b",
  "brainstem"="#f4d578", "cecum"="#daaa6c", "cerebellum"="#efd80b",
  "cerebral-cortex"="#dcd71a", "cervix"="#b6d7a9", "colon"="#f2c063",
  "cornua-uteri"="#69d683", "corpus-uteri"="#80d897", "duodenum"="#eb9d63",
  "epididymis"="#70d24b", "heart"="#bc58e3", "hippocampus"="#f1cb05",
  "hypothalamus"="#825e19", "ileum"="#ce9639", "jejunum"="#eb951c",
  "kidney"="#4f3136", "liver"="#ad8c8b", "lung"="#36b5f1",
  "lymph-node"="#c33a11", "mammary-gland"="#fed9d0",
  "medulla-oblongata"="#d0b35b", "midbrain"="#f9ed19",
  "muscle"="#a180ca", "omasum"="#f8ae81", "optic-chiasm"="#fece01",
  "ovary"="#69d28c", "oviduct"="#79ffaa", "pineal"="#807120",
  "pituitary"="#f1d95d", "pons"="#fcd222", "rectum"="#efe0a8",
  "reticulum"="#d97c68", "rumen"="#fc9891", "skin"="#d09dc5",
  "soft-horn"="#a25d73", "spleen"="#962932", "splenium"="#7c6919",
  "testis"="#7ef351", "thymus"="#ff3a32", "thyroid"="#f359d1"
)

missing_colors <- setdiff(tissue_order, names(tissue_cols))
if (length(missing_colors) > 0L) {
  stop("Colors are missing for these tissues: ", paste(missing_colors, collapse = ", "))
}
tissue_cols <- tissue_cols[tissue_order]

# ================================
# Processing note.
# ================================
system_shapes <- c(
  Nervous_System = 10,
  Digestive_System = 7,
  Reproductive_System = 12,
  Immune_System = 5,
  Endocrine_System = 0,
  Renal_System = 6,
  Respiratory_System = 2,
  Cardiovascular_System = 8,
  Muscular_System = 1,
  Adipose = 3,
  Skin = 4
)

if (!identical(names(system_shapes), system_order)) {
  stop("system_shapes and system_tissue_order have inconsistent names or order")
}

system_shape_map <- data.frame(
  system = system_order,
  shape = unname(system_shapes[system_order]),
  stringsAsFactors = FALSE
)

tissue2system <- data.frame(
  tissue1 = tissue_order,
  system = rep(system_order, lengths(system_tissue_order)),
  stringsAsFactors = FALSE
)
tissue2system$system <- factor(tissue2system$system, levels = system_order)

plotmap_raw <- read.table(
  plotmap_file,
  sep="\t",
  header=FALSE,
  stringsAsFactors=FALSE,
  check.names=FALSE,
  fill=TRUE,
  quote="",
  comment.char="",
  colClasses="character"
)

if (ncol(plotmap_raw) < 4L) {
  stop("plotmap must contain at least four columns; found: ", ncol(plotmap_raw))
}

colnames(plotmap_raw) <- paste0("V", seq_len(ncol(plotmap_raw)))
plotmap_raw[] <- lapply(plotmap_raw, function(x) trimws(as.character(x)))

number_test <- suppressWarnings(as.numeric(plotmap_raw[[1]]))
score_test  <- suppressWarnings(as.numeric(plotmap_raw[[3]]))
if (
  nrow(plotmap_raw) >= 2L &&
  (is.na(number_test[1]) || is.na(score_test[1])) &&
  sum(is.finite(number_test[-1]) & is.finite(score_test[-1])) >= 1L
) {
  message("Detected and removed the plotmap header row")
  plotmap_raw <- plotmap_raw[-1, , drop=FALSE]
  number_test <- suppressWarnings(as.numeric(plotmap_raw[[1]]))
  score_test  <- suppressWarnings(as.numeric(plotmap_raw[[3]]))
}

if (any(!is.finite(number_test))) {
  bad_rows <- which(!is.finite(number_test))
  stop(
    "Non-numeric values occur in plotmap column 1 (Number); example rows: ",
    paste(head(bad_rows, 10), collapse=", ")
  )
}
if (any(!is.finite(score_test))) {
  bad_rows <- which(!is.finite(score_test))
  stop(
    "Non-numeric values occur in plotmap column 3 (Score); example rows: ",
    paste(head(bad_rows, 10), collapse=", ")
  )
}

candidate_idx <- seq.int(4L, ncol(plotmap_raw))
candidate_stats <- data.frame(
  column = colnames(plotmap_raw)[candidate_idx],
  index = candidate_idx,
  matched_rows = vapply(
    candidate_idx,
    function(j) sum(plotmap_raw[[j]] %in% tissue_order, na.rm=TRUE),
    integer(1)
  ),
  matched_tissues = vapply(
    candidate_idx,
    function(j) length(unique(plotmap_raw[[j]][plotmap_raw[[j]] %in% tissue_order])),
    integer(1)
  ),
  stringsAsFactors=FALSE
)

message("Candidate plotmap column match statistics:")
for (i in seq_len(nrow(candidate_stats))) {
  message(sprintf(
    "  %s:matched rows=%d；matched tissues=%d",
    candidate_stats$column[i],
    candidate_stats$matched_rows[i],
    candidate_stats$matched_tissues[i]
  ))
}

best_rows <- candidate_stats[
  candidate_stats$matched_rows == max(candidate_stats$matched_rows),
  , drop=FALSE
]
best_rows <- best_rows[
  best_rows$matched_tissues == max(best_rows$matched_tissues),
  , drop=FALSE
]

if (nrow(best_rows) == 0L || best_rows$matched_rows[1] == 0L) {
  stop("No plotmap column from column 4 onward matches the 43-tissue mapping table")
}

if (nrow(best_rows) > 1L) {
  warning(
    "Multiple columns have the same tissue match count; using the first: ",
    paste(best_rows$column, collapse=", ")
  )
}

tissue_col_idx <- best_rows$index[1]
tissue_col_name <- best_rows$column[1]
message("Automatically detected tissue column: ", tissue_col_name, " (original column ", tissue_col_idx, ")")

# Processing note.
for (j in candidate_idx) {
  vals <- unique(plotmap_raw[[j]])
  vals <- vals[!is.na(vals) & nzchar(vals)]
  message(
    "  ", colnames(plotmap_raw)[j], "  examples: ",
    paste(head(vals, 6), collapse=" | ")
  )
}

plotmap_full <- data.frame(
  Number  = number_test,
  ID      = plotmap_raw[[2]],
  Score   = score_test,
  tissue1 = plotmap_raw[[tissue_col_idx]],
  stringsAsFactors=FALSE,
  check.names=FALSE
)

plotmap_full$tissue <- plotmap_full$tissue1

plotmap_tissues <- unique(plotmap_full$tissue1)
plotmap_tissues <- plotmap_tissues[!is.na(plotmap_tissues) & nzchar(plotmap_tissues)]
unmatched_plotmap <- setdiff(plotmap_tissues, tissue_order)
missing_plotmap <- setdiff(tissue_order, plotmap_tissues)

if (length(unmatched_plotmap) > 0L) {
  stop(
    "Values in the detected tissue column are absent from the mapping table: ",
    paste(unmatched_plotmap, collapse=", ")
  )
}
if (length(missing_plotmap) > 0L) {
  warning(
    "Tissues in the mapping table are absent from plotmap: ",
    paste(missing_plotmap, collapse=", ")
  )
}

plotmap_full$tissue1 <- factor(plotmap_full$tissue1, levels=tissue_order)
plotmap_full$tissue  <- factor(plotmap_full$tissue,  levels=tissue_order)
plotmap_full$x_plot <- compress_x(
  plotmap_full$Number, x1=x1, x2=x2, k0=k0, k1=k1, k2=k2
)

# Processing note.
plotmap_draw <- plotmap_full %>%
  group_by(tissue1) %>%
  arrange(Number, .by_group=TRUE) %>%
  slice(round(seq(1, n(), length.out=min(N_PER_TISSUE, n())))) %>%
  ungroup()

# Processing note.
peak_data <- plotmap_full %>%
  group_by(tissue1) %>%
  slice_max(Score, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  left_join(tissue2system, by="tissue1") %>%
  mutate(system=factor(system, levels=system_order))

if (anyNA(peak_data$system)) {
  stop("Some peak tissues could not be assigned to a system")
}

# ================================
# Processing note.
# ================================
x_breaks_raw <- c(
  0, 2500, 50000, 60000, 70000, 80000,
  90000, 100000, 110000, 120000, 130000
)
x_breaks_pos <- compress_x(
  x_breaks_raw, x1=x1, x2=x2, k0=k0, k1=k1, k2=k2
)

line_layer <- geom_line(linewidth=0.6)
if (USE_RASTER_LINE) {
  line_layer <- ggrastr::rasterise(
    line_layer, dpi=RASTER_DPI, dev=RASTER_DEV, scale=1
  )
}

p_curve <- ggplot(
  plotmap_draw,
  aes(x=x_plot, y=Score, color=tissue1, group=tissue1)
) +
  line_layer +
  geom_point(
    data=peak_data,
    aes(x=x_plot, y=Score, shape=system),
    inherit.aes=FALSE,
    color="black", fill=NA, size=3
  ) +
  scale_shape_manual(
    values=system_shapes,
    breaks=system_order,
    labels=system_order,
    name="System",
    drop=FALSE
  ) +
  scale_color_manual(
    values=tissue_cols,
    breaks=tissue_order,
    guide="none",
    drop=FALSE
  ) +
  scale_x_continuous(
    breaks=x_breaks_pos,
    labels=as.character(x_breaks_raw),
    limits=range(compress_x(
      c(0, 130000), x1=x1, x2=x2, k0=k0, k1=k1, k2=k2
    )),
    expand=c(0, 0)
  ) +
  theme_minimal() +
  theme(
    panel.grid=element_blank(),
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    axis.text.x=element_text(angle=30, hjust=1),
    legend.position="bottom",
    legend.direction="horizontal",
    plot.margin=margin(t=6, r=6+shrink_right_mm, b=6, l=6, unit="mm"),
    axis.ticks.x=element_line(color="black", linewidth=0.35),
    axis.ticks.y=element_line(color="black", linewidth=0.35),
    axis.ticks.length=unit(1.2, "mm")
  ) +
  guides(shape=guide_legend(nrow=3, byrow=TRUE)) +
  xlab("Enhancer rank by H3K27ac signal") +
  ylab("Score")

# ================================
# Processing note.
# ================================
sum_df <- read.table(
  summary_file, sep="\t", header=FALSE,
  stringsAsFactors=FALSE, fill=TRUE, check.names=FALSE
)
base_names <- c("Tissue", "Number", "Size", "H3K27ac_signal")
num_extra <- ncol(sum_df) - length(base_names)
extra_names <- if (num_extra > 0) paste0("Extra", seq_len(num_extra)) else character(0)
colnames(sum_df) <- c(base_names, extra_names)
sum_df$Tissue <- trimws(as.character(sum_df$Tissue))

summary_tissues <- unique(sum_df$Tissue)
summary_tissues <- summary_tissues[!is.na(summary_tissues) & nzchar(summary_tissues)]
unmatched_summary <- setdiff(summary_tissues, tissue_order)
missing_summary <- setdiff(tissue_order, summary_tissues)

if (length(unmatched_summary) > 0L) {
  stop("Summary contains tissues absent from the mapping table: ",
       paste(unmatched_summary, collapse = ", "))
}
if (length(missing_summary) > 0L) {
  warning("Tissues in the mapping table are absent from summary: ",
          paste(missing_summary, collapse = ", "))
}

sum_df <- sum_df %>%
  select(Tissue, Number) %>%
  filter(Tissue %in% tissue_order) %>%
  mutate(
    Number=as.numeric(Number),
    Tissue=factor(Tissue, levels=rev(tissue_order))
  )

if (anyNA(sum_df$Number)) {
  stop("The summary Number column contains non-numeric or missing values")
}

# ================================
# Processing note.
# ================================
flat_win <- plotmap_full %>% filter(Number >= 30000, Number <= 35000)
if (nrow(flat_win) < 50) {
  flat_win <- plotmap_full %>% filter(Number <= 2500)
}
if (nrow(flat_win) == 0) stop("No data are available for baseline calculation")

baseline_q <- 0.05
baseline_y <- as.numeric(quantile(flat_win$Score, probs=baseline_q, na.rm=TRUE))

b <- ggplot_build(p_curve)
pp <- b$layout$panel_params[[1]]
yr <- if (!is.null(pp$y.range)) {
  pp$y.range
} else if (!is.null(pp$y$range)) {
  pp$y$range
} else {
  range(plotmap_full$Score, na.rm=TRUE)
}

ymin <- yr[1]
ymax <- yr[2]
if (!is.finite(ymin) || !is.finite(ymax) || ymax <= ymin) {
  stop("The main-panel y-axis range is invalid")
}
bottom_frac <- (baseline_y - ymin) / (ymax - ymin)
bottom_frac <- max(0, min(0.95, bottom_frac))

# ================================
# Processing note.
# ================================
p_bar_inset <- ggplot(sum_df, aes(x=Number, y=Tissue, fill=Tissue)) +
  geom_col(width=0.72, alpha=0.98) +
  scale_fill_manual(values=tissue_cols, guide="none", drop=FALSE) +
  scale_x_continuous(
    position="top",
    breaks=pretty_breaks(n=4),
    labels=label_number(big.mark=","),
    expand=expansion(mult=c(0, 0.03))
  ) +
  labs(x="Super enhancer number", y=NULL) +
  theme_minimal(base_size=7) +
  theme(
    panel.grid=element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.line.x=element_line(color="black", linewidth=0.35),
    axis.ticks.x=element_line(color="black", linewidth=0.35),
    axis.ticks.length=unit(1.2, "mm"),
    panel.border=element_blank(),
    plot.background=element_rect(fill=NA, colour=NA),
    panel.background=element_rect(fill=NA, colour=NA),
    plot.margin=margin(0, 0, 0, 0)
  )

# ================================
# Processing note.
# ================================
p_main <- p_curve +
  inset_element(
    p_bar_inset,
    left=unit(0.5, "mm"),
    right=unit(bar_w_mm, "mm"),
    bottom=unit(bottom_frac, "npc"),
    top=unit(0.995, "npc"),
    align_to="panel"
  )

# ================================
# Processing note.
# ================================
tag <- if (USE_RASTER_LINE) {
  sprintf("sample%d_rast%d", N_PER_TISSUE, RASTER_DPI)
} else {
  sprintf("sample%d_vector", N_PER_TISSUE)
}

out_pdf <- paste0(
  "super_enhancer_plotmap_insetbar_piecewiseX_43tissues_",
  tag, ".pdf"
)
pdf(out_pdf, width=OUT_W, height=OUT_H, useDingbats=FALSE, compress=TRUE)
print(p_main)
dev.off()
message("Done! Output PDF: ", out_pdf)

out_jpg <- paste0(
  "super_enhancer_plotmap_insetbar_piecewiseX_43tissues_",
  tag, "_preview.jpg"
)
if (requireNamespace("ragg", quietly=TRUE)) {
  ragg::agg_jpeg(
    filename=out_jpg, width=OUT_W, height=OUT_H,
    units="in", res=200, quality=85
  )
  print(p_main)
  dev.off()
} else {
  ggsave(
    filename=out_jpg, plot=p_main,
    width=OUT_W, height=OUT_H, units="in",
    dpi=200, device="jpeg", bg="white", quality=85
  )
}
message("Done! Output JPG: ", out_jpg)
message("System order: ", paste(system_order, collapse=" | "))
message("Tissue number: ", length(tissue_order))
