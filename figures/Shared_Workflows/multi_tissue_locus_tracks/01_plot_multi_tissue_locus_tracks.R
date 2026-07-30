#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(glue)
  library(ggplot2)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(Rsamtools)
  library(Biostrings)
  library(ggseqlogo)
  library(patchwork)
  library(grid)   #   mm
})

# ============================================================
# ============================================================
opt_list <- list(
  make_option(c("--tissues"), type="character", default="",
              help="[ ] /  id: , file (  id).  6  (obs,pred,nuc,corr,fp,logo) ."),
  make_option(c("--hits_tissues"), type="character", default="",
              help="Hits   id  : file .  --tissues  ."),

  # region mode (mutually exclusive)
  make_option(c("--only_region"), type="character", default="",
              help=" :chr:start-end.  --region_big/--region_small  ."),
  make_option(c("--region_big"), type="character", default="",
              help=" :  chr:start-end."),
  make_option(c("--region_small"), type="character", default="",
              help=" :  chr:start-end; ."),

  # track selection (6 basic metrics only; hits always drawn)
  make_option(c("--tracks"), type="character", default="obs,pred,nuc,corr,fp,logo",
              help=" ( ). :obs,pred,nuc,corr,fp,logo."),
  make_option(c("--big_tracks"), type="character", default="obs,pred,nuc,corr,fp",
              help=" ( ) ( ). :obs,pred,nuc,corr,fp,logo."),
  make_option(c("--small_tracks"), type="character", default="logo",
              help=" ( ) ( ). :obs,pred,nuc,corr,fp,logo."),

  make_option(c("--highlight"), type="logical", default=TRUE,
              help=" : ."),
  make_option(c("--highlight_alpha"), type="double", default=0.25,
              help=" (0~1, )."),

  make_option(c("-o","--outpdf"), type="character", default="track_multi.pdf",
              help="output PDF  ."),
  make_option(c("--base_pt"), type="double", default=6,
              help="  pt(  5~7)."),

  make_option(c("--tile_big"),  type="integer", default=25,
              help=" (obs/pred/nuc/fp) (bp)."),
  make_option(c("--tile_corr"), type="integer", default=5,
              help="corr(signed) (bp)."),

  make_option(c("--bin_stat"), type="character", default="mean",
              help=" statistics :mean|max|sum."),
  make_option(c("--pos_style"), type="character", default="area",
              help=" :area(  ribbon+line)  bars( )."),

  make_option(c("--hit_arrow_len"), type="integer", default=20,
              help="Hits   V   chevron  (  bp)."),
  make_option(c("--hit_arrow_step"), type="integer", default=60,
              help="Hits   bp   V   chevron(  bp)."),

  make_option(c("--out_w"), type="double", default=12,
              help="PDF  ( )."),
  make_option(c("--out_h"), type="double", default=8,
              help="PDF  ( );  --auto_height=TRUE  ."),

  make_option(c("--yq"), type="double", default=0.999,
              help=" (  *_q  )."),

  make_option(c("--yscale"), type="character", default="track",
              help=paste0(
                "y ( UCSC autoScale):",
                "track= / ;",
                "metric= y (max/min );",
                "track_q= (--yq) ;",
                "metric_q= (--yq) ."
              )),

  make_option(c("--auto_height"), type="logical", default=TRUE,
              help="  PDF  .TRUE:PDF  = ;FALSE:  --out_h."),
  make_option(c("--track_h"), type="double", default=0.95,
              help=" “ ”(obs/pred/nuc/corr/fp/logo)  ( )."),
  make_option(c("--track_h_scale"), type="double", default=0.55,
              help="Scale bar  ( )."),
  make_option(c("--track_h_genes"), type="double", default=0.95,
              help="Genes  ( )."),
  make_option(c("--track_h_hits"), type="double", default=0.95,
              help="Hits( )  ( )."),
  make_option(c("--track_h_map"), type="character", default="",
              help=" ( ),  key=value.key  :obs,pred,nuc,corr,fp,logo,scale,genes,hits. :'logo=1.6,hits=1.3'."),

  make_option(c("--logo_pad"), type="double", default=1.25,
              help="Logo   y  (>1  , , )."),

  make_option(c("--gene_arrow_step_bp"), type="integer", default=250,
              help="Genes  (bp), ."),
  make_option(c("--gene_arrow_seg_bp"), type="integer", default=30,
              help="Genes  (bp,  arrow  )."),
  make_option(c("--gene_arrow_mm"), type="double", default=1.6,
              help="Genes  ( , , ).")
)

opt <- parse_args(OptionParser(option_list = opt_list))
opt$bin_stat  <- tolower(opt$bin_stat)
opt$pos_style <- tolower(opt$pos_style)
opt$yscale    <- tolower(str_trim(opt$yscale))

if (!opt$bin_stat %in% c("mean","max","sum")) stop("--bin_stat  :mean|max|sum")
if (!opt$pos_style %in% c("area","bars")) stop("--pos_style  :area|bars")
if (!is.finite(opt$logo_pad) || opt$logo_pad <= 0) stop("--logo_pad  ")
if (!opt$yscale %in% c("track","metric","track_q","metric_q")) {
  stop("--yscale  :track|metric|track_q|metric_q")
}

allowed_metrics <- c("obs","pred","nuc","corr","fp","logo")

parse_id_list <- function(x) {
  x <- str_trim(x)
  if (x == "") return(character(0))
  if (file.exists(x)) {
    ids <- readLines(x, warn = FALSE)
    ids <- str_trim(ids)
    ids <- ids[ids != ""]
    return(ids[!duplicated(ids)])
  }
  ids <- unlist(strsplit(x, ",", fixed = TRUE))
  ids <- str_trim(ids)
  ids <- ids[ids != ""]
  ids[!duplicated(ids)]
}

parse_metric_list <- function(x) {
  x <- str_trim(x)
  if (x == "") return(character(0))
  m <- unlist(strsplit(x, ",", fixed = TRUE))
  m <- tolower(str_trim(m))
  m <- m[m != ""]
  m <- m[!duplicated(m)]
  bad <- setdiff(m, allowed_metrics)
  if (length(bad) > 0) stop(" :", paste(bad, collapse = ", "),
                            "\n :", paste(allowed_metrics, collapse = ","))
  m
}

tissues <- parse_id_list(opt$tissues)
if (length(tissues) == 0) stop("--tissues  ( file )")

hits_tissues <- if (str_trim(opt$hits_tissues) == "") tissues else parse_id_list(opt$hits_tissues)
if (length(hits_tissues) == 0) stop("--hits_tissues  , check")

allowed_height_keys <- c(allowed_metrics, "scale", "genes", "hits")

parse_height_map <- function(x) {
  x <- str_trim(x)
  if (x == "") return(setNames(numeric(0), character(0)))
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- numeric(0)
  for (p in parts) {
    p <- str_trim(p)
    if (p == "") next
    kv <- unlist(strsplit(p, "=", fixed = TRUE))
    if (length(kv) != 2) stop("--track_h_map token  :", p, "(  key=value)")
    k <- tolower(str_trim(kv[1]))
    v <- suppressWarnings(as.numeric(str_trim(kv[2])))
    if (!k %in% allowed_height_keys) {
      stop("  key:", k, "\n :", paste(allowed_height_keys, collapse = ","))
    }
    if (!is.finite(v) || v <= 0) stop(" ( ):", p)
    out[k] <- v
  }
  out
}

hmap <- parse_height_map(opt$track_h_map)
H <- function(key, default_h) {
  key <- tolower(key)
  if (key %in% names(hmap)) return(unname(hmap[[key]]))
  default_h
}

# ============================================================
# ============================================================
tissue_color_file <- "/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv"

gtf_gz  <- "/data/home/sczd644/run/zsw_chrombpnet/track/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz"
fa_file <- "/data/home/sczd644/run/zsw_chrombpnet/track/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna"

hits_bed_path <- function(id) glue("/data/home/sczd644/run/zsw_chrombpnet/track/hits/{id}_hits_tf.bed")
obs_bw_path   <- function(id) glue("/data/home/sczd644/run/zsw_chrombpnet/track/atac/{id}_obs_merged.bw")
pred_bw_path  <- function(id) glue("/data/home/sczd644/run/zsw_chrombpnet/track/atac/{id}_peaks_chrombpnet_nobias.bw")

nuc_bgz_path    <- function(id) glue("/data/home/sczd644/run/zsw_chrombpnet/region_ann/nucleoatac/{id}/{id}.occ.bedgraph.gz")
contrib_bw_path <- function(id) glue("/data/home/sczd644/run/zsw_chrombpnet/chrombpnet_contribs/{id}_chrombpnet_contribs/{id}_chrombpnet_contribs.counts_scores.bw")

tobias_corr_bw_path <- function(id) glue("/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/tobias/{id}_corrected.bw")
tobias_fp_bw_path   <- function(id) glue("/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/tobias/{id}_footprint.bw")

resolve_hits_bed <- function(id) {
  p <- hits_bed_path(id)
  if (!file.exists(p) && file.exists(paste0(p, ".gz"))) p <- paste0(p, ".gz")
  p
}

need0 <- c(tissue_color_file, gtf_gz, fa_file, paste0(fa_file, ".fai"))
miss0 <- need0[!file.exists(need0)]
if (length(miss0) > 0) stop(" file:\n", paste(miss0, collapse = "\n"))

# ============================================================
# 3) helper:region / theme / import
# ============================================================
str_to_gr <- function(region) {
  m <- str_match(region, "^(\\S+):(\\d+)-(\\d+)$")
  if (any(is.na(m))) stop(" :", region, "(  chr:start-end)")
  GRanges(seqnames = m[,2], ranges = IRanges(as.integer(m[,3]), as.integer(m[,4])))
}

# region mode
if (str_trim(opt$only_region) != "") {
  if (str_trim(opt$region_big) != "" || str_trim(opt$region_small) != "") {
    stop(" :  --only_region   --region_big/--region_small")
  }
  region_big_gr   <- str_to_gr(opt$only_region)
  region_small_gr <- region_big_gr
} else {
  if (str_trim(opt$region_big) == "") stop("  --only_region   --region_big")
  region_big_gr   <- str_to_gr(opt$region_big)
  region_small_gr <- if (str_trim(opt$region_small) == "") region_big_gr else str_to_gr(opt$region_small)
}

same_region <- identical(as.character(seqnames(region_big_gr)), as.character(seqnames(region_small_gr))) &&
  start(region_big_gr) == start(region_small_gr) && end(region_big_gr) == end(region_small_gr)

mode_dual <- (!same_region) && (str_trim(opt$only_region) == "")

# track lists
tracks_single <- parse_metric_list(opt$tracks)
big_tracks    <- parse_metric_list(opt$big_tracks)
small_tracks  <- parse_metric_list(opt$small_tracks)

if (!mode_dual) {
  if (length(tracks_single) == 0) stop("  --tracks  ")
} else {
  if (length(big_tracks) == 0) stop("  --big_tracks  ")
  if (length(small_tracks) == 0) stop("  --small_tracks  ")
}

# tissue color map
tc <- fread(tissue_color_file)
stopifnot(all(c("tissue","color") %in% colnames(tc)))
missing_t <- setdiff(tissues, tc$tissue)
if (length(missing_t) > 0) stop("tissue_colors.tsv  :", paste(missing_t, collapse = ", "))

tissue_color <- function(id) tc$color[match(id, tc$tissue)]

theme_track <- function(base_pt = 6, show_x = FALSE) {
  theme_classic(base_size = base_pt) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1,
                                  size = base_pt, margin = margin(r = 6)),
      axis.title.x = element_blank(),
      axis.text.x  = if (show_x) element_text(size = base_pt) else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.25) else element_blank(),
      axis.line.x  = if (show_x) element_line(linewidth = 0.25) else element_blank(),
      plot.margin  = margin(1, 2.2, 1, 1.8),
      legend.position = "none"
    )
}

highlight_layer <- function(hl_gr, for_region_gr) {
  if (!isTRUE(opt$highlight)) return(NULL)
  if (same_region) return(NULL)
  if (is.null(hl_gr) || length(hl_gr) == 0) return(NULL)
  if (!(hl_gr %over% for_region_gr)) return(NULL)

  xmin <- max(start(hl_gr), start(for_region_gr))
  xmax <- min(end(hl_gr),   end(for_region_gr))
  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin) return(NULL)

  annotate("rect",
           xmin = xmin, xmax = xmax,
           ymin = -Inf, ymax = Inf,
           fill = "yellow", alpha = opt$highlight_alpha, color = NA)
}

import_signal <- function(path, which_gr) {
  gr <- rtracklayer::import(path, which = which_gr)
  n <- length(gr)

  if (!"score" %in% colnames(mcols(gr))) {
    if ("value" %in% colnames(mcols(gr))) {
      v <- suppressWarnings(as.numeric(mcols(gr)$value))
      if (length(v) != n) v <- rep(0, n)
      mcols(gr)$score <- v
    } else {
      mcols(gr)$score <- rep(0, n)
    }
  } else {
    v <- suppressWarnings(as.numeric(mcols(gr)$score))
    if (length(v) != n) v <- rep(0, n)
    mcols(gr)$score <- v
  }
  mcols(gr)$score[!is.finite(mcols(gr)$score)] <- 0
  gr
}

.sig_cache <- new.env(parent = emptyenv())
cache_key <- function(tag, id, region_gr) {
  paste(tag, id, as.character(seqnames(region_gr)), start(region_gr), end(region_gr), sep="|")
}
get_signal_cached <- function(tag, id, path, region_gr) {
  k <- cache_key(tag, id, region_gr)
  if (exists(k, envir = .sig_cache, inherits = FALSE)) {
    return(get(k, envir = .sig_cache, inherits = FALSE))
  }
  gr <- import_signal(path, region_gr)
  assign(k, gr, envir = .sig_cache)
  gr
}

# GRanges(signal) -> per-bp Rle within region
gr_to_region_rle <- function(sig_gr, region_gr) {
  w <- width(region_gr)
  if (length(sig_gr) == 0) return(S4Vectors::Rle(0.0, w))
  gr <- sig_gr[sig_gr %over% region_gr]
  if (length(gr) == 0) return(S4Vectors::Rle(0.0, w))

  pstart <- pmax(start(gr), start(region_gr))
  pend   <- pmin(end(gr),   end(region_gr))

  gr2 <- GRanges(
    seqnames = "region",
    ranges   = IRanges(start = pstart - start(region_gr) + 1L,
                       end   = pend   - start(region_gr) + 1L),
    strand   = "*"
  )
  mcols(gr2)$score <- mcols(gr)$score
  mcols(gr2)$score[!is.finite(mcols(gr2)$score)] <- 0

  cov <- coverage(gr2, weight = "score", width = w)
  cov[["region"]]
}

binned_rect_df <- function(sig_gr, region_gr, tile_width, stat = c("mean","max","sum")) {
  stat <- match.arg(stat)
  w <- width(region_gr)
  rle <- gr_to_region_rle(sig_gr, region_gr)

  tile_width <- as.integer(tile_width)
  if (!is.finite(tile_width) || tile_width <= 0) tile_width <- 25L

  starts <- seq(1L, w, by = tile_width)
  ends   <- pmin(w, starts + tile_width - 1L)
  v <- Views(rle, IRanges(starts, ends))

  y <- switch(
    stat,
    mean = as.numeric(viewMeans(v)),
    max  = as.numeric(viewMaxs(v)),
    sum  = as.numeric(viewSums(v))
  )
  y[!is.finite(y)] <- 0

  xmin <- start(region_gr) + (starts - 1L)
  xmax <- start(region_gr) + (ends   - 1L)
  xmid <- (xmin + xmax) / 2

  tibble(xmin = xmin, xmax = xmax, x = xmid, y = y)
}

# ============================================================
# ============================================================
annot_tissue_right <- function(region_gr, y_top, tissue, col) {
  annotate("text",
           x = end(region_gr), y = y_top,
           label = tissue,
           hjust = 0.98, vjust = 1.2,
           size = (opt$base_pt - 0.5)/ggplot2::.pt,
           color = col, fontface = "bold")
}

plot_pos_area <- function(sig_gr, region_gr, ylab, color,
                          tile_width, stat = "mean",
                          show_x = FALSE, hl_gr = NULL,
                          tissue = NULL, tissue_col = NULL,
                          ymax_fixed = NULL) {
  df <- binned_rect_df(sig_gr, region_gr, tile_width, stat = stat)

  ymax <- if (!is.null(ymax_fixed)) {
    as.numeric(ymax_fixed)[1]
  } else if (opt$yscale == "track") {
    max(df$y, na.rm = TRUE)
  } else {
    suppressWarnings(quantile(df$y, opt$yq, na.rm = TRUE))
  }
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(df$y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  range_lab <- sprintf("[0-%.3g]", ymax)

  ggplot(df, aes(x = x, y = y)) +
    highlight_layer(hl_gr, region_gr) +
    geom_ribbon(aes(ymin = 0, ymax = pmin(y, ymax)), fill = color, color = NA) +
    geom_line(color = color, linewidth = 0.18) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)),
                    ylim = c(0, ymax), expand = FALSE, clip = "on") +
    annotate("text", x = start(region_gr), y = ymax,
             label = range_lab, hjust = 0, vjust = 1.2,
             size = (opt$base_pt - 1)/ggplot2::.pt) +
    { if (!is.null(tissue) && !is.null(tissue_col)) annot_tissue_right(region_gr, ymax, tissue, tissue_col) } +
    labs(y = ylab) +
    theme_track(opt$base_pt, show_x = show_x)
}

plot_pos_bars <- function(sig_gr, region_gr, ylab, color,
                          tile_width, stat = "mean",
                          show_x = FALSE, hl_gr = NULL,
                          tissue = NULL, tissue_col = NULL,
                          ymax_fixed = NULL) {
  df <- binned_rect_df(sig_gr, region_gr, tile_width, stat = stat)

  ymax <- if (!is.null(ymax_fixed)) {
    as.numeric(ymax_fixed)[1]
  } else if (opt$yscale == "track") {
    max(df$y, na.rm = TRUE)
  } else {
    suppressWarnings(quantile(df$y, opt$yq, na.rm = TRUE))
  }
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(df$y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  range_lab <- sprintf("[0-%.3g]", ymax)

  ggplot(df) +
    highlight_layer(hl_gr, region_gr) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = pmin(y, ymax)),
              fill = color, color = NA) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)),
                    ylim = c(0, ymax), expand = FALSE, clip = "on") +
    annotate("text", x = start(region_gr), y = ymax,
             label = range_lab, hjust = 0, vjust = 1.2,
             size = (opt$base_pt - 1)/ggplot2::.pt) +
    { if (!is.null(tissue) && !is.null(tissue_col)) annot_tissue_right(region_gr, ymax, tissue, tissue_col) } +
    labs(y = ylab) +
    theme_track(opt$base_pt, show_x = show_x)
}

plot_signed_bars <- function(sig_gr, region_gr, ylab, color,
                             tile_width, stat = "mean",
                             show_x = FALSE, hl_gr = NULL,
                             tissue = NULL, tissue_col = NULL,
                             ylim_fixed = NULL) {
  df <- binned_rect_df(sig_gr, region_gr, tile_width, stat = stat)

  if (!is.null(ylim_fixed) && length(ylim_fixed) >= 2) {
    ymin <- as.numeric(ylim_fixed)[1]
    ymax <- as.numeric(ylim_fixed)[2]
  } else if (opt$yscale == "track") {
    ymin <- min(df$y, na.rm = TRUE)
    ymax <- max(df$y, na.rm = TRUE)
  } else {
    yy <- suppressWarnings(quantile(abs(df$y), opt$yq, na.rm = TRUE))
    if (!is.finite(yy) || yy <= 0) yy <- max(abs(df$y), na.rm = TRUE)
    if (!is.finite(yy) || yy <= 0) yy <- 1
    ymin <- -yy
    ymax <-  yy
  }

  ymin <- min(ymin, 0)
  ymax <- max(ymax, 0)

  if (!is.finite(ymin) || !is.finite(ymax) || (ymin == 0 && ymax == 0)) { ymin <- -1; ymax <- 1 }
  if (ymax < ymin) { tmp <- ymin; ymin <- ymax; ymax <- tmp }
  if (ymin == ymax) { ymin <- ymin - 1; ymax <- ymax + 1 }

  range_lab <- sprintf("[%.3g-%.3g]", ymin, ymax)

  ggplot(df) +
    highlight_layer(hl_gr, region_gr) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = pmax(ymin, pmin(0, y)),
                  ymax = pmin(ymax, pmax(0, y))),
              fill = color, color = NA) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)),
                    ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    annotate("text", x = start(region_gr), y = ymax,
             label = range_lab, hjust = 0, vjust = 1.2,
             size = (opt$base_pt - 1)/ggplot2::.pt) +
    { if (!is.null(tissue) && !is.null(tissue_col)) annot_tissue_right(region_gr, ymax, tissue, tissue_col) } +
    labs(y = ylab) +
    theme_track(opt$base_pt, show_x = show_x)
}

plot_scalebar <- function(region_gr) {
  w <- width(region_gr)
  candidates <- c(10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000)
  len <- candidates[which.min(abs(candidates - w/4))]
  lab <- if (len >= 1000) paste0(len/1000, " kb") else paste0(len, " bp")

  x2 <- end(region_gr); x1 <- x2 - len; y <- 0.5
  ggplot() +
    annotate("text", x = start(region_gr), y = y,
             label = paste0(as.character(seqnames(region_gr)), ":", start(region_gr), "-", end(region_gr)),
             hjust = 0, size = opt$base_pt/ggplot2::.pt) +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.4) +
    annotate("segment", x = x1, xend = x1, y = y-0.08, yend = y+0.08, linewidth = 0.4) +
    annotate("segment", x = x2, xend = x2, y = y-0.08, yend = y+0.08, linewidth = 0.4) +
    annotate("text", x = (x1+x2)/2, y = y+0.12, label = lab,
             size = opt$base_pt/ggplot2::.pt) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)), ylim = c(0,1),
                    expand = FALSE, clip = "on") +
    theme_void() +
    theme(plot.margin = margin(0.5, 2.2, 0.5, 1.8))
}

# ============================================================
# ============================================================
make_gene_arrow_segments <- function(x0, x1, y, strand, step_bp = 250, seg_bp = 30) {
  if (!is.finite(x0) || !is.finite(x1)) return(tibble())
  if (x1 <= x0) return(tibble())

  step_bp <- as.integer(step_bp); if (!is.finite(step_bp) || step_bp <= 0) step_bp <- 250L
  seg_bp  <- as.integer(seg_bp);  if (!is.finite(seg_bp)  || seg_bp  <= 0) seg_bp  <- 30L

  st <- strand
  if (is.na(st) || st == "" || st == "*") st <- "+"

  if (st == "+") {
    tips <- seq.int(x0 + step_bp, x1 - step_bp, by = step_bp)
    if (length(tips) == 0) return(tibble())
    xend <- tips
    x    <- pmax(x0, tips - seg_bp)
  } else {
    tips <- seq.int(x1 - step_bp, x0 + step_bp, by = -step_bp)
    if (length(tips) == 0) return(tibble())
    xend <- tips
    x    <- pmin(x1, tips + seg_bp)
  }

  tibble(x = x, xend = xend, y = y, yend = y)
}

plot_one_gene <- function(gtf_gz, region_gr, hl_gr = NULL, show_x = TRUE) {
  gtf <- rtracklayer::import(gtf_gz, which = region_gr)

  feat <- NULL
  if ("type" %in% colnames(mcols(gtf))) feat <- as.character(mcols(gtf)$type)
  if (is.null(feat) && "feature" %in% colnames(mcols(gtf))) feat <- as.character(mcols(gtf)$feature)
  if (is.null(feat)) {
    return(ggplot() + highlight_layer(hl_gr, region_gr) +
             coord_cartesian(xlim=c(start(region_gr), end(region_gr)), expand=FALSE, clip="on") +
             labs(y="Genes") + theme_track(opt$base_pt, show_x))
  }

  tx <- gtf[tolower(feat) %in% c("transcript","mrna","rna")]  #   transcript/mRNA
  if (length(tx) == 0) {
    return(ggplot() + highlight_layer(hl_gr, region_gr) +
             coord_cartesian(xlim=c(start(region_gr), end(region_gr)), expand=FALSE, clip="on") +
             labs(y="Genes") + theme_track(opt$base_pt, show_x))
  }

  # exon / CDS / UTR
  ex  <- gtf[tolower(feat) == "exon"]
  cds <- gtf[tolower(feat) %in% c("cds","start_codon","stop_codon")]  #   coding
  utr <- gtf[tolower(feat) %in% c("five_prime_utr","three_prime_utr","utr")]

  tx_id <- if ("transcript_id" %in% colnames(mcols(tx))) as.character(mcols(tx)$transcript_id) else NA_character_
  gene_name <- NULL
  if ("gene" %in% colnames(mcols(tx))) gene_name <- as.character(mcols(tx)$gene)
  if (is.null(gene_name) && "gene_name" %in% colnames(mcols(tx))) gene_name <- as.character(mcols(tx)$gene_name)
  if (is.null(gene_name) && "gene_id" %in% colnames(mcols(tx))) gene_name <- as.character(mcols(tx)$gene_id)
  if (is.null(gene_name)) gene_name <- rep("Gene", length(tx))

  olw <- pmax(0L, pmin(end(tx), end(region_gr)) - pmax(start(tx), start(region_gr)) + 1L)
  best_gene <- names(which.max(tapply(olw, gene_name, sum)))
  tx_sel <- tx[gene_name == best_gene]
  tx_id_sel <- if ("transcript_id" %in% colnames(mcols(tx_sel))) as.character(mcols(tx_sel)$transcript_id) else rep(NA_character_, length(tx_sel))
  strand_sel <- as.character(strand(tx_sel))
  strand_sel[is.na(strand_sel) | strand_sel=="" | strand_sel=="*"] <- "+"

  if (length(tx_sel) > 8) {
    o <- order(olw[gene_name == best_gene], decreasing = TRUE)
    tx_sel <- tx_sel[o[1:8]]
    tx_id_sel <- tx_id_sel[o[1:8]]
    strand_sel <- strand_sel[o[1:8]]
  }

  tx_sel <- tx_sel[order(start(tx_sel), end(tx_sel))]
  tx_id_sel <- tx_id_sel[order(start(tx_sel), end(tx_sel))]
  strand_sel <- strand_sel[order(start(tx_sel), end(tx_sel))]

  nlanes <- length(tx_sel)
  lane_y <- seq(from = nlanes, to = 1)  #   UCSC: ,

  base_df <- tibble(
    lane = seq_len(nlanes),
    y    = lane_y,
    x0   = pmax(start(tx_sel), start(region_gr)),
    x1   = pmin(end(tx_sel),   end(region_gr)),
    strand = strand_sel,
    txid   = tx_id_sel
  ) %>% filter(x1 > x0)

  arrow_df <- bind_rows(lapply(seq_len(nrow(base_df)), function(i){
    arr <- make_gene_arrow_segments(
      x0 = base_df$x0[i], x1 = base_df$x1[i], y = base_df$y[i], strand = base_df$strand[i],
      step_bp = opt$gene_arrow_step_bp,
      seg_bp  = opt$gene_arrow_seg_bp
    )
    if (nrow(arr) == 0) return(tibble())
    arr$lane <- base_df$lane[i]
    arr$y    <- base_df$y[i]
    arr
  }))

  boxes_cds <- list()
  boxes_utr <- list()

  infer_utr <- function(ex_gr, cds_gr) {
    if (length(ex_gr) == 0) return(GRanges())
    ex_r <- reduce(ex_gr, ignore.strand = TRUE)
    if (length(cds_gr) == 0) return(ex_r)
    cds_r <- reduce(cds_gr, ignore.strand = TRUE)
    psetdiff(ex_r, cds_r)
  }

  for (i in seq_len(nlanes)) {
    tid <- tx_id_sel[i]
    if (!is.na(tid) && tid != "" && "transcript_id" %in% colnames(mcols(gtf))) {
      ex_i  <- ex[as.character(mcols(ex)$transcript_id) == tid]
      cds_i <- cds[as.character(mcols(cds)$transcript_id) == tid]
      utr_i <- if (length(utr) > 0) utr[as.character(mcols(utr)$transcript_id) == tid] else GRanges()
    } else {
      ex_i  <- ex[ex %over% tx_sel[i]]
      cds_i <- cds[cds %over% tx_sel[i]]
      utr_i <- if (length(utr) > 0) utr[utr %over% tx_sel[i]] else GRanges()
    }

    ex_i  <- trim(ex_i[ex_i %over% region_gr])
    cds_i <- trim(cds_i[cds_i %over% region_gr])
    utr_i <- trim(utr_i[utr_i %over% region_gr])

    if (length(utr_i) == 0) {
      utr_i <- infer_utr(ex_i, cds_i)
      utr_i <- utr_i[utr_i %over% region_gr]
    }

    if (length(cds_i) > 0) {
      boxes_cds[[i]] <- tibble(
        xmin = start(cds_i), xmax = end(cds_i),
        lane = i, y = lane_y[i]
      )
    }
    if (length(utr_i) > 0) {
      boxes_utr[[i]] <- tibble(
        xmin = start(utr_i), xmax = end(utr_i),
        lane = i, y = lane_y[i]
      )
    }
  }

  cds_df <- bind_rows(boxes_cds)
  utr_df <- bind_rows(boxes_utr)

  y_cds <- 0.22
  y_utr <- 0.12

  p <- ggplot() +
    highlight_layer(hl_gr, region_gr) +
    # intron baseline
    geom_segment(data = base_df,
                 aes(x = x0, xend = x1, y = y, yend = y),
                 linewidth = 0.25, lineend = "round") +
    # direction arrows on baseline (UCSC full mode)
    { if (nrow(arrow_df) > 0)
      geom_segment(
        data = arrow_df,
        aes(x = x, xend = xend, y = y, yend = yend),
        linewidth = 0.25,
        lineend = "round",
        arrow = grid::arrow(length = grid::unit(opt$gene_arrow_mm, "mm"), type = "closed")
      )
    } +
    # UTR thin boxes
    { if (nrow(utr_df) > 0)
      geom_rect(data = utr_df,
                aes(xmin = xmin, xmax = xmax, ymin = y - y_utr, ymax = y + y_utr),
                fill = "black", color = "black", linewidth = 0)
    } +
    # CDS thick boxes
    { if (nrow(cds_df) > 0)
      geom_rect(data = cds_df,
                aes(xmin = xmin, xmax = xmax, ymin = y - y_cds, ymax = y + y_cds),
                fill = "black", color = "black", linewidth = 0)
    } +
    # gene label (left)
    annotate("text", x = start(region_gr), y = max(lane_y) + 0.55,
             label = best_gene, hjust = 0,
             size = opt$base_pt/ggplot2::.pt) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)),
                    ylim = c(0.5, max(lane_y) + 1.0),
                    expand = FALSE, clip = "on") +
    labs(y = "Genes") +
    theme_track(opt$base_pt, show_x = show_x)

  p
}


# ============================================================
# ============================================================
plot_contrib_logo_local <- function(contrib_gr, region_gr, tissue, tissue_col,
                                    hl_gr=NULL, show_x=FALSE, ylim_fixed=NULL) {
  fa <- FaFile(fa_file)
  on.exit(try(close(fa), silent = TRUE), add = TRUE)

  seq <- as.character(getSeq(fa, region_gr)[[1]])
  bases <- strsplit(seq, "")[[1]]
  w <- length(bases)

  if (!is.finite(w) || w <= 0) {
    return(ggplot() + labs(y="Logo") + theme_track(opt$base_pt, show_x = show_x))
  }

  v <- as.numeric(gr_to_region_rle(contrib_gr, region_gr))
  v[!is.finite(v)] <- 0
  if (length(v) < w) v <- c(v, rep(0, w - length(v)))
  if (length(v) > w) v <- v[seq_len(w)]

  if (sum(abs(v)) < 1e-12) {
    if (!is.null(ylim_fixed) && length(ylim_fixed) >= 2) {
      ymin <- as.numeric(ylim_fixed)[1]
      ymax <- as.numeric(ylim_fixed)[2]
    } else {
      ymin <- -1; ymax <- 1
    }
    ymin <- min(ymin, 0); ymax <- max(ymax, 0)
    if (!is.finite(ymin) || !is.finite(ymax) || (ymin==0 && ymax==0)) { ymin <- -1; ymax <- 1 }
    if (ymax < ymin) { tmp <- ymin; ymin <- ymax; ymax <- tmp }
    if (ymin == ymax) { ymin <- ymin - 1; ymax <- ymax + 1 }

    ymin_plot <- ymin * opt$logo_pad
    ymax_plot <- ymax * opt$logo_pad
    range_lab <- sprintf("[%.3g-%.3g]", ymin, ymax)

    br <- pretty(c(1, w), n=4); br <- br[br>=1 & br<=w]
    lab <- start(region_gr) + br - 1

    return(
      ggplot() +
        highlight_layer(hl_gr, region_gr) +
        geom_hline(yintercept = 0, linewidth = 0.2) +
        coord_cartesian(xlim=c(1, w), ylim=c(ymin_plot, ymax_plot), expand=FALSE, clip="on") +
        scale_x_continuous(breaks=br, labels=lab) +
        annotate("text", x = 1, y = ymax_plot, label = range_lab,
                 hjust = 0, vjust = 1.2,
                 size = (opt$base_pt - 1)/ggplot2::.pt) +
        annotate("text", x = w, y = ymax_plot, label = tissue,
                 hjust = 0.98, vjust = 1.2,
                 size = (opt$base_pt - 0.5)/ggplot2::.pt,
                 color = tissue_col, fontface = "bold") +
        labs(y="Logo") +
        theme_track(opt$base_pt, show_x = show_x)
    )
  }

  mat <- matrix(0, nrow=4, ncol=w)
  rownames(mat) <- c("A","C","G","T")
  for (i in seq_len(w)) {
    b <- toupper(bases[i])
    if (b %in% rownames(mat)) mat[b, i] <- v[i]
  }

  if (!is.null(ylim_fixed) && length(ylim_fixed) >= 2) {
    ymin <- as.numeric(ylim_fixed)[1]
    ymax <- as.numeric(ylim_fixed)[2]
  } else if (opt$yscale == "track") {
    ymin <- min(v, na.rm = TRUE)
    ymax <- max(v, na.rm = TRUE)
  } else {
    yy <- suppressWarnings(quantile(abs(v), opt$yq, na.rm = TRUE))
    if (!is.finite(yy) || yy <= 0) yy <- max(abs(v), na.rm = TRUE)
    if (!is.finite(yy) || yy <= 0) yy <- 1
    ymin <- -yy
    ymax <-  yy
  }

  ymin <- min(ymin, 0); ymax <- max(ymax, 0)
  if (!is.finite(ymin) || !is.finite(ymax) || (ymin==0 && ymax==0)) { ymin <- -1; ymax <- 1 }
  if (ymax < ymin) { tmp <- ymin; ymin <- ymax; ymax <- tmp }
  if (ymin == ymax) { ymin <- ymin - 1; ymax <- ymax + 1 }

  ymin_plot <- ymin * opt$logo_pad
  ymax_plot <- ymax * opt$logo_pad
  range_lab <- sprintf("[%.3g-%.3g]", ymin, ymax)

  p <- ggseqlogo(mat, method="custom", seq_type="dna") +
    scale_x_continuous(limits = c(1, w), expand = c(0, 0)) +
    coord_cartesian(ylim = c(ymin_plot, ymax_plot), expand = FALSE, clip = "on") +
    labs(y = "Logo") +
    theme_track(opt$base_pt, show_x = show_x)

  p +
    annotate("text", x = 1, y = ymax_plot, label = range_lab,
             hjust = 0, vjust = 1.2,
             size = (opt$base_pt - 1)/ggplot2::.pt) +
    annotate("text", x = w, y = ymax_plot, label = tissue,
             hjust = 0.98, vjust = 1.2,
             size = (opt$base_pt - 0.5)/ggplot2::.pt,
             color = tissue_col, fontface = "bold")
}

# ============================================================
# ============================================================
pack_lanes <- function(df) {
  df <- df[order(df$start, df$end), ]
  lane_end <- numeric(0)
  lane <- integer(nrow(df))
  for (i in seq_len(nrow(df))) {
    placed <- FALSE
    for (l in seq_along(lane_end)) {
      if (df$start[i] > lane_end[l]) {
        lane[i] <- l
        lane_end[l] <- df$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      lane_end <- c(lane_end, df$end[i])
      lane[i] <- length(lane_end)
    }
  }
  df$lane <- lane
  df
}

make_chevrons_local_segments <- function(start_bp, end_bp, y0, y1, strand,
                                         len_bp = 20, step_bp = 60) {
  if (!is.finite(start_bp) || !is.finite(end_bp) || end_bp <= start_bp) return(tibble())
  st <- strand
  if (is.na(st) || st == "" || st == "*") st <- "+"

  len_bp  <- abs(as.integer(len_bp));  if (!is.finite(len_bp)  || len_bp  <= 0) len_bp  <- 20L
  step_bp <- abs(as.integer(step_bp)); if (!is.finite(step_bp) || step_bp <= 0) step_bp <- 60L

  h <- (y1 - y0) * 0.70
  ymid <- (y0 + y1) / 2
  yb0 <- ymid - h/2
  yb1 <- ymid + h/2

  if (st == "+") {
    from <- start_bp + step_bp
    to   <- end_bp - 1
    if (from > to) return(tibble())
    tips  <- seq.int(from, to, by = step_bp)
    bases <- pmax(start_bp, tips - len_bp)
  } else {
    from <- end_bp - step_bp
    to   <- start_bp + 1
    if (from < to) return(tibble())
    tips  <- seq.int(from, to, by = -step_bp)
    bases <- pmin(end_bp, tips + len_bp)
  }
  if (length(tips) == 0) return(tibble())

  out <- vector("list", length(tips))
  for (i in seq_along(tips)) {
    out[[i]] <- tibble(
      g    = i,
      seg  = c(1L, 2L),
      x    = c(bases[i], bases[i]),
      y    = c(yb0, yb1),
      xend = c(tips[i], tips[i]),
      yend = c(ymid, ymid)
    )
  }
  bind_rows(out)
}

plot_hits_blocks_local <- function(hits_gr, region_gr, show_x=TRUE) {
  w <- width(region_gr)

  h <- hits_gr[hits_gr %over% region_gr]
  if (length(h) == 0) {
    return(ggplot() +
             coord_cartesian(xlim=c(1, w), expand=FALSE, clip="on") +
             labs(y="Instances") +
             theme_track(opt$base_pt, show_x))
  }

  df <- tibble(
    start = start(h) - start(region_gr) + 1L,
    end   = end(h)   - start(region_gr) + 1L,
    strand= as.character(strand(h)),
    name  = if ("name" %in% colnames(mcols(h))) as.character(mcols(h)$name) else ""
  ) %>%
    mutate(
      start = pmax(1L, start),
      end   = pmin(as.integer(w), end),
      strand = ifelse(is.na(strand) | strand=="*" | strand=="", "+", strand)
    ) %>%
    filter(end > start) %>%
    arrange(start, end)

  if (nrow(df) == 0) {
    return(ggplot() +
             coord_cartesian(xlim=c(1, w), expand=FALSE, clip="on") +
             labs(y="Instances") +
             theme_track(opt$base_pt, show_x))
  }

  df <- pack_lanes(df) %>%
    mutate(
      y0 = lane - 0.35, y1 = lane + 0.35,
      fill = ifelse(strand=="+", "#d73027", "#4575b4")
    )

  chev <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    chevi <- make_chevrons_local_segments(
      start_bp = df$start[i], end_bp = df$end[i],
      y0 = df$y0[i], y1 = df$y1[i], strand = df$strand[i],
      len_bp = opt$hit_arrow_len, step_bp = opt$hit_arrow_step
    )
    if (nrow(chevi) > 0) chevi <- chevi %>% mutate(hit_i = i, lane = df$lane[i])
    chev[[i]] <- chevi
  }
  chev_df <- bind_rows(chev)
  nlanes <- max(df$lane)

  br <- pretty(c(1, w), n=4)
  br <- br[br>=1 & br<=w]
  lab <- start(region_gr) + br - 1

  ggplot() +
    geom_rect(data=df, aes(xmin=start, xmax=end, ymin=y0, ymax=y1),
              fill=df$fill, color=df$fill, linewidth=0) +
    { if (nrow(chev_df) > 0)
      geom_segment(
        data = chev_df,
        aes(x = x, y = y, xend = xend, yend = yend,
            group = interaction(hit_i, g, seg, lane)),
        color = "white", linewidth = 0.45, lineend = "round"
      )
    } +
    geom_text(data=df, aes(x=(start+end)/2, y=y1+0.18, label=name),
              size=opt$base_pt/ggplot2::.pt, color="black") +
    coord_cartesian(xlim=c(1, w), ylim=c(0.5, nlanes+0.9), expand=FALSE, clip="on") +
    scale_x_continuous(breaks=br, labels=lab) +
    labs(y="Instances") +
    theme_track(opt$base_pt, show_x=show_x)
}

read_hits_union <- function(tissue_ids, region_gr) {
  all <- GRanges()
  for (id in tissue_ids) {
    p <- resolve_hits_bed(id)
    if (!file.exists(p)) stop("  hits bed:tissue=", id, "\n", p)
    gr <- rtracklayer::import(p, which = region_gr)
    if (length(gr) == 0) next

    if (!"name" %in% colnames(mcols(gr))) {
      if (ncol(mcols(gr)) >= 1) mcols(gr)$name <- as.character(mcols(gr)[[1]]) else mcols(gr)$name <- ""
    }
    all <- c(all, gr)
  }
  if (length(all) == 0) return(all)

  nm <- if ("name" %in% colnames(mcols(all))) as.character(mcols(all)$name) else ""
  key <- paste0(as.character(seqnames(all)), "|", start(all), "|", end(all), "|", as.character(strand(all)), "|", nm)
  all[!duplicated(key)]
}

# ============================================================
# ============================================================
metrics_needed <- if (!mode_dual) tracks_single else unique(c(big_tracks, small_tracks))
for (t in tissues) {
  for (m in metrics_needed) {
    p <- switch(m,
      obs  = obs_bw_path(t),
      pred = pred_bw_path(t),
      nuc  = nuc_bgz_path(t),
      corr = tobias_corr_bw_path(t),
      fp   = tobias_fp_bw_path(t),
      logo = contrib_bw_path(t)
    )
    if (!file.exists(p)) stop(" file:tissue=", t, " metric=", m, "\n", p)
  }
}
for (t in hits_tissues) {
  p <- resolve_hits_bed(t)
  if (!file.exists(p)) stop("  hits file:tissue=", t, "\n", p)
}

# ============================================================
# ============================================================
region_key <- function(region_gr) {
  paste0(as.character(seqnames(region_gr)), ":", start(region_gr), "-", end(region_gr))
}
is_signed_metric <- function(metric) metric %in% c("corr","logo")

compute_metric_limits <- function(metric, region_gr, tissues_vec) {
  use_q <- (opt$yscale == "metric_q")

  if (!is_signed_metric(metric)) {
    ymax <- -Inf
    tile_w <- if (metric == "corr") opt$tile_corr else opt$tile_big # corr
    tag <- switch(metric,
      obs  = "obs",
      pred = "pred",
      nuc  = "nuc",
      fp   = "fp",
      stop("unknown metric: ", metric)
    )
    path_fun <- switch(metric,
      obs  = obs_bw_path,
      pred = pred_bw_path,
      nuc  = nuc_bgz_path,
      fp   = tobias_fp_bw_path
    )

    vals_all <- numeric(0)
    for (t in tissues_vec) {
      gr <- get_signal_cached(tag, t, path_fun(t), region_gr)
      df <- binned_rect_df(gr, region_gr, tile_w, stat = opt$bin_stat)
      if (use_q) {
        vals_all <- c(vals_all, df$y)
      } else {
        ymax <- max(ymax, max(df$y, na.rm = TRUE))
      }
    }

    if (use_q) {
      vals_all <- vals_all[is.finite(vals_all)]
      if (length(vals_all) == 0) vals_all <- 0
      ymax <- suppressWarnings(quantile(vals_all, opt$yq, na.rm = TRUE))
    }
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    return(as.numeric(ymax)[1])
  } else {
    ymin <- Inf
    ymax <- -Inf

    vals_all <- numeric(0)

    if (metric == "logo") {
      for (t in tissues_vec) {
        gr <- get_signal_cached("contrib", t, contrib_bw_path(t), region_gr)
        v  <- as.numeric(gr_to_region_rle(gr, region_gr))
        v[!is.finite(v)] <- 0
        if (use_q) {
          vals_all <- c(vals_all, v)
        } else {
          ymin <- min(ymin, min(v, na.rm = TRUE))
          ymax <- max(ymax, max(v, na.rm = TRUE))
        }
      }
    } else if (metric == "corr") {
      for (t in tissues_vec) {
        gr <- get_signal_cached("corr", t, tobias_corr_bw_path(t), region_gr)
        df <- binned_rect_df(gr, region_gr, opt$tile_corr, stat = opt$bin_stat)
        if (use_q) {
          vals_all <- c(vals_all, df$y)
        } else {
          ymin <- min(ymin, min(df$y, na.rm = TRUE))
          ymax <- max(ymax, max(df$y, na.rm = TRUE))
        }
      }
    } else {
      stop("unknown signed metric: ", metric)
    }

    if (use_q) {
      vals_all <- vals_all[is.finite(vals_all)]
      if (length(vals_all) == 0) vals_all <- 0
      ymin <- suppressWarnings(quantile(vals_all, 1 - opt$yq, na.rm = TRUE))
      ymax <- suppressWarnings(quantile(vals_all, opt$yq,     na.rm = TRUE))
    }

    ymin <- min(ymin, 0)
    ymax <- max(ymax, 0)

    if (!is.finite(ymin) || !is.finite(ymax) || (ymin == 0 && ymax == 0)) { ymin <- -1; ymax <- 1 }
    if (ymax < ymin) { tmp <- ymin; ymin <- ymax; ymax <- tmp }
    if (ymin == ymax) { ymin <- ymin - 1; ymax <- ymax + 1 }

    return(c(as.numeric(ymin)[1], as.numeric(ymax)[1]))
  }
}

.global_ylim <- list()
if (opt$yscale %in% c("metric","metric_q")) {
  combos <- list()

  add_combo <- function(metric, region_gr) {
    k <- paste(metric, region_key(region_gr), sep="|")
    if (!k %in% names(combos)) combos[[k]] <<- list(metric=metric, region=region_gr)
  }

  if (!mode_dual) {
    for (m in tracks_single) add_combo(m, region_big_gr)
  } else {
    for (m in big_tracks)   add_combo(m, region_big_gr)
    for (m in small_tracks) add_combo(m, region_small_gr)
  }

  for (k in names(combos)) {
    m  <- combos[[k]]$metric
    rg <- combos[[k]]$region
    .global_ylim[[k]] <- compute_metric_limits(m, rg, tissues)
  }
}

get_global_ylim <- function(metric, region_gr) {
  if (!opt$yscale %in% c("metric","metric_q")) return(NULL)
  k <- paste(metric, region_key(region_gr), sep="|")
  .global_ylim[[k]]
}

# ============================================================
# ============================================================
plot_metric_one_tissue <- function(metric, tissue, region_gr, hl_gr=NULL, show_x=FALSE) {
  col <- tissue_color(tissue)
  plot_pos <- if (opt$pos_style == "bars") plot_pos_bars else plot_pos_area

  yl <- get_global_ylim(metric, region_gr)  # metric/metric_q

  if (metric == "obs") {
    gr <- get_signal_cached("obs", tissue, obs_bw_path(tissue), region_gr)
    return(plot_pos(gr, region_gr, "Obs", col, opt$tile_big, opt$bin_stat, show_x, hl_gr, tissue, col,
                    ymax_fixed = if (!is.null(yl) && length(yl)==1) yl else NULL))
  }
  if (metric == "pred") {
    gr <- get_signal_cached("pred", tissue, pred_bw_path(tissue), region_gr)
    return(plot_pos(gr, region_gr, "Pred", col, opt$tile_big, opt$bin_stat, show_x, hl_gr, tissue, col,
                    ymax_fixed = if (!is.null(yl) && length(yl)==1) yl else NULL))
  }
  if (metric == "nuc") {
    gr <- get_signal_cached("nuc", tissue, nuc_bgz_path(tissue), region_gr)
    return(plot_pos(gr, region_gr, "NucleoATAC", col, opt$tile_big, opt$bin_stat, show_x, hl_gr, tissue, col,
                    ymax_fixed = if (!is.null(yl) && length(yl)==1) yl else NULL))
  }
  if (metric == "corr") {
    gr <- get_signal_cached("corr", tissue, tobias_corr_bw_path(tissue), region_gr)
    return(plot_signed_bars(gr, region_gr, "TOBIAS corrected", col, opt$tile_corr, opt$bin_stat, show_x, hl_gr, tissue, col,
                            ylim_fixed = if (!is.null(yl) && length(yl)==2) yl else NULL))
  }
  if (metric == "fp") {
    gr <- get_signal_cached("fp", tissue, tobias_fp_bw_path(tissue), region_gr)
    return(plot_pos(gr, region_gr, "TOBIAS footprint", col, opt$tile_big, opt$bin_stat, show_x, hl_gr, tissue, col,
                    ymax_fixed = if (!is.null(yl) && length(yl)==1) yl else NULL))
  }
  if (metric == "logo") {
    gr <- get_signal_cached("contrib", tissue, contrib_bw_path(tissue), region_gr)
    return(plot_contrib_logo_local(gr, region_gr, tissue, col, hl_gr = hl_gr, show_x = show_x,
                                   ylim_fixed = if (!is.null(yl) && length(yl)==2) yl else NULL))
  }
  stop(" :metric=", metric)
}

# ============================================================
# ============================================================
hl_big <- if (mode_dual && isTRUE(opt$highlight)) region_small_gr else NULL

all_tracks <- list()
heights <- numeric(0)

# 1) scale bar(always top;based on big)
all_tracks <- c(all_tracks, list(plot_scalebar(region_big_gr)))
heights <- c(heights, H("scale", opt$track_h_scale))

if (!mode_dual) {
  for (m in tracks_single) {
    for (t in tissues) {
      all_tracks <- c(all_tracks, list(plot_metric_one_tissue(m, t, region_big_gr, hl_gr=NULL, show_x=FALSE)))
      heights <- c(heights, H(m, opt$track_h))
    }
  }

  all_tracks <- c(all_tracks, list(plot_one_gene(gtf_gz, region_big_gr, hl_gr=NULL, show_x=FALSE)))
  heights <- c(heights, H("genes", opt$track_h_genes))

  hits_union <- read_hits_union(hits_tissues, region_big_gr)
  all_tracks <- c(all_tracks, list(plot_hits_blocks_local(hits_union, region_big_gr, show_x=TRUE)))
  heights <- c(heights, H("hits", opt$track_h_hits))

} else {
  for (m in big_tracks) {
    for (t in tissues) {
      all_tracks <- c(all_tracks, list(plot_metric_one_tissue(m, t, region_big_gr, hl_gr=hl_big, show_x=FALSE)))
      heights <- c(heights, H(m, opt$track_h))
    }
  }

  all_tracks <- c(all_tracks, list(plot_one_gene(gtf_gz, region_big_gr, hl_gr=hl_big, show_x=TRUE)))
  heights <- c(heights, H("genes", opt$track_h_genes))

  for (m in small_tracks) {
    for (t in tissues) {
      all_tracks <- c(all_tracks, list(plot_metric_one_tissue(m, t, region_small_gr, hl_gr=NULL, show_x=FALSE)))
      heights <- c(heights, H(m, opt$track_h))
    }
  }

  hits_union <- read_hits_union(hits_tissues, region_small_gr)
  all_tracks <- c(all_tracks, list(plot_hits_blocks_local(hits_union, region_small_gr, show_x=TRUE)))
  heights <- c(heights, H("hits", opt$track_h_hits))
}

title_txt <- if (!mode_dual) {
  glue("tissues={paste(tissues, collapse=',')} | region={as.character(seqnames(region_big_gr))}:{start(region_big_gr)}-{end(region_big_gr)} | yscale={opt$yscale}")
} else {
  glue("tissues={paste(tissues, collapse=',')} | BIG={opt$region_big} | SMALL={opt$region_small} | yscale={opt$yscale}")
}

final_plot <- patchwork::wrap_plots(all_tracks, ncol = 1, heights = heights) +
  patchwork::plot_annotation(title = title_txt)

final_plot <- final_plot & theme(plot.title = element_text(size = min(8, opt$base_pt + 2), hjust = 0.5))

pdf_h <- if (isTRUE(opt$auto_height)) sum(heights) else opt$out_h

pdf.options(useDingbats = FALSE)
ggsave(opt$outpdf, final_plot, device = "pdf",
       width = opt$out_w, height = pdf_h, units = "in")

message("Saved: ", opt$outpdf)
message("PDF size: ", opt$out_w, " x ", pdf_h, " inches (auto_height=", opt$auto_height, ")")
message("yscale: ", opt$yscale)
message("Track heights (in): ", paste(sprintf("%.2f", heights), collapse = ", "))
