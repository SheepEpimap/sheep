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
})

# ============================================================
# ============================================================
opt_list <- list(
  make_option(c("-i","--id"), type="character",
              help="tissue id (must match tissue_colors.tsv tissue column)"),
  make_option(c("--region_big"), type="character", help="big region chr:start-end"),
  make_option(c("--region_small"), type="character", default="",
              help="small region chr:start-end; if empty -> same as big"),
  make_option(c("--highlight"), type="logical", default=TRUE,
              help="highlight small region on big tracks when region_small is inside region_big"),
  make_option(c("--highlight_alpha"), type="double", default=0.25, help="highlight alpha"),

  make_option(c("-o","--outpdf"), type="character", default="track_one.pdf", help="output PDF path"),
  make_option(c("--base_pt"), type="double", default=6, help="font size (pt), 5-7 recommended"),

  make_option(c("--tile_big"),  type="integer", default=25, help="bin width for positive tracks (Obs/Pred/Nuc/Footprint)"),
  make_option(c("--tile_corr"), type="integer", default=5,  help="bin width for TOBIAS corrected (signed) track"),

  make_option(c("--bin_stat"), type="character", default="mean",
              help="bin summary: mean|max|sum (default mean)"),

  make_option(c("--pos_style"), type="character", default="area",
              help="positive-track style for Obs/Pred/Nuc/Footprint: area|bars (default area)"),

  make_option(c("--hit_arrow_len"), type="integer", default=20,
              help="V-chevron arm length in bp (local coord)"),
  make_option(c("--hit_arrow_step"), type="integer", default=60,
              help="one V-chevron every N bp (local coord)"),

  make_option(c("--out_w"), type="double", default=12, help="output PDF width (in)"),
  make_option(c("--out_h"), type="double", default=8,  help="output PDF height (in)"),

  make_option(c("--yq"), type="double", default=0.999,
              help="quantile for y-limit (use abs(y) for signed). default 0.999")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$id), !is.null(opt$region_big))
opt$bin_stat  <- tolower(opt$bin_stat)
opt$pos_style <- tolower(opt$pos_style)
if (!opt$bin_stat %in% c("mean","max","sum")) stop("--bin_stat must be one of: mean|max|sum")
if (!opt$pos_style %in% c("area","bars")) stop("--pos_style must be one of: area|bars")

# ============================================================
# ============================================================
tissue_color_file <- "/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv"

hits_bed   <- glue("/data/home/sczd644/run/zsw_chrombpnet/track/hits/{opt$id}_hits_tf.bed")
obs_bw     <- glue("/data/home/sczd644/run/zsw_chrombpnet/track/atac/{opt$id}_obs_merged.bw")
pred_bw    <- glue("/data/home/sczd644/run/zsw_chrombpnet/track/atac/{opt$id}_peaks_chrombpnet_nobias.bw")

nuc_bgz    <- glue("/data/home/sczd644/run/zsw_chrombpnet/region_ann/nucleoatac/{opt$id}/{opt$id}.occ.bedgraph.gz")
contrib_bw <- glue("/data/home/sczd644/run/zsw_chrombpnet/chrombpnet_contribs/{opt$id}_chrombpnet_contribs/{opt$id}_chrombpnet_contribs.counts_scores.bw")

tobias_corr_bw <- glue("/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/tobias/{opt$id}_corrected.bw")
tobias_fp_bw   <- glue("/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/tobias/{opt$id}_footprint.bw")

gtf_gz  <- "/data/home/sczd644/run/zsw_chrombpnet/track/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz"
fa_file <- "/data/home/sczd644/run/zsw_chrombpnet/track/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna"

if (!file.exists(hits_bed) && file.exists(paste0(hits_bed, ".gz"))) hits_bed <- paste0(hits_bed, ".gz")

need <- c(tissue_color_file, hits_bed, obs_bw, pred_bw, nuc_bgz, contrib_bw,
          tobias_corr_bw, tobias_fp_bw, gtf_gz, fa_file, paste0(fa_file, ".fai"))
miss <- need[!file.exists(need)]
if (length(miss) > 0) stop("Missing files:\n", paste(miss, collapse = "\n"))

# ============================================================
# 3) helper:region / theme / import
# ============================================================
str_to_gr <- function(region) {
  m <- str_match(region, "^(\\S+):(\\d+)-(\\d+)$")
  if (any(is.na(m))) stop("Bad region format: ", region, " (expect chr:start-end)")
  GRanges(seqnames = m[,2], ranges = IRanges(as.integer(m[,3]), as.integer(m[,4])))
}

region_big_gr   <- str_to_gr(opt$region_big)
region_small_gr <- if (nchar(opt$region_small) == 0) region_big_gr else str_to_gr(opt$region_small)

same_region <- identical(as.character(seqnames(region_big_gr)), as.character(seqnames(region_small_gr))) &&
  start(region_big_gr) == start(region_small_gr) && end(region_big_gr) == end(region_small_gr)

# tissue color
tc <- fread(tissue_color_file)
stopifnot(all(c("tissue","color") %in% colnames(tc)))
if (!opt$id %in% tc$tissue) stop("id not in tissue_colors.tsv: ", opt$id)
tissue_col <- tc$color[match(opt$id, tc$tissue)]

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
      plot.margin  = margin(1, 1.5, 1, 1.5),
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
  if (!"score" %in% colnames(mcols(gr))) {
    if ("value" %in% colnames(mcols(gr))) mcols(gr)$score <- mcols(gr)$value else mcols(gr)$score <- 0
  }
  mcols(gr)$score <- suppressWarnings(as.numeric(mcols(gr)$score))
  mcols(gr)$score[!is.finite(mcols(gr)$score)] <- 0
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

# bin summary -> rect df (xmin/xmax for bar-like tracks)
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

plot_pos_area <- function(sig_gr, region_gr, ylab, color,
                          tile_width, stat = "mean",
                          show_x = FALSE, hl_gr = NULL) {
  df <- binned_rect_df(sig_gr, region_gr, tile_width, stat = stat)

  ymax <- suppressWarnings(quantile(df$y, opt$yq, na.rm = TRUE))
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(df$y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  range_lab <- sprintf("[0-%.3g]", ymax)

  ggplot(df, aes(x = x, y = y)) +
    highlight_layer(hl_gr, region_gr) +
    geom_area(fill = color, color = NA, linewidth = 0) +
    geom_line(color = color, linewidth = 0.18) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)),
                    ylim = c(0, ymax), expand = FALSE) +
    annotate("text", x = start(region_gr), y = ymax,
             label = range_lab, hjust = 0, vjust = 1.2,
             size = (opt$base_pt - 1)/ggplot2::.pt) +
    labs(y = ylab) +
    theme_track(opt$base_pt, show_x = show_x)
}

plot_pos_bars <- function(sig_gr, region_gr, ylab, color,
                          tile_width, stat = "mean",
                          show_x = FALSE, hl_gr = NULL) {
  df <- binned_rect_df(sig_gr, region_gr, tile_width, stat = stat)

  ymax <- suppressWarnings(quantile(df$y, opt$yq, na.rm = TRUE))
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(df$y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  range_lab <- sprintf("[0-%.3g]", ymax)

  ggplot(df) +
    highlight_layer(hl_gr, region_gr) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = pmin(y, ymax)),
              fill = color, color = NA) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)),
                    ylim = c(0, ymax), expand = FALSE) +
    annotate("text", x = start(region_gr), y = ymax,
             label = range_lab, hjust = 0, vjust = 1.2,
             size = (opt$base_pt - 1)/ggplot2::.pt) +
    labs(y = ylab) +
    theme_track(opt$base_pt, show_x = show_x)
}

plot_signed_bars <- function(sig_gr, region_gr, ylab, color,
                             tile_width, stat = "mean",
                             show_x = FALSE, hl_gr = NULL) {
  df <- binned_rect_df(sig_gr, region_gr, tile_width, stat = stat)

  ymax <- suppressWarnings(quantile(abs(df$y), opt$yq, na.rm = TRUE))
  if (!is.finite(ymax) || ymax <= 0) ymax <- max(abs(df$y), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  range_lab <- sprintf("[%.3g-%.3g]", -ymax, ymax)

  ggplot(df) +
    highlight_layer(hl_gr, region_gr) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = pmax(-ymax, pmin(0, y)),
                  ymax = pmin( ymax, pmax(0, y))),
              fill = color, color = NA) +
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)),
                    ylim = c(-ymax, ymax), expand = FALSE) +
    annotate("text", x = start(region_gr), y = ymax,
             label = range_lab, hjust = 0, vjust = 1.2,
             size = (opt$base_pt - 1)/ggplot2::.pt) +
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
    coord_cartesian(xlim = c(start(region_gr), end(region_gr)), ylim = c(0,1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0.5, 1.5, 0.5, 1.5))
}

# ============================================================
# ============================================================
make_gene_arrows <- function(x0, x1, y, strand, step = 250, aw = 60, ah = 0.16) {
  if (!is.finite(x0) || !is.finite(x1)) return(tibble())
  if (x1 <= x0) return(tibble())

  step <- as.integer(step)
  aw   <- as.integer(aw)
  if (!is.finite(step) || step <= 0) step <- 250L
  if (!is.finite(aw)   || aw   <= 0) aw   <- 60L

  tips <- if (strand == "+") seq.int(x0 + step, x1 - step, by = step) else seq.int(x1 - step, x0 + step, by = -step)
  if (length(tips) == 0) return(tibble())

  out <- vector("list", length(tips))
  for (i in seq_along(tips)) {
    tip <- tips[i]
    if (strand == "+") {
      base <- max(x0, tip - aw)
      out[[i]] <- tibble(g=i, x=c(base, base, tip), y=c(y-ah, y+ah, y))
    } else {
      base <- min(x1, tip + aw)
      out[[i]] <- tibble(g=i, x=c(base, base, tip), y=c(y-ah, y+ah, y))
    }
  }
  bind_rows(out)
}

plot_one_gene <- function(gtf_gz, region_gr, hl_gr = NULL, show_x = TRUE) {
  gtf <- rtracklayer::import(gtf_gz, which = region_gr)

  feat <- NULL
  if ("type" %in% colnames(mcols(gtf))) feat <- mcols(gtf)$type
  if (is.null(feat) && "feature" %in% colnames(mcols(gtf))) feat <- mcols(gtf)$feature
  if (is.null(feat)) {
    return(ggplot() + highlight_layer(hl_gr, region_gr) +
             coord_cartesian(xlim=c(start(region_gr), end(region_gr)), expand=FALSE) +
             labs(y="Genes") + theme_track(opt$base_pt, show_x))
  }

  tx <- gtf[feat == "transcript"]
  ex <- gtf[feat == "exon"]
  if (length(tx) == 0 || length(ex) == 0) {
    return(ggplot() + highlight_layer(hl_gr, region_gr) +
             coord_cartesian(xlim=c(start(region_gr), end(region_gr)), expand=FALSE) +
             labs(y="Genes") + theme_track(opt$base_pt, show_x))
  }

  ov <- findOverlaps(tx, region_gr, ignore.strand = TRUE)
  if (length(ov) == 0) {
    return(ggplot() + highlight_layer(hl_gr, region_gr) +
             coord_cartesian(xlim=c(start(region_gr), end(region_gr)), expand=FALSE) +
             labs(y="Genes") + theme_track(opt$base_pt, show_x))
  }

  qh <- queryHits(ov)
  tx_s <- start(tx[qh]); tx_e <- end(tx[qh])
  rg_s <- start(region_gr); rg_e <- end(region_gr)
  ow <- pmax(0L, pmin(tx_e, rg_e) - pmax(tx_s, rg_s) + 1L)
  pick_tx <- tx[qh[which.max(ow)]]

  tx_id <- if ("transcript_id" %in% colnames(mcols(pick_tx))) mcols(pick_tx)$transcript_id else NA_character_
  gene  <- if ("gene" %in% colnames(mcols(pick_tx))) as.character(mcols(pick_tx)$gene) else as.character(mcols(pick_tx)$gene_id)
  st    <- as.character(strand(pick_tx))
  if (is.na(st) || st == "" || st == "*") st <- "+"

  ex_pick <- if (!is.na(tx_id)) ex[mcols(ex)$transcript_id == tx_id] else ex[ex %over% pick_tx]
  exon_df <- tibble(start = start(ex_pick), end = end(ex_pick)) %>% arrange(start, end)

  y <- 1; eh <- 0.22
  arr <- make_gene_arrows(start(pick_tx), end(pick_tx), y, st)

  ggplot() +
    highlight_layer(hl_gr, region_gr) +
    geom_segment(aes(x=start(pick_tx), xend=end(pick_tx), y=y, yend=y), linewidth=0.25) +
    geom_rect(data=exon_df, aes(xmin=start, xmax=end, ymin=y-eh, ymax=y+eh),
              fill="black", color="black", linewidth=0) +
    { if (nrow(arr) > 0)
        geom_polygon(data=arr, aes(x=x, y=y, group=g),
                     fill="black", color="black", linewidth=0)
    } +
    annotate("text", x=start(region_gr), y=y+0.45, label=gene,
             hjust=0, size=opt$base_pt/ggplot2::.pt) +
    coord_cartesian(xlim=c(start(region_gr), end(region_gr)),
                    ylim=c(0.3, 1.7), expand=FALSE) +
    labs(y="Genes") +
    theme_track(opt$base_pt, show_x=show_x)
}

# ============================================================
# ============================================================
plot_contrib_logo_local <- function(contrib_gr, region_gr) {
  fa <- FaFile(fa_file)
  on.exit(try(close(fa), silent = TRUE), add = TRUE)

  seq <- as.character(getSeq(fa, region_gr)[[1]])
  bases <- strsplit(seq, "")[[1]]
  w <- length(bases)

  if (!is.finite(w) || w <= 0) {
    return(ggplot() + labs(y="Score") + theme_track(opt$base_pt, show_x = FALSE))
  }

  v <- as.numeric(gr_to_region_rle(contrib_gr, region_gr))
  v[!is.finite(v)] <- 0
  if (length(v) < w) v <- c(v, rep(0, w - length(v)))
  if (length(v) > w) v <- v[seq_len(w)]

  mat <- matrix(0, nrow=4, ncol=w); rownames(mat) <- c("A","C","G","T")
  for (i in seq_len(w)) {
    b <- toupper(bases[i])
    if (b %in% rownames(mat)) mat[b, i] <- v[i]
  }

  ggseqlogo(mat, method="custom", seq_type="dna") +
    scale_x_continuous(limits = c(1, w), expand = c(0, 0)) +
    labs(y = "Score") +
    theme_track(opt$base_pt, show_x = FALSE)
}

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
             coord_cartesian(xlim=c(1, w), expand=FALSE) +
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
             coord_cartesian(xlim=c(1, w), expand=FALSE) +
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
    coord_cartesian(xlim=c(1, w), ylim=c(0.5, nlanes+0.9), expand=FALSE) +
    scale_x_continuous(breaks=br, labels=lab) +
    labs(y="Instances") +
    theme_track(opt$base_pt, show_x=show_x)
}

# ============================================================
# 7) load tracks
# ============================================================
hl <- if (isTRUE(opt$highlight) && !same_region) region_small_gr else NULL

obs_big  <- import_signal(obs_bw,         region_big_gr)
pred_big <- import_signal(pred_bw,        region_big_gr)
nuc_big  <- import_signal(nuc_bgz,        region_big_gr)
corr_big <- import_signal(tobias_corr_bw, region_big_gr)   # signed
fp_big   <- import_signal(tobias_fp_bw,   region_big_gr)

contrib_small <- import_signal(contrib_bw, region_small_gr)
hits_small    <- rtracklayer::import(hits_bed, which = region_small_gr)

# ============================================================
# 8) build plots
# ============================================================
p_scale <- plot_scalebar(region_big_gr)

plot_pos <- if (opt$pos_style == "bars") plot_pos_bars else plot_pos_area

p_obs <- plot_pos(
  sig_gr = obs_big, region_gr = region_big_gr,
  ylab = "Obs", color = tissue_col,
  tile_width = opt$tile_big, stat = opt$bin_stat,
  show_x = FALSE, hl_gr = hl
)

p_pred <- plot_pos(
  sig_gr = pred_big, region_gr = region_big_gr,
  ylab = "Pred", color = tissue_col,
  tile_width = opt$tile_big, stat = opt$bin_stat,
  show_x = FALSE, hl_gr = hl
)

p_nuc <- plot_pos(
  sig_gr = nuc_big, region_gr = region_big_gr,
  ylab = "NucleoATAC", color = tissue_col,
  tile_width = opt$tile_big, stat = opt$bin_stat,
  show_x = FALSE, hl_gr = hl
)

p_cor <- plot_signed_bars(
  sig_gr = corr_big, region_gr = region_big_gr,
  ylab = "TOBIAS corrected", color = tissue_col,
  tile_width = opt$tile_corr, stat = opt$bin_stat,
  show_x = FALSE, hl_gr = hl
)

p_fp <- plot_pos(
  sig_gr = fp_big, region_gr = region_big_gr,
  ylab = "TOBIAS footprint", color = tissue_col,
  tile_width = opt$tile_big, stat = opt$bin_stat,
  show_x = FALSE, hl_gr = hl
)

p_gene <- plot_one_gene(gtf_gz, region_big_gr, hl_gr = hl, show_x = !same_region)

p_logo <- plot_contrib_logo_local(contrib_small, region_small_gr)
p_hits <- plot_hits_blocks_local(hits_small, region_small_gr, show_x = TRUE)

# ============================================================
# ============================================================
all_tracks <- list(
  p_scale,
  p_obs, p_pred, p_nuc,
  p_cor, p_fp,
  p_gene,
  p_logo, p_hits
)

heights <- c(
  0.55,              # scale
  1, 1, 1,           # obs/pred/nuc
  1, 1,              # corrected/footprint
  0.95,              # genes
  0.95, 0.95         # score/instances( :  2)
)

final_plot <- patchwork::wrap_plots(
  all_tracks, ncol = 1,
  heights = heights
) +
  patchwork::plot_annotation(
    title = glue("{opt$id} | BIG: {opt$region_big} | SMALL: {if (same_region) opt$region_big else opt$region_small}")
  )

final_plot <- final_plot & theme(
  plot.title = element_text(size = min(8, opt$base_pt + 2), hjust = 0.5)
)

pdf.options(useDingbats = FALSE)
ggsave(opt$outpdf, final_plot, device = "pdf",
       width = opt$out_w, height = opt$out_h, units = "in")
message("Saved: ", opt$outpdf)
