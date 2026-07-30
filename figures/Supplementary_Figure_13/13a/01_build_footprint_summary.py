#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

from matplotlib import font_manager as fm


def configure_editable_text(font_ttf: str = "") -> str:
    """
      TrueType  ,  PDF/SVG output :
      - pdf.fonttype=42: TrueType/CID TrueType(  AI/InkScape  )
      - svg.fonttype="none": SVG   <text>   path( )
      PDF  , output  SVG   PDF.
    """
    matplotlib.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "pdf.use14corefonts": False,
        "text.usetex": False,
        "svg.fonttype": "none",
    })

    if font_ttf:
        if not os.path.exists(font_ttf):
            raise FileNotFoundError(f"--font_ttf not found: {font_ttf}")
        fm.fontManager.addfont(font_ttf)
        family = fm.FontProperties(fname=font_ttf).get_name()
    else:
        dejavu = os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans.ttf")
        if os.path.exists(dejavu):
            fm.fontManager.addfont(dejavu)
            family = "DejaVu Sans"
        else:
            family = matplotlib.rcParams.get("font.family", "sans-serif")

    matplotlib.rcParams["font.family"] = family

    print(f"[INFO] font.family   = {matplotlib.rcParams['font.family']}")
    print(f"[INFO] pdf.fonttype  = {matplotlib.rcParams['pdf.fonttype']}")
    print(f"[INFO] svg.fonttype  = {matplotlib.rcParams['svg.fonttype']}")

    return family


import matplotlib.pyplot as plt


# -------------------------
# Interval helpers (promoter/enhancer overlap)
# -------------------------
def load_bed_intervals(bed_path: str):
    """
      bed: chr start end
      dict: chrom -> (starts, ends, prefix_max_end)

    prefix_max_end + searchsorted  :
      query [qs, qe),  start < qe   idx,
      idx   prefix_max_end > qs, .
    """
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        usecols=[0, 1, 2], names=["chr", "start", "end"]
    )
    df["chr"] = df["chr"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    out = {}
    for chrom, sub in df.groupby("chr", sort=False):
        sub = sub.sort_values("start")
        s = sub["start"].to_numpy(dtype=np.int64)
        e = sub["end"].to_numpy(dtype=np.int64)
        pref = np.maximum.accumulate(e)
        out[chrom] = (s, e, pref)
    return out


def any_overlap(chrom_intervals, chrom: str, qs: np.ndarray, qe: np.ndarray) -> np.ndarray:
    """vectorized overlap for one chromosome"""
    if chrom not in chrom_intervals:
        return np.zeros(qs.shape[0], dtype=bool)

    starts, ends, pref = chrom_intervals[chrom]
    idx = np.searchsorted(starts, qe, side="left") - 1

    ok = idx >= 0
    out = np.zeros(qs.shape[0], dtype=bool)
    if not np.any(ok):
        return out

    idx_ok = idx[ok]
    out[ok] = pref[idx_ok] > qs[ok]
    return out


# -------------------------
# TSS nearest distance per chrom
# -------------------------
def nearest_tss_distance_per_chrom(centers: np.ndarray, tss_sorted: np.ndarray) -> np.ndarray:
    """
    centers: int64
    tss_sorted: sorted int64
    return float distances; use inf for invalid side to avoid int overflow
    """
    if tss_sorted is None or len(tss_sorted) == 0:
        return np.full(centers.shape[0], np.nan, dtype=float)

    idx = np.searchsorted(tss_sorted, centers, side="left")

    idx_r = np.clip(idx, 0, len(tss_sorted) - 1)
    dist_r = np.abs(tss_sorted[idx_r] - centers).astype(float)

    idx_l = np.clip(idx - 1, 0, len(tss_sorted) - 1)
    dist_l = np.abs(centers - tss_sorted[idx_l]).astype(float)

    # boundary fix
    dist_l[idx == 0] = np.inf
    dist_r[idx == len(tss_sorted)] = np.inf

    return np.minimum(dist_l, dist_r)


def read_tissues_from_file(path: str):
    with open(path) as f:
        return [x.strip() for x in f if x.strip() and not x.startswith("#")]


def main():
    ap = argparse.ArgumentParser(
        description="All-tissues pooled B/C/D plot with motifs on Y-axis (vertical ordering)."
    )
    ap.add_argument("--hits_root", required=True,
                    help="Root dir: {hits_root}/{tissue}/{tissue}_footprint_hits.tsv")
    ap.add_argument("--annotation_dir", required=True,
                    help="Dir containing AAGs_38/{tissue}_promoter.bed and {tissue}_enhancer.bed")
    ap.add_argument("--tss_bed", required=True,
                    help="TSS bed, at least 2 cols: chr, tss_start")

    ap.add_argument("--tissues", nargs="*", default=None,
                    help="Explicit tissue list, e.g. --tissues abomasum rumen cecum")
    ap.add_argument("--tissue_list", default=None,
                    help="A txt file, one tissue per line (preferred for many tissues)")

    ap.add_argument("--p_cut", type=float, default=0.05,
                    help="Footprint pvalue cutoff (default 0.05)")

    ap.add_argument("--exclude_unresolved", action="store_true",
                    help="Drop motifs starting with 'unresolved#'")

    ap.add_argument("--out_prefix", default="ALL_TISSUES",
                    help="Output prefix for pdf/png/svg/tsv")

    ap.add_argument("--ytick_step", type=int, default=1,
                    help="Show every Nth motif label on y-axis (default 1 = show all)")

    ap.add_argument("--height_per_motif", type=float, default=0.18,
                    help="Figure height per motif (inch). Increase if labels overlap.")

    ap.add_argument("--top_n", type=int, default=0,
                    help="Plot only top N motifs after sorting (0 = all).")

    ap.add_argument("--font_ttf", default="",
                    help="Optional: path to a .ttf font file to embed as editable text in PDF.")

    args = ap.parse_args()

    configure_editable_text(args.font_ttf)

    # ---- tissue list
    if args.tissues:
        tissues = args.tissues
    elif args.tissue_list:
        tissues = read_tissues_from_file(args.tissue_list)
    else:
        tissues = sorted([d for d in os.listdir(args.hits_root)
                          if os.path.isdir(os.path.join(args.hits_root, d))])

    if len(tissues) == 0:
        raise ValueError("No tissues found. Provide --tissues/--tissue_list or check --hits_root.")

    print(f"[INFO] tissues = {len(tissues)}")

    # ---- bins for B/C
    bins = [0, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, np.inf]
    bin_labels = [
        "0-200bp", "200-500bp", "500bp-1kb", "1-2kb", "2-5kb",
        "5-10kb", "10-20kb", "20-50kb", "50-100kb", ">100kb"
    ]

    # ---- load TSS once
    tss = pd.read_csv(args.tss_bed, sep="\t", header=None, usecols=[0, 1], names=["chr", "tss_start"])
    tss["chr"] = tss["chr"].astype(str)
    tss["tss_start"] = tss["tss_start"].astype(int)

    tss_groups = {}
    for chrom, sub in tss.groupby("chr", sort=False):
        tss_groups[chrom] = np.sort(sub["tss_start"].to_numpy(dtype=np.int64))

    print(f"[INFO] TSS loaded: {tss.shape[0]}")

    # ---- global accumulators
    hits_bin_all = None
    pass_bin_all = None
    pass_tot_all = None
    region_cnt_all = None

    # ---- process each tissue
    for tissue in tissues:
        hits_file = os.path.join(args.hits_root, tissue, f"{tissue}_footprint_hits.tsv")
        promoter_bed = os.path.join(args.annotation_dir, "AAGs_38", f"{tissue}_promoter.bed")
        enhancer_bed = os.path.join(args.annotation_dir, "AAGs_38", f"{tissue}_enhancer.bed")

        if not os.path.exists(hits_file):
            print(f"[WARN] missing hits: {hits_file} (skip)")
            continue
        if not (os.path.exists(promoter_bed) and os.path.exists(enhancer_bed)):
            print(f"[WARN] missing promoter/enhancer bed for {tissue} (skip)")
            continue

        prom = load_bed_intervals(promoter_bed)
        enh = load_bed_intervals(enhancer_bed)

        df = pd.read_csv(
            hits_file, sep="\t", header=0,
            usecols=["chr", "start", "end", "motif_name", "pvalue"]
        )
        df["chr"] = df["chr"].astype(str)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["center"] = ((df["start"] + df["end"]) // 2).astype(int)

        if args.exclude_unresolved:
            df = df[~df["motif_name"].astype(str).str.startswith("unresolved#")].copy()

        if df.shape[0] == 0:
            print(f"[INFO] {tissue}: no rows after filtering (skip)")
            continue

        # nearest TSS distance per chrom
        dist = np.full(df.shape[0], np.nan, dtype=float)
        for chrom, idxs in df.groupby("chr", sort=False).indices.items():
            if chrom not in tss_groups:
                continue
            idxs = np.asarray(idxs, dtype=np.int64)
            centers = df.iloc[idxs]["center"].to_numpy(dtype=np.int64)
            dist[idxs] = nearest_tss_distance_per_chrom(centers, tss_groups[chrom])
        df["dist_to_tss"] = dist
        df = df.dropna(subset=["dist_to_tss"]).copy()

        # binning
        df["bin"] = pd.cut(
            df["dist_to_tss"],
            bins=bins,
            labels=bin_labels,
            right=False,
            include_lowest=True
        )
        df = df.dropna(subset=["bin"]).copy()

        if df.shape[0] == 0:
            print(f"[INFO] {tissue}: all rows dropped after dist/bin (skip)")
            continue

        # hits_bin counts (denominator for C)
        hits_bin = df.groupby(["motif_name", "bin"], observed=True).size()
        hits_bin_all = hits_bin if hits_bin_all is None else hits_bin_all.add(hits_bin, fill_value=0)

        # footprint hits (numerator for B/C/D)
        df_pass = df[df["pvalue"] < args.p_cut].copy().reset_index(drop=True)
        if df_pass.shape[0] == 0:
            print(f"[INFO] {tissue}: footprint_hits=0 under p<{args.p_cut}")
            continue

        pass_bin = df_pass.groupby(["motif_name", "bin"], observed=True).size()
        pass_tot = df_pass.groupby("motif_name").size()
        pass_bin_all = pass_bin if pass_bin_all is None else pass_bin_all.add(pass_bin, fill_value=0)
        pass_tot_all = pass_tot if pass_tot_all is None else pass_tot_all.add(pass_tot, fill_value=0)

        # region3 annotation on df_pass
        region = np.full(df_pass.shape[0], "distal", dtype=object)

        for chrom, idxs in df_pass.groupby("chr", sort=False).groups.items():
            idxs = np.asarray(idxs, dtype=np.int64)
            qs = df_pass.loc[idxs, "start"].to_numpy(dtype=np.int64)
            qe = df_pass.loc[idxs, "end"].to_numpy(dtype=np.int64)

            is_prom = any_overlap(prom, chrom, qs, qe)
            region[idxs[is_prom]] = "promoter"

            remaining = ~is_prom
            if np.any(remaining):
                is_enh = any_overlap(enh, chrom, qs[remaining], qe[remaining])
                region[idxs[remaining][is_enh]] = "enhancer"

        df_pass["region3"] = region
        region_cnt = df_pass.groupby(["motif_name", "region3"], observed=True).size()
        region_cnt_all = region_cnt if region_cnt_all is None else region_cnt_all.add(region_cnt, fill_value=0)

        print(f"[OK] {tissue}: hits={df.shape[0]} pass={df_pass.shape[0]} motifs_pass={pass_tot.shape[0]}")

    if pass_tot_all is None or pass_tot_all.shape[0] == 0:
        raise ValueError("No footprint hits aggregated. Check p_cut, paths, and tissues.")

    # ---- build B/C/D matrices (ALL motifs)
    B = (
        pass_bin_all.div(pass_tot_all, level=0)
        .unstack("bin")
        .reindex(columns=bin_labels)
        .fillna(0.0)
    )

    C = (
        (pass_bin_all / hits_bin_all)
        .unstack("bin")
        .reindex(columns=bin_labels)
        .astype(float)
    )

    D_cnt = region_cnt_all.unstack("region3").fillna(0.0)
    for col in ["promoter", "enhancer", "distal"]:
        if col not in D_cnt.columns:
            D_cnt[col] = 0.0
    D_cnt = D_cnt[["promoter", "enhancer", "distal"]]
    D = (D_cnt.T / D_cnt.sum(axis=1)).T.fillna(0.0)

    # ---- motif order (vertical sorting): by B[0-200bp], if 0 then B[200-500bp]
    key1 = B["0-200bp"].to_numpy()
    key2 = B["200-500bp"].to_numpy()
    key_final = np.where(key1 == 0, key2, key1)

    motif_order = B.index.to_numpy()
    motif_order = motif_order[np.argsort(-key_final)]  # descending

    # align full tables
    B = B.reindex(index=motif_order).fillna(0.0)
    C = C.reindex(index=motif_order)
    D = D.reindex(index=motif_order).fillna(0.0)

    B.to_csv(f"{args.out_prefix}.B_share.tsv", sep="\t")
    C.to_csv(f"{args.out_prefix}.C_pass_over_hits.tsv", sep="\t")
    D.to_csv(f"{args.out_prefix}.D_region_share.tsv", sep="\t")

    # ---- decide topN for plotting
    if args.top_n and args.top_n > 0:
        top_n = min(args.top_n, len(motif_order))
        motif_order_plot = motif_order[:top_n]
    else:
        motif_order_plot = motif_order

    Bp = B.loc[motif_order_plot]
    Cp = C.loc[motif_order_plot]
    Dp = D.loc[motif_order_plot]

    # ---- plotting: motifs on Y (vertical), three panels in one row
    n_m = len(motif_order_plot)
    y = np.arange(n_m)

    fig_h = max(8.0, n_m * args.height_per_motif)
    fig = plt.figure(figsize=(22, fig_h))

    gs = fig.add_gridspec(nrows=1, ncols=4, width_ratios=[4.2, 3.0, 2.2, 0.25], wspace=0.25)

    axB = fig.add_subplot(gs[0, 0])
    axC = fig.add_subplot(gs[0, 1])
    axD = fig.add_subplot(gs[0, 2])
    axCbar = fig.add_subplot(gs[0, 3])

    ylow, yhigh = -0.5, n_m - 0.5

    # ----- Panel B: horizontal stacked bars -----
    left = np.zeros(n_m, dtype=float)
    handles_B = []
    for lab in bin_labels:
        vals = Bp[lab].to_numpy(dtype=float)
        h = axB.barh(y, vals, left=left, height=1.0, align="center", label=lab)
        left += vals
        handles_B.append(h[0])

    axB.set_xlim(0, 1.0)
    axB.set_xlabel("Fraction (sum=1 per motif)")
    axB.set_title(f"B: Footprint distribution across TSS distance bins (pooled, p<{args.p_cut})")

    step = max(1, args.ytick_step)
    yt = y[::step]
    axB.set_yticks(yt)
    axB.set_yticklabels(motif_order_plot[::step])

    axB.set_ylim(yhigh, ylow)
    axB.margins(y=0)
    axB.legend(handles_B, bin_labels, title="TSS bin", loc="upper right", frameon=True)

    # ----- Panel C: heatmap -----
    C_mat = Cp[bin_labels].to_numpy(dtype=float)  # n_m x n_bins
    C_masked = np.ma.masked_invalid(C_mat)

    cmap = plt.cm.Reds.copy()
    cmap.set_bad(color="white")

    finite_vals = C_mat[np.isfinite(C_mat)]
    vmax = max(np.quantile(finite_vals, 0.95), 1e-6) if finite_vals.size > 0 else 1.0

    im = axC.imshow(
        C_masked,
        aspect="auto",
        origin="upper",
        interpolation="nearest",
        vmin=0.0, vmax=vmax,
        cmap=cmap,
        extent=(-0.5, len(bin_labels) - 0.5, yhigh, ylow)
    )

    axC.set_title(f"C: Footprint proportion within each bin (pass/hits, pooled, p<{args.p_cut})")
    axC.set_xticks(np.arange(len(bin_labels)))
    axC.set_xticklabels(bin_labels, rotation=45, ha="right")
    axC.set_yticks([])
    axC.set_ylim(yhigh, ylow)

    cb = fig.colorbar(im, cax=axCbar)
    cb.set_label("Footprint proportion (p-cut)")

    # ----- Panel D: horizontal stacked bars -----
    left = np.zeros(n_m, dtype=float)
    for lab in ["promoter", "enhancer", "distal"]:
        vals = Dp[lab].to_numpy(dtype=float)
        axD.barh(y, vals, left=left, height=1.0, align="center", label=lab)
        left += vals

    axD.set_xlim(0, 1.0)
    axD.set_xlabel("Fraction (sum=1 per motif)")
    axD.set_title(f"D: Region distribution of footprint hits (pooled, p<{args.p_cut})")
    axD.set_yticks([])
    axD.set_ylim(yhigh, ylow)
    axD.margins(y=0)
    axD.legend(title="Region", loc="upper right", frameon=True)

    fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.06, wspace=0.25)

    # ---- outputs
    tag = f".top{len(motif_order_plot)}" if (args.top_n and args.top_n > 0) else ""
    out_png = f"{args.out_prefix}{tag}.BCD_vertical.p{args.p_cut}.png"
    out_pdf = f"{args.out_prefix}{tag}.BCD_vertical.p{args.p_cut}.pdf"
    out_svg = f"{args.out_prefix}{tag}.BCD_vertical.p{args.p_cut}.svg"

    configure_editable_text(args.font_ttf)

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")   # ✅   PDF( : )
    fig.savefig(out_svg, bbox_inches="tight")   # ✅  :SVG
    plt.close(fig)

    print(f"[OK] saved figure: {out_png}")
    print(f"[OK] saved figure: {out_pdf}")
    print(f"[OK] saved figure: {out_svg}")
    print(f"[OK] saved tables: {args.out_prefix}.B_share.tsv / .C_pass_over_hits.tsv / .D_region_share.tsv")


if __name__ == "__main__":
    main()
