#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import re
import argparse
import numpy as np
import pandas as pd

# ============================================================
# ============================================================
import matplotlib as mpl
mpl.use("Agg")

mpl.rcParams["pdf.fonttype"] = 42  #
mpl.rcParams["ps.fonttype"]  = 42  #

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

import matplotlib.pyplot as plt
from scipy.stats import pearsonr


# ============================================================
# ============================================================

def parse_region(region: str):
    """Parse 'chr:start-end'."""
    region_main = region.strip().split()[0]
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", region_main)
    if not m:
        raise ValueError(f"--region must be chr:start-end, got: {region}")
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    if start > end:
        start, end = end, start
    return chrom, start, end, region_main


def read_tissues(tissue_file: str):
    tissues = []
    with open(tissue_file, "r") as f:
        for line in f:
            t = line.strip()
            if t:
                tissues.append(t)
    if not tissues:
        raise RuntimeError(f"No tissues found in {tissue_file}")
    return tissues


def read_tissue_colors(color_tsv: str):
    """
    tissue_colors.tsv:
      tissue  color
      cerebral-cortex #dcd71a
    """
    df = pd.read_csv(color_tsv, sep=r"\s+", engine="python")
    df.columns = [c.lower() for c in df.columns]
    if "tissue" not in df.columns or "color" not in df.columns:
        raise RuntimeError("tissue_colors.tsv must contain columns: tissue, color")
    return dict(zip(df["tissue"].astype(str), df["color"].astype(str)))


def parse_peak_coords(name_series: pd.Series):
    """
    Accept:
      chr1_123_456...
      chr1:123-456...
    Convert ':','-' -> '_' then take first 3 parts.
    """
    s = (
        name_series.astype(str)
        .str.replace(":", "_", regex=False)
        .str.replace("-", "_", regex=False)
    )
    parts = s.str.split("_", n=3, expand=True)
    chrom = parts[0]
    start = pd.to_numeric(parts[1], errors="coerce")
    end = pd.to_numeric(parts[2], errors="coerce")
    start2 = np.minimum(start, end)
    end2 = np.maximum(start, end)
    return chrom, start2, end2


def format_p_x10(p: float):
    """Format p-value like 6.84×10^-10."""
    if p < 2.2e-16:
        return "<2.2×10^-16"
    expn = int(np.floor(np.log10(p)))
    base = p / (10 ** expn)
    base = float(f"{base:.3g}")
    return f"{base}×10^{expn}"


# ============================================================
# ============================================================

def compute_region_counts_sum(
    tissue: str,
    chrom: str,
    region_start: int,
    region_end: int,
    real_base: str,
    pred_base: str,
    real_suffix: str,
    pred_suffix: str
):
    """
    For one tissue:
      observed = sum(counts) over peaks overlapping region
      predicted = sum(counts) over peaks overlapping region
    """
    real_file = os.path.join(real_base, tissue, real_suffix.lstrip("/"))
    pred_file = os.path.join(pred_base, tissue, f"{tissue}{pred_suffix}")

    if not os.path.exists(real_file) or not os.path.exists(pred_file):
        return None

    # read name (col0) + counts (col3)
    real = pd.read_csv(real_file, sep=r"\s+", header=None, usecols=[0, 3], engine="python")
    pred = pd.read_csv(pred_file, sep=r"\s+", header=None, usecols=[0, 3], engine="python")
    real.columns = ["name", "observed"]
    pred.columns = ["name", "predicted"]

    df = pd.merge(real, pred, on="name", how="inner")
    if df.empty:
        return None

    c, s, e = parse_peak_coords(df["name"])
    df["chr"] = c
    df["start"] = s
    df["end"] = e
    df = df.dropna(subset=["start", "end"])
    if df.empty:
        return None

    # overlap filter
    df = df[
        (df["chr"] == chrom) &
        (df["end"] >= region_start) &
        (df["start"] <= region_end)
    ]
    if df.empty:
        return None

    return {
        "Tissue": tissue,
        "observed": float(df["observed"].sum()),
        "predicted": float(df["predicted"].sum()),
        "n_peaks": int(df.shape[0])
    }


# ============================================================
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description="Counts+SUM per region, one point per tissue, label top-right N tissues."
    )
    ap.add_argument("--region", required=True, help="chr:start-end")
    ap.add_argument("--feature", default="", help="Optional gene/feature name for display")

    ap.add_argument("--label-top-n", type=int, default=1,
                    help="Label top N tissues with highest signal (top-right points).")

    ap.add_argument("--tissue-file", default="/data/home/sczd644/run/zsw_chrombpnet/tissue.txt")
    ap.add_argument("--tissue-colors",
                    default="/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv")

    ap.add_argument("--real-base", default="/data/home/sczd644/run/zsw_chrombpnet/ATAC_bams")
    ap.add_argument("--pred-base", default="/data/home/sczd644/run/zsw_chrombpnet/pred_bw")
    ap.add_argument("--real-suffix", default="/data/peaks_no_blacklist.observed.out")
    ap.add_argument("--pred-suffix", default="_peaks_no_blacklist.predicted.out")

    ap.add_argument("--pseudocount", type=float, default=1.0)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--point-size", type=float, default=60.0)

    args = ap.parse_args()

    chrom, rstart, rend, region_main = parse_region(args.region)
    tissues = read_tissues(args.tissue_file)
    color_map = read_tissue_colors(args.tissue_colors)

    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for t in tissues:
        out = compute_region_counts_sum(
            tissue=t,
            chrom=chrom,
            region_start=rstart,
            region_end=rend,
            real_base=args.real_base,
            pred_base=args.pred_base,
            real_suffix=args.real_suffix,
            pred_suffix=args.pred_suffix
        )
        if out is not None:
            rows.append(out)

    if len(rows) < 3:
        raise RuntimeError(f"Too few tissues with data in this region (need >=3). Got: {len(rows)}")

    df = pd.DataFrame(rows)

    # log transform for plotting/correlation
    pc = args.pseudocount
    df["log_obs"] = np.log10(df["observed"] + pc)
    df["log_pred"] = np.log10(df["predicted"] + pc)

    # Pearson on log scale (match axes)
    r_val, p_val = pearsonr(df["log_obs"].to_numpy(), df["log_pred"].to_numpy())

    a, b = np.polyfit(df["log_obs"], df["log_pred"], 1)

    df["color"] = df["Tissue"].map(color_map).fillna("grey50")

    # top-right N: score = x + y
    df["score_xy"] = df["log_obs"] + df["log_pred"]
    nlab = max(0, int(args.label_top_n))
    nlab = min(nlab, df.shape[0])
    label_df = df.sort_values("score_xy", ascending=False).head(nlab)

    region_key = re.sub(r"[:\-\s\(\)]+", "_", region_main)
    prefix = os.path.join(args.outdir, f"region_{region_key}")

    df.to_csv(prefix + ".per_tissue.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(4.5, 4.0), dpi=300)

    ax.scatter(df["log_obs"], df["log_pred"],
               s=args.point_size,
               c=df["color"].tolist(),
               marker="o")

    xx = np.linspace(df["log_obs"].min(), df["log_obs"].max(), 200)
    ax.plot(xx, a * xx + b, color="red", lw=2)

    ax.set_xlabel(f"log10(observed counts + {pc:g})")
    ax.set_ylabel(f"log10(predicted counts + {pc:g})")

    title_left = region_main if not args.feature else f"{region_main} ({args.feature})"
    ax.text(df["log_obs"].min(), df["log_pred"].max(),
            f"{title_left}\nPearson r = {r_val:.3f}\nP = {format_p_x10(p_val)}",
            ha="left", va="top", fontsize=10)

    for _, row in label_df.iterrows():
        ax.annotate(
            row["Tissue"],
            xy=(row["log_obs"], row["log_pred"]),
            xytext=(row["log_obs"] + 0.03, row["log_pred"] + 0.03),
            arrowprops=dict(arrowstyle="-", lw=1),
            fontsize=10
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out_pdf = prefix + ".obs_vs_pred.pdf"
    fig.savefig(out_pdf, format="pdf")  #  save  PDF
    plt.close(fig)

    print("DONE")
    print("Plot :", out_pdf)
    print("Data :", prefix + ".per_tissue.tsv")


if __name__ == "__main__":
    main()
