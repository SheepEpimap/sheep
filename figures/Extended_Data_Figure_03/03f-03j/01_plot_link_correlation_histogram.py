#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

"""
plot_E5_sig_r_hist.py

 :
1) read <in_dir>/<prefix>_output_E5.tsv
2)   qval <= 0.05
3)   pearson_r  ( )
4)   Nature  :
   - r <= 0.3 :
   - r > 0.3  : #BFC67D
   -   x = 0.3
5)  :
   - r <= 0.3 (n=xxx)
   - r > 0.3  (n=xxx)
   - All pairs shown: qval <= 0.05

inputfile :
  pearson_r
  qval

 :
  n_samples
"""

import os
import argparse
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def iter_chunks(path: str, chunksize: int):
    """
     read TSV, read , .
    """
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")

    usecols = ["pearson_r", "qval"]
    if "n_samples" in header:
        usecols.append("n_samples")

    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=chunksize):
        yield chunk


def safe_numeric(s: pd.Series) -> pd.Series:
    """
     ,  NaN.
    """
    return pd.to_numeric(s, errors="coerce")


def accumulate_sig_hist_counts(
    file_path: str,
    bins: int,
    xmin: float,
    xmax: float,
    qmax: float,
    require_n86: bool,
    chunksize: int,
):
    """
      qval <= qmax   pearson_r  .

     :
      counts      :   bin
      edges       : bin
      total_n     :   pair
      n_le_thr    : r <= 0.3   pair
      n_gt_thr    : r > 0.3   pair
      median_r    :   r
      mean_r      :   r
    """
    edges = np.linspace(xmin, xmax, bins + 1)
    counts = np.zeros(bins, dtype=np.int64)

    all_r = []

    total_n = 0
    total_sum = 0.0
    n_le_thr = 0
    n_gt_thr = 0
    thr = 0.3

    for chunk in iter_chunks(file_path, chunksize):
        chunk["pearson_r"] = safe_numeric(chunk["pearson_r"])
        chunk["qval"] = safe_numeric(chunk["qval"])

        if "n_samples" in chunk.columns:
            chunk["n_samples"] = safe_numeric(chunk["n_samples"])

        m = np.isfinite(chunk["pearson_r"].to_numpy()) & np.isfinite(chunk["qval"].to_numpy())
        if not np.any(m):
            continue

        r = chunk.loc[m, "pearson_r"].to_numpy(dtype=np.float64)
        q = chunk.loc[m, "qval"].to_numpy(dtype=np.float64)

        if require_n86 and "n_samples" in chunk.columns:
            n = chunk.loc[m, "n_samples"].to_numpy(dtype=np.float64)
            keep_n = (n == 86)
            r = r[keep_n]
            q = q[keep_n]

        if r.size == 0:
            continue

        sig = (q <= qmax)
        if not np.any(sig):
            continue

        r = r[sig]

        keep_range = (r >= xmin) & (r <= xmax)
        r = r[keep_range]
        if r.size == 0:
            continue

        n_le_thr += int(np.sum(r <= thr))
        n_gt_thr += int(np.sum(r > thr))

        c, _ = np.histogram(r, bins=edges)
        counts += c

        total_n += r.size
        total_sum += np.sum(r)
        all_r.append(r)

    if total_n == 0:
        raise ValueError("No significant Pearson r values (qval <= qmax) found for plotting.")

    all_r = np.concatenate(all_r)
    median_r = float(np.median(all_r))
    mean_r = float(total_sum / total_n)

    return counts, edges, total_n, n_le_thr, n_gt_thr, median_r, mean_r


def plot_hist_nature(
    counts: np.ndarray,
    edges: np.ndarray,
    total_n: int,
    n_le_thr: int,
    n_gt_thr: int,
    median_r: float,
    mean_r: float,
    threshold: float,
    qmax: float,
    out_pdf: str,
):
    """
      Nature  .
    - r <= threshold :
    - r > threshold  : #BFC67D
    """
    centers = (edges[:-1] + edges[1:]) / 2.0
    widths = np.diff(edges)

    density = counts / (total_n * widths)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(3.1, 2.35))
    ax = fig.add_subplot(111)

    bar_colors = np.where(centers <= threshold, "#9A9A9A", "#BFC67D")

    ax.bar(
        centers,
        density,
        width=widths * 0.95,
        color=bar_colors,
        edgecolor="none",
        align="center"
    )

    ax.axvline(
        threshold,
        color="black",
        linestyle=(0, (4, 4)),
        linewidth=0.9
    )

    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel("Pearson r")
    ax.set_ylabel("Density")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", which="both", direction="out", length=3, width=0.8)
    ax.set_xticks([-1, -0.5, 0, 0.3, 0.5, 1.0])

    legend_handles = [
        Patch(facecolor="#BFC67D", edgecolor="none", label=f"r > 0.3 (n={n_gt_thr:,})"),
        Patch(facecolor="#9A9A9A", edgecolor="none", label=f"r ≤ 0.3 (n={n_le_thr:,})"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=False,
        handlelength=1.2,
        handleheight=1.2,
        borderpad=0.2,
        labelspacing=0.4,
        handletextpad=0.4,
    )

    ax.text(
        0.98, 0.70,
        f"All pairs shown:\nqval ≤ {qmax}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7
    )

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", required=True, help="directory containing <prefix>_output_E5.tsv")
    ap.add_argument("--prefix", required=True, help="e.g. H3K27ac")
    ap.add_argument("--out_dir", default="", help="output directory (default: in_dir)")
    ap.add_argument("--qmax", type=float, default=0.05, help="only keep qval <= qmax")
    ap.add_argument("--bins", type=int, default=80, help="number of histogram bins")
    ap.add_argument("--xmin", type=float, default=-1.0, help="minimum x for histogram")
    ap.add_argument("--xmax", type=float, default=1.0, help="maximum x for histogram")
    ap.add_argument("--threshold", type=float, default=0.3, help="color threshold")
    ap.add_argument("--chunksize", type=int, default=1000000, help="rows per chunk")
    ap.add_argument("--require_n86", action="store_true", help="keep only rows with n_samples == 86 if available")
    args = ap.parse_args()

    out_dir = args.out_dir if args.out_dir else args.in_dir
    os.makedirs(out_dir, exist_ok=True)

    file_path = os.path.join(args.in_dir, f"{args.prefix}_output_E5.tsv")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    counts, edges, total_n, n_le_thr, n_gt_thr, median_r, mean_r = accumulate_sig_hist_counts(
        file_path=file_path,
        bins=args.bins,
        xmin=args.xmin,
        xmax=args.xmax,
        qmax=args.qmax,
        require_n86=args.require_n86,
        chunksize=args.chunksize,
    )

    out_pdf = os.path.join(out_dir, f"{args.prefix}.E5.q{args.qmax}.sig_r_hist.pdf")

    plot_hist_nature(
        counts=counts,
        edges=edges,
        total_n=total_n,
        n_le_thr=n_le_thr,
        n_gt_thr=n_gt_thr,
        median_r=median_r,
        mean_r=mean_r,
        threshold=args.threshold,
        qmax=args.qmax,
        out_pdf=out_pdf,
    )

    print(f"[INFO] input file         : {file_path}")
    print(f"[INFO] qval cutoff        : <= {args.qmax}")
    print(f"[INFO] significant pairs  : {total_n:,}")
    print(f"[INFO] r <= 0.3 pairs     : {n_le_thr:,}")
    print(f"[INFO] r > 0.3 pairs      : {n_gt_thr:,}")
    print(f"[INFO] median_r           : {median_r:.6f}")
    print(f"[INFO] mean_r             : {mean_r:.6f}")
    print(f"[OK] wrote: {out_pdf}")


if __name__ == "__main__":
    main()
