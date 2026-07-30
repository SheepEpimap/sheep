#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import math
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def mm_to_in(mm: float) -> float:
    return mm / 25.4


def read_tissue_colors(colors_tsv: str) -> dict:
    """
    tissue_colors.tsv  :tissue  color
      Tab/
    """
    cdf = pd.read_csv(colors_tsv, sep=r"\s+", engine="python")
    cdf.columns = [c.lower() for c in cdf.columns]
    if "tissue" not in cdf.columns or "color" not in cdf.columns:
        raise ValueError("colors file must contain columns: tissue, color")
    return dict(zip(cdf["tissue"].astype(str), cdf["color"].astype(str)))


def set_nature_style(font: str = "Arial", base_size: float = 7.0, line_w: float = 0.8):
    """
      Nature  : , , ,TrueType
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font, "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": base_size,
        "axes.linewidth": line_w,
        "xtick.major.width": line_w,
        "ytick.major.width": line_w,
        "xtick.minor.width": line_w,
        "ytick.minor.width": line_w,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "pdf.fonttype": 42,   # TrueType( )
        "ps.fonttype": 42,
    })


def draw_one_metric(ax, metric_name: str, values: np.ndarray,
                    colors: np.ndarray, seed: int = 1, point_size: float = 10.0):
    """
     :
      -  ( )
      -  ( )
    """
    rng = np.random.default_rng(seed)

    ax.boxplot(
        [values],
        positions=[1],
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=0.9),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
        boxprops=dict(edgecolor="black", linewidth=0.8),
    )

    for b in ax.artists:
        b.set_facecolor("white")

    jitter = rng.uniform(-0.18, 0.18, size=len(values))
    x = 1.0 + jitter
    ax.scatter(x, values, s=point_size, c=colors, edgecolors="none", alpha=0.95)

    ax.set_xlim(0.5, 1.5)
    ax.set_xticks([])
    ax.set_title(metric_name, pad=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    ap = argparse.ArgumentParser(
        description="ChromBPNet QC metrics: one boxplot per metric, overlay tissue-colored points."
    )
    ap.add_argument("--metrics", required=True, help="chrombpnet_metrics_summary.tsv")
    ap.add_argument("--colors", required=True, help="tissue_colors.tsv (tissue, color)")
    ap.add_argument("--outpdf", required=True, help="output PDF path")

    ap.add_argument("--exclude", default="pearsonr",
                    help="Comma-separated metric names to exclude (case-insensitive). Default: pearsonr")

    ap.add_argument("--ncols", type=int, default=3)
    ap.add_argument("--width-mm", type=float, default=183.0)   # Nature   183mm
    ap.add_argument("--height-mm", type=float, default=120.0)
    ap.add_argument("--font", type=str, default="Arial")
    ap.add_argument("--font-size", type=float, default=7.0)
    ap.add_argument("--line-width", type=float, default=0.8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--point-size", type=float, default=10.0)
    args = ap.parse_args()

    set_nature_style(font=args.font, base_size=args.font_size, line_w=args.line_width)

    df = pd.read_csv(args.metrics, sep=r"\s+", engine="python")
    if "tissue" not in df.columns:
        raise ValueError("metrics file must contain column: tissue")
    df["tissue"] = df["tissue"].astype(str)

    exclude_set = {x.strip().lower() for x in args.exclude.split(",") if x.strip()}
    exclude_set.add("tissue")  #   tissue

    metric_cols = [c for c in df.columns if c.lower() not in exclude_set]
    if len(metric_cols) == 0:
        raise ValueError(f"No metric columns left after excluding: {sorted(exclude_set)}")

    for c in metric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    color_map = read_tissue_colors(args.colors)
    df["color"] = df["tissue"].map(color_map).fillna("grey50")

    n_metrics = len(metric_cols)
    ncols = max(1, args.ncols)
    nrows = int(math.ceil(n_metrics / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=(mm_to_in(args.width_mm), mm_to_in(args.height_mm)),
        dpi=300
    )
    axes = np.array(axes).reshape(-1)

    for i, m in enumerate(metric_cols):
        ax = axes[i]
        sub = df[["tissue", "color", m]].dropna()
        if sub.shape[0] < 3:
            ax.text(0.5, 0.5, f"{m}\n<3 points", ha="center", va="center")
            ax.set_axis_off()
            continue

        values = sub[m].to_numpy()
        colors = sub["color"].to_numpy()

        draw_one_metric(
            ax=ax,
            metric_name=m,
            values=values,
            colors=colors,
            seed=args.seed,
            point_size=args.point_size
        )

    for j in range(n_metrics, len(axes)):
        axes[j].set_axis_off()

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.outpdf), exist_ok=True)
    fig.savefig(args.outpdf)
    plt.close(fig)

    print("DONE:", args.outpdf)
    print("Plotted metrics:", ", ".join(metric_cols))
    print("Excluded metrics:", ", ".join(sorted(exclude_set - {'tissue'})))


if __name__ == "__main__":
    main()
