#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

import matplotlib.pyplot as plt


def read_table_flexible(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", engine="python")


def _normalize_mpl_color(s: str) -> str:
    if s is None:
        return "#7f7f7f"
    s = str(s).strip()
    m = re.fullmatch(r"(grey|gray)(\d{1,3})", s, flags=re.IGNORECASE)
    if m:
        v = max(0, min(100, int(m.group(2))))
        return str(v / 100.0)
    return s


def load_tissue_colors(colors_tsv: str) -> dict:
    c = read_table_flexible(colors_tsv)
    c.columns = [x.lower() for x in c.columns]
    if "tissue" not in c.columns or "color" not in c.columns:
        raise ValueError("colors file must contain columns: tissue, color")
    return {str(t): _normalize_mpl_color(col) for t, col in zip(c["tissue"], c["color"])}


def add_break_marks(ax_top, ax_bot, d=0.008, lw=0.8):
    kwargs = dict(transform=ax_top.transAxes, color="black", clip_on=False, linewidth=lw)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

    kwargs = dict(transform=ax_bot.transAxes, color="black", clip_on=False, linewidth=lw)
    ax_bot.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_bot.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)


def main():
    ap = argparse.ArgumentParser(
        description="Per-rep correlation barplot with a y-axis break keeping 0 tick; editable PDF/SVG."
    )
    ap.add_argument("--corfile", required=True)
    ap.add_argument("--colors", required=True)
    ap.add_argument("--outpdf", required=True)
    ap.add_argument("--outsvg", default="")

    ap.add_argument("--mode", choices=["tissue", "tissuepair"], default="tissuepair")

    ap.add_argument("--label-col", choices=["pair", "tissuepair"], default="pair")

    ap.add_argument("--break-low", type=float, default=0.05)   #
    ap.add_argument("--break-high", type=float, default=0.50)  #
    ap.add_argument("--bottom-base", type=float, default=0.0)  #  ,  0   0

    ap.add_argument("--sort", choices=["desc", "asc", "none"], default="desc")
    ap.add_argument("--bar-edge", action="store_true")
    ap.add_argument("--rotate", type=int, default=90)

    ap.add_argument("--gap", type=float, default=0.04)
    ap.add_argument("--top-height", type=float, default=6.0)
    ap.add_argument("--bottom-height", type=float, default=0.7)

    args = ap.parse_args()

    low, high = args.break_low, args.break_high
    if not (low < high):
        raise ValueError("--break-low must be < --break-high")
    if args.bottom_base != 0.0:
        raise ValueError("For showing 0 tick, please keep --bottom-base 0.0")

    df = read_table_flexible(args.corfile)
    df.columns = [c.lower() for c in df.columns]
    if "tissue" not in df.columns or "cor" not in df.columns:
        raise ValueError("corfile must contain columns: Tissue and cor")

    df["tissue"] = df["tissue"].astype(str)
    df["cor"] = pd.to_numeric(df["cor"], errors="coerce")
    df = df.dropna(subset=["cor"]).copy()

    if args.mode == "tissue":
        plot_df = df.groupby("tissue", as_index=False)["cor"].mean()
        tissues_for_color = plot_df["tissue"].tolist()
        corvals = plot_df["cor"].to_numpy()
        x_labels = plot_df["tissue"].astype(str).tolist()

    else:
        label_col = args.label_col.lower()
        if label_col not in df.columns:
            if label_col == "pair" and "tissuepair" in df.columns:
                label_col = "tissuepair"
            else:
                raise ValueError(f"mode=tissuepair requires column: {args.label_col}")

        plot_df = df.groupby(["tissue", label_col], as_index=False)["cor"].mean()
        plot_df = plot_df.rename(columns={label_col: "label"})

        tissues_for_color = plot_df["tissue"].astype(str).tolist()
        corvals = plot_df["cor"].to_numpy()

        if label_col == "pair":
            x_labels = (plot_df["tissue"].astype(str) + "\n" + plot_df["label"].astype(str)).tolist()
        else:
            x_labels = plot_df["label"].astype(str).tolist()

    if args.sort == "desc":
        order = np.argsort(-corvals)
    elif args.sort == "asc":
        order = np.argsort(corvals)
    else:
        order = np.arange(len(corvals))

    corvals = corvals[order]
    x_labels = [x_labels[i] for i in order]
    tissues_for_color = [tissues_for_color[i] for i in order]

    in_gap = (corvals > low) & (corvals < high)
    if np.any(in_gap):
        raise ValueError(
            f"Break interval ({low:.3f}, {high:.3f}) cuts real data. "
            f"Pick break-low smaller (e.g. 0.02~0.08) and break-high=0.50."
        )

    cmap = load_tissue_colors(args.colors)
    colors = [cmap.get(t, "#7f7f7f") for t in tissues_for_color]

    y_max = float(np.nanmax(corvals))
    bot_lo, bot_hi = args.bottom_base, low          # ✅   0
    top_lo, top_hi = high, y_max + 0.01

    n = len(corvals)
    fig_w = max(10, 0.28 * n)
    fig_h = 6.5

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True,
        figsize=(fig_w, fig_h),
        gridspec_kw={"height_ratios": [args.top_height, args.bottom_height]}
    )
    fig.subplots_adjust(hspace=args.gap)

    x = np.arange(n)
    edgecolor = "black" if args.bar_edge else None
    linewidth = 0.4 if args.bar_edge else 0.0

    ax_top.bar(x, corvals, color=colors, edgecolor=edgecolor, linewidth=linewidth)
    ax_bot.bar(x, corvals, color=colors, edgecolor=edgecolor, linewidth=linewidth)

    ax_top.set_ylim(top_lo, top_hi)
    ax_bot.set_ylim(bot_lo, bot_hi)

    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_bot.spines["right"].set_visible(False)

    ax_bot.set_yticks([0.0])
    ax_bot.set_yticklabels(["0"])

    top_ticks = np.arange(high, top_hi + 1e-9, 0.1)
    ax_top.set_yticks(np.round(top_ticks, 2))

    add_break_marks(ax_top, ax_bot, d=0.008, lw=0.8)

    ax_top.set_title("Observed vs Predicted accessibility correlation (per-rep)")
    fig.text(0.02, 0.5, "Pearson correlation (r)", va="center", rotation=90)

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(x_labels, rotation=args.rotate, ha="right", fontsize=8)

    ax_top.grid(False)
    ax_bot.grid(False)

    # ---------- save ----------
    outdir = os.path.dirname(args.outpdf)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    fig.tight_layout(rect=[0.04, 0.02, 1, 0.98])
    fig.savefig(args.outpdf, format="pdf")
    if args.outsvg:
        svgdir = os.path.dirname(args.outsvg)
        if svgdir:
            os.makedirs(svgdir, exist_ok=True)
        fig.savefig(args.outsvg, format="svg")

    plt.close(fig)

    print("DONE:")
    print("  PDF:", args.outpdf)
    if args.outsvg:
        print("  SVG:", args.outsvg)
    print(f"  y break: bottom=[{bot_lo},{bot_hi}]  top=[{top_lo},{top_hi:.3f}]")


if __name__ == "__main__":
    main()
