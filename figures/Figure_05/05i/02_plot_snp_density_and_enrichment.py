#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
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
mpl.rcParams["axes.linewidth"] = 1.1
mpl.rcParams["xtick.major.width"] = 1.0
mpl.rcParams["ytick.major.width"] = 1.0
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.colors import TwoSlopeNorm
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

from scipy.stats import fisher_exact


# ============================================================
# 1. inputoutput
# ============================================================

BASE_DIR = "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/enhancer_5group_snp_overlap_python"

SUMMARY_FILE = os.path.join(
    BASE_DIR,
    "05_group_summary.no_enhancer_merge.tsv"
)

OUTDIR = os.path.join(
    BASE_DIR,
    "zutu",
    "snp_density_AB_combined_EnhA_greater_blue"
)

FIGDIR = os.path.join(OUTDIR, "figures")
STATDIR = os.path.join(OUTDIR, "statistics")

for d in [OUTDIR, FIGDIR, STATDIR]:
    os.makedirs(d, exist_ok=True)


# ============================================================
# 2. grouping
# ============================================================

GROUP_ORDER = [
    "G1_1_tissue",
    "G2_2_5_tissues",
    "G3_6_10_tissues",
    "G4_11_20_tissues",
    "G5_21_43_tissues",
]

GROUP_LABEL = {
    "G1_1_tissue": "TS-EnhA",
    "G2_2_5_tissues": "NS-EnhA",
    "G3_6_10_tissues": "MS-EnhA",
    "G4_11_20_tissues": "BS-EnhA",
    "G5_21_43_tissues": "ES-EnhA",
}

LINE_COLOR = "#0072B2"

LEFT_YLABEL_X = -0.13


# ============================================================
# ============================================================

def bh_adjust(pvalues):
    """
    Benjamini-Hochberg FDR correction.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    out = np.full(len(pvalues), np.nan)

    valid = ~np.isnan(pvalues)

    if valid.sum() == 0:
        return out

    p = pvalues[valid]
    n = len(p)

    order = np.argsort(p)
    ranked = p[order]

    adj = np.empty(n, dtype=float)
    prev = 1.0

    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[i] = prev

    tmp = np.empty(n, dtype=float)
    tmp[order] = np.minimum(adj, 1.0)

    out[valid] = tmp

    return out


def q_to_star(q):
    """
    q value  .
    """
    if pd.isna(q):
        return ""

    if q < 0.001:
        return "***"

    if q < 0.01:
        return "**"

    if q < 0.05:
        return "*"

    return ""


def safe_log2_or(a, b, c, d):
    """
      0.5 Haldane-Anscombe correction,  0   OR  .
    """
    return np.log2(
        ((a + 0.5) / (b + 0.5)) /
        ((c + 0.5) / (d + 0.5))
    )


def nature_style_ax(ax):
    """
     .
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        labelsize=10,
        length=4,
        width=1.0
    )

    ax.grid(False)


# ============================================================
# 4. read summary
# ============================================================

def read_summary():
    """
    read SNP overlap summary.
    """
    if not os.path.exists(SUMMARY_FILE):
        raise FileNotFoundError(
            "Summary file not found:{}\nRun the SNP-overlap script first.".format(SUMMARY_FILE)
        )

    df = pd.read_csv(SUMMARY_FILE, sep="\t")

    required = [
        "group",
        "group_label",
        "enhancer_n",
        "total_enhancer_bp_no_merge",
        "overlap_unique_snp_n",
        "overlap_unique_snp_per_mb",
        "enhancers_with_snp_n",
        "enhancers_with_snp_pct",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError("summary fileMissing required columns: {}".format(missing))

    df = df[df["group"].isin(GROUP_ORDER)].copy()

    for col in [
        "enhancer_n",
        "total_enhancer_bp_no_merge",
        "overlap_unique_snp_n",
        "overlap_unique_snp_per_mb",
        "enhancers_with_snp_n",
        "enhancers_with_snp_pct",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["group"] = pd.Categorical(
        df["group"],
        categories=GROUP_ORDER,
        ordered=True
    )

    df = df.sort_values("group").copy()
    df["group"] = df["group"].astype(str)

    df["group_label"] = df["group"].map(GROUP_LABEL)

    return df


# ============================================================
# ============================================================

def fisher_snp_density_enrichment(summary):
    """
      enhancer bp   SNP  .

    2x2:
                    SNP positions    non-SNP enhancer bp
    current group        a                    b
    other groups         c                    d

     :
        alternative="greater"

     :
          SNP density odds  .
    """
    rows = []

    for group in GROUP_ORDER:
        this = summary[summary["group"] == group].copy()
        rest = summary[summary["group"] != group].copy()

        a = int(this["overlap_unique_snp_n"].sum())
        bp_this = int(this["total_enhancer_bp_no_merge"].sum())
        b = bp_this - a

        c = int(rest["overlap_unique_snp_n"].sum())
        bp_rest = int(rest["total_enhancer_bp_no_merge"].sum())
        d = bp_rest - c

        if b < 0 or d < 0:
            raise ValueError(
                "  SNP   enhancer bp,  SNP density Fisher  ."
            )

        oddsratio, p = fisher_exact(
            [[a, b], [c, d]],
            alternative="greater"
        )

        log2_or = safe_log2_or(a, b, c, d)

        density_this = a / bp_this * 1_000_000 if bp_this > 0 else np.nan
        density_rest = c / bp_rest * 1_000_000 if bp_rest > 0 else np.nan

        rows.append({
            "test_type": "SNP_density_per_bp_enrichment_greater",
            "group": group,
            "group_label": GROUP_LABEL[group],
            "a_current_snp": a,
            "b_current_non_snp_bp": b,
            "c_other_snp": c,
            "d_other_non_snp_bp": d,
            "current_enhancer_bp": bp_this,
            "other_enhancer_bp": bp_rest,
            "current_snp_per_mb": density_this,
            "other_snp_per_mb": density_rest,
            "fisher_or": oddsratio,
            "log2_or": log2_or,
            "fisher_p_enrichment": p,
        })

    out = pd.DataFrame(rows)

    out["fisher_q_BH_enrichment"] = bh_adjust(out["fisher_p_enrichment"].values)
    out["significance_enrichment"] = out["fisher_q_BH_enrichment"].apply(q_to_star)

    return out


# ============================================================
# ============================================================

def draw_vector_log2or_heatmap(ax_heat, mat, qmat, group_labels, x, norm, cmap, vmax):
    """
      Rectangle   heatmap.

     :
    1.   imshow output  image/raster  ;
    2.   PDF   Illustrator   heatmap  ;
    3.   heatmap   Illustrator  .
    """
    mat_values = mat.values.astype(float)

    ax_heat.set_xlim(0.5, len(group_labels) + 0.5)
    ax_heat.set_ylim(0.5, -0.5)

    for idx, xpos in enumerate(x):
        value = mat_values[0, idx]

        if pd.isna(value) or not np.isfinite(value):
            face_color = "#FFFFFF"
        else:
            face_color = mcolors.to_hex(
                cmap(norm(value)),
                keep_alpha=False
            )

        rect = Rectangle(
            (xpos - 0.5, -0.5),
            1.0,
            1.0,
            facecolor=face_color,
            edgecolor="none",
            linewidth=0,
            antialiased=False
        )

        ax_heat.add_patch(rect)

        q = qmat.iloc[0, idx]

        if pd.isna(value) or pd.isna(q):
            continue

        if value <= 0:
            continue

        star = q_to_star(q)

        if star == "":
            continue

        if abs(value) >= 0.55 * vmax:
            text_color = "white"
            stroke_color = "black"
        else:
            text_color = "black"
            stroke_color = "white"

        ax_heat.text(
            xpos,
            0,
            star,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=text_color,
            path_effects=[
                pe.withStroke(
                    linewidth=1.3,
                    foreground=stroke_color
                )
            ],
        )

    ax_heat.set_xticks(x)
    ax_heat.set_xticklabels(
        group_labels,
        rotation=35,
        ha="right",
        fontsize=10
    )

    ax_heat.tick_params(
        axis="x",
        bottom=False,
        labelbottom=True,
        top=False,
        labeltop=False,
        length=0,
        width=0,
        pad=6
    )

    ax_heat.set_xlabel(
        "EnhA tissue-breadth group",
        fontsize=10,
        labelpad=7
    )

    ax_heat.set_yticks([0])
    ax_heat.set_yticklabels(
        ["SNP density"],
        fontsize=9
    )

    ax_heat.tick_params(
        axis="y",
        left=False,
        right=False,
        labelleft=True,
        pad=7
    )

    ax_heat.set_ylabel("")

    for spine in ax_heat.spines.values():
        spine.set_visible(False)


def draw_vector_horizontal_colorbar(ax_cbar, cmap, norm, vmin, vmax, n_steps=256):
    """
      Rectangle   colorbar.

    n_steps  ,colorbar  ;
    n_steps  ,Illustrator  .
      256, .
    """
    edges = np.linspace(vmin, vmax, n_steps + 1)

    for i in range(n_steps):
        left = edges[i]
        right = edges[i + 1]
        mid = (left + right) / 2.0

        face_color = mcolors.to_hex(
            cmap(norm(mid)),
            keep_alpha=False
        )

        rect = Rectangle(
            (left, 0),
            right - left,
            1,
            facecolor=face_color,
            edgecolor="none",
            linewidth=0,
            antialiased=False
        )

        ax_cbar.add_patch(rect)

    ax_cbar.set_xlim(vmin, vmax)
    ax_cbar.set_ylim(0, 1)

    ax_cbar.set_yticks([])
    ax_cbar.xaxis.set_major_locator(MaxNLocator(nbins=5))

    ax_cbar.set_xlabel(
        "log2 OR",
        fontsize=9,
        labelpad=4
    )

    ax_cbar.tick_params(
        axis="x",
        labelsize=8,
        length=3,
        width=0.8
    )

    for spine in ax_cbar.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)


# ============================================================
# ============================================================

def plot_combined_AB_snp_density(summary, density_enrich):
    """
      a/b  :

     :
        Observed SNP density

     :
        SNP density enrichment heatmap

     :
        1.   heatmap   x  ,x = 1,2,3,4,5;
        2. heatmap   0.5-5.5;
        3.   x  ;
        4.   c/d;
        5.  ;
        6.   #0072B2;
        7. heatmap   log2OR  , ;
        8.   log2OR > 0   one-sided enrichment q < 0.05  ;
        9. grouping  heatmap  ;
        10. colorbar  ;
        11. heatmap   colorbar   Rectangle  ,  imshow/fig.colorbar.
    """
    group_labels = [GROUP_LABEL[g] for g in GROUP_ORDER]
    x = np.arange(1, len(GROUP_ORDER) + 1)

    summary_idx = summary.set_index("group")
    density_idx = density_enrich.set_index("group").loc[GROUP_ORDER].reset_index()

    fig = plt.figure(figsize=(5.8, 5.0))

    left_x0 = 0.16
    left_w = 0.78

    cbar_y0 = 0.10
    cbar_h = 0.035

    gap_heat_cbar = 0.145

    heat_y0 = cbar_y0 + cbar_h + gap_heat_cbar

    fig_w, fig_h = fig.get_size_inches()
    heat_h = (left_w * fig_w / len(GROUP_ORDER)) / fig_h

    gap_line_heat = 0.015

    line_y0 = heat_y0 + heat_h + gap_line_heat
    line_h = 0.90 - line_y0

    if line_h <= 0:
        raise ValueError(
            " .  figsize  ,  gap_heat_cbar / heat_h."
        )

    ax_line = fig.add_axes([left_x0, line_y0, left_w, line_h])
    ax_heat = fig.add_axes([left_x0, heat_y0, left_w, heat_h], sharex=ax_line)
    ax_cbar = fig.add_axes([left_x0, cbar_y0, left_w, cbar_h])

    # ========================================================
    # ========================================================

    y_density = [
        summary_idx.loc[g, "overlap_unique_snp_per_mb"]
        if g in summary_idx.index else np.nan
        for g in GROUP_ORDER
    ]

    ax_line.plot(
        x,
        y_density,
        marker="o",
        linewidth=2.2,
        markersize=5,
        color=LINE_COLOR,
    )

    ax_line.set_xlim(0.5, len(GROUP_ORDER) + 0.5)

    ax_line.set_xticks(x)
    ax_line.tick_params(
        axis="x",
        bottom=False,
        labelbottom=False,
        top=False,
        labeltop=False,
        length=0,
        width=0
    )

    ax_line.set_xlabel("")

    ax_line.set_ylabel(
        "Unique SNPs per Mb enhancer",
        fontsize=10
    )

    ax_line.yaxis.set_label_coords(LEFT_YLABEL_X, 0.5)

    nature_style_ax(ax_line)


    # ========================================================
    # ========================================================

    mat = pd.DataFrame(
        index=["SNP density"],
        columns=group_labels,
        dtype=float
    )

    qmat = pd.DataFrame(
        index=["SNP density"],
        columns=group_labels,
        dtype=float
    )

    for _, row in density_idx.iterrows():
        group_label = row["group_label"]
        mat.loc["SNP density", group_label] = row["log2_or"]
        qmat.loc["SNP density", group_label] = row["fisher_q_BH_enrichment"]

    mat_values = mat.values.astype(float)

    vmax = np.nanmax(np.abs(mat_values))

    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1

    norm = TwoSlopeNorm(
        vmin=-vmax,
        vcenter=0,
        vmax=vmax
    )

    cmap = plt.get_cmap("RdBu_r")

    draw_vector_log2or_heatmap(
        ax_heat=ax_heat,
        mat=mat,
        qmat=qmat,
        group_labels=group_labels,
        x=x,
        norm=norm,
        cmap=cmap,
        vmax=vmax
    )

    # ========================================================
    # ========================================================

    draw_vector_horizontal_colorbar(
        ax_cbar=ax_cbar,
        cmap=cmap,
        norm=norm,
        vmin=-vmax,
        vmax=vmax,
        n_steps=256
    )

    out_pdf = os.path.join(
        FIGDIR,
        "combined_AB_snp_density_and_enrichment_EnhA.greater_blue.pdf"
    )

    out_png = os.path.join(
        FIGDIR,
        "combined_AB_snp_density_and_enrichment_EnhA.greater_blue.png"
    )

    out_svg = os.path.join(
        FIGDIR,
        "combined_AB_snp_density_and_enrichment_EnhA.greater_blue.svg"
    )

    plt.savefig(
        out_pdf,
        bbox_inches="tight",
        dpi=600,
        facecolor="white",
        transparent=False
    )

    plt.savefig(
        out_png,
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        transparent=False
    )

    plt.savefig(
        out_svg,
        bbox_inches="tight",
        facecolor="white",
        transparent=False
    )

    plt.close()

    return out_pdf, out_png, out_svg


# ============================================================
# ============================================================

def main():
    print("[INFO] Reading summary:")
    print(SUMMARY_FILE)

    summary = read_summary()

    summary_out = os.path.join(
        STATDIR,
        "input_group_summary.used.tsv"
    )

    summary.to_csv(
        summary_out,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    print("[INFO] Running one-sided Fisher enrichment test for SNP density...")
    density_enrich = fisher_snp_density_enrichment(summary)

    density_file = os.path.join(
        STATDIR,
        "snp_density_enrichment.current_group_vs_others.fisher_greater.tsv"
    )

    density_enrich.to_csv(
        density_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    print("[INFO] Plotting combined A/B SNP density figure only...")
    out_pdf, out_png, out_svg = plot_combined_AB_snp_density(
        summary=summary,
        density_enrich=density_enrich
    )

    print("")
    print(" . output a/b  , output c/d.")
    print("Output directory:", OUTDIR)
    print("Input summary:", summary_out)
    print("SNP density enrichment greater:", density_file)
    print("Figure PDF:", out_pdf)
    print("Figure PNG:", out_png)
    print("Figure SVG:", out_svg)


if __name__ == "__main__":
    main()
