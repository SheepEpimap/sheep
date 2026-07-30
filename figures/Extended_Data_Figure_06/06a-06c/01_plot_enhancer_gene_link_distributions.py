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
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["axes.linewidth"] = 1.2
mpl.rcParams["xtick.major.width"] = 1.1
mpl.rcParams["ytick.major.width"] = 1.1
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"

import matplotlib.pyplot as plt


# ============================================================
# ============================================================

INFILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_linked_enhancer_gene_pairs.tsv"
)

OUTDIR = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "zutu/final_three_cdf_plots_EnhA_combined"
)

os.makedirs(OUTDIR, exist_ok=True)

OUTPDF = os.path.join(
    OUTDIR,
    "combined_genes_per_enhancer_enhancers_per_gene_distance.cdf.pdf"
)


# ============================================================
# ============================================================

GROUP_ORDER = [
    "G1_1_tissue",
    "G2_2_5_tissues",
    "G3_6_10_tissues",
    "G4_11_20_tissues",
    "G5_21_43_tissues",
]

GROUP_LABELS = {
    "G1_1_tissue": "TS-EnhA",
    "G2_2_5_tissues": "NS-EnhA",
    "G3_6_10_tissues": "MS-EnhA",
    "G4_11_20_tissues": "BS-EnhA",
    "G5_21_43_tissues": "ES-EnhA",
}

GROUP_FULL_NAMES = {
    "G1_1_tissue": "Tissue-specific enhancer activity",
    "G2_2_5_tissues": "Narrowly shared enhancer activity",
    "G3_6_10_tissues": "Moderately shared enhancer activity",
    "G4_11_20_tissues": "Broadly shared enhancer activity",
    "G5_21_43_tissues": "Extensively shared enhancer activity",
}

GROUP_RANGES = {
    "G1_1_tissue": "1 tissue",
    "G2_2_5_tissues": "2–5 tissues",
    "G3_6_10_tissues": "6–10 tissues",
    "G4_11_20_tissues": "11–20 tissues",
    "G5_21_43_tissues": "21–43 tissues",
}

GROUP_COLORS = {
    "G1_1_tissue": "#0072B2",
    "G2_2_5_tissues": "#E69F00",
    "G3_6_10_tissues": "#009E73",
    "G4_11_20_tissues": "#D55E00",
    "G5_21_43_tissues": "#CC79A7",
}


# ============================================================
# 3. helper functions
# ============================================================

def ecdf(values):
    """
      ECDF.

    input:
        values:

    output:
        x:
        y:   x
    """
    values = pd.Series(values).dropna().astype(float).values

    if len(values) == 0:
        return np.array([]), np.array([])

    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)

    return x, y


def format_mean(value, metric):
    """
      legend  .
    """
    if pd.isna(value):
        return "NA"

    if metric == "distance_bp":
        return f"{value:.0f}"

    return f"{value:.1f}"


def nature_style(ax):
    """
     .
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=10,
        length=4,
        width=1.1
    )

    ax.grid(False)


def plot_cdf_on_ax(
    ax,
    data_dict,
    metric_name,
    xlabel,
    panel_label,
    log_x=False
):
    """
      ax   CDF  .

     :
        ax: matplotlib
        data_dict:
        metric_name:  ,  mean
        xlabel: x
        panel_label: panel  ,  a/b/c
        log_x:   log10 x
    """
    for group in GROUP_ORDER:
        values = data_dict[group]
        x, y = ecdf(values)

        if len(x) == 0:
            continue

        mean_value = pd.Series(values).dropna().astype(float).mean()
        mean_text = format_mean(mean_value, metric_name)

        ax.step(
            x,
            y,
            where="post",
            linewidth=2.0,
            color=GROUP_COLORS[group],
            label=f"{GROUP_LABELS[group]}, mean={mean_text}",
        )

    nature_style(ax)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Cumulative fraction", fontsize=12)
    ax.set_ylim(0, 1.02)

    if log_x:
        ax.set_xscale("log")

    ax.legend(
        frameon=False,
        fontsize=8.2,
        loc="lower right",
        handlelength=2.2
    )

    ax.text(
        -0.15,
        1.06,
        panel_label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left"
    )


# ============================================================
# 4. read enhancer-gene pair file
# ============================================================

df = pd.read_csv(INFILE, sep="\t", header=0, dtype=str)

required_cols = [
    "tissue_count_group",
    "tissue_count_range",
    "n_tissues",
    "enhancer_id",
    "gene",
    "distance",
]

for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"inputfileMissing required columns: {col}")

df["distance"] = pd.to_numeric(df["distance"], errors="coerce")
df["gene"] = df["gene"].replace(["NA", "nan", "None", ""], np.nan)

df = df[df["tissue_count_group"].isin(GROUP_ORDER)].copy()

df = df.dropna(subset=["gene"]).copy()

df = df.drop_duplicates(
    subset=["tissue_count_group", "enhancer_id", "gene"],
    keep="first"
).copy()

df["distance_bp"] = df["distance"].abs()


# ============================================================
# ============================================================

genes_per_enhancer_dict = {}
enhancers_per_gene_dict = {}
distance_dict = {}

for group in GROUP_ORDER:

    sub = df[df["tissue_count_group"] == group].copy()

    # Genes per enhancer:
    genes_per_enhancer = (
        sub.groupby("enhancer_id")["gene"]
        .nunique()
        .reset_index(name="genes_per_enhancer")
    )

    # Enhancers per gene:
    enhancers_per_gene = (
        sub.groupby("gene")["enhancer_id"]
        .nunique()
        .reset_index(name="enhancers_per_gene")
    )

    # Distance:
    distance_values = sub["distance_bp"].dropna().astype(float)

    genes_per_enhancer_values = genes_per_enhancer["genes_per_enhancer"].astype(float)
    enhancers_per_gene_values = enhancers_per_gene["enhancers_per_gene"].astype(float)

    genes_per_enhancer_dict[group] = genes_per_enhancer_values
    enhancers_per_gene_dict[group] = enhancers_per_gene_values
    distance_dict[group] = distance_values


# ============================================================
# ============================================================

fig, axes = plt.subplots(
    1,
    3,
    figsize=(14.8, 4.2)
)

plot_cdf_on_ax(
    ax=axes[0],
    data_dict=genes_per_enhancer_dict,
    metric_name="genes_per_enhancer",
    xlabel="Genes per enhancer",
    panel_label="a",
    log_x=False
)

plot_cdf_on_ax(
    ax=axes[1],
    data_dict=enhancers_per_gene_dict,
    metric_name="enhancers_per_gene",
    xlabel="Enhancers per gene",
    panel_label="b",
    log_x=False
)

plot_cdf_on_ax(
    ax=axes[2],
    data_dict=distance_dict,
    metric_name="distance_bp",
    xlabel="Distance (bp)",
    panel_label="c",
    log_x=True
)

plt.tight_layout(w_pad=2.2)

plt.savefig(OUTPDF, bbox_inches="tight")
plt.close()


# ============================================================
# ============================================================

print(" . output1 Combined figure.")
print(f"Output directory: {OUTDIR}")
print(f"Combined figure: {OUTPDF}")
