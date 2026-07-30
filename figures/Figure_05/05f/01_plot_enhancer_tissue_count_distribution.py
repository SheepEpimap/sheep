#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["axes.linewidth"] = 1.2
mpl.rcParams["xtick.major.width"] = 1.2
mpl.rcParams["ytick.major.width"] = 1.2
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# ============================================================
# ============================================================

INFILE = "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/all_E5_Gs_one_count.csv"

BASE_OUTDIR = "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution"
OUTDIR = os.path.join(BASE_OUTDIR, "zutu")
PLOT_DIR = os.path.join(OUTDIR, "plots")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


# ============================================================
# ============================================================

GROUP_COLORS = {
    "G1_1_tissue": "#0072B2",
    "G2_2_5_tissues": "#E69F00",
    "G3_6_10_tissues": "#009E73",
    "G4_11_20_tissues": "#D55E00",
    "G5_21_43_tissues": "#CC79A7",
}

GROUP_LABELS = {
    "G1_1_tissue": "TS-EnhA",
    "G2_2_5_tissues": "NS-EnhA",
    "G3_6_10_tissues": "MS-EnhA",
    "G4_11_20_tissues": "BS-EnhA",
    "G5_21_43_tissues": "ES-EnhA",
}

GROUP_ORDER = [
    "G1_1_tissue",
    "G2_2_5_tissues",
    "G3_6_10_tissues",
    "G4_11_20_tissues",
    "G5_21_43_tissues",
]


def assign_group(n):
    n = int(n)
    if n == 1:
        return "G1_1_tissue"
    elif 2 <= n <= 5:
        return "G2_2_5_tissues"
    elif 6 <= n <= 10:
        return "G3_6_10_tissues"
    elif 11 <= n <= 20:
        return "G4_11_20_tissues"
    elif 21 <= n <= 43:
        return "G5_21_43_tissues"
    else:
        return "Unknown"


def comma_formatter(x, pos):
    return f"{int(x):,}"


# ============================================================
# ============================================================

header = pd.read_csv(INFILE, sep="\t", nrows=0)
cols = [x.strip() for x in header.columns]

if len(cols) < 4:
    raise ValueError("inputfile , checkfile .")

count_col = cols[-1]
print(f"[INFO]  read : {count_col}")

df_count = pd.read_csv(
    INFILE,
    sep="\t",
    usecols=[count_col],
)

df_count.columns = ["n_tissues"]
df_count["n_tissues"] = pd.to_numeric(df_count["n_tissues"], errors="coerce")
df_count = df_count.dropna(subset=["n_tissues"])
df_count["n_tissues"] = df_count["n_tissues"].astype(int)

df_count = df_count[
    (df_count["n_tissues"] >= 1) &
    (df_count["n_tissues"] <= 43)
].copy()

total_enha = df_count.shape[0]
print(f"[INFO]   EnhA  : {total_enha:,}")


# ============================================================
# ============================================================

summary = (
    df_count["n_tissues"]
    .value_counts()
    .reindex(range(1, 44), fill_value=0)
    .reset_index()
)

summary.columns = ["n_tissues", "enhancer_count"]
summary["percent"] = summary["enhancer_count"] / summary["enhancer_count"].sum() * 100
summary["tissue_count_group"] = summary["n_tissues"].apply(assign_group)
summary["group_label"] = summary["tissue_count_group"].map(GROUP_LABELS)
summary["color"] = summary["tissue_count_group"].map(GROUP_COLORS)

summary_file = os.path.join(OUTDIR, "enhancer_tissue_count_distribution.summary.only_plot.tsv")
summary.to_csv(summary_file, sep="\t", index=False)


# ============================================================
# ============================================================

x = summary["n_tissues"].astype(int)
y = summary["enhancer_count"].astype(int)
colors = summary["color"].tolist()

Y_BOTTOM_MIN = 0
Y_BOTTOM_MAX = 75000
Y_TOP_MIN = 150000
Y_TOP_MAX = max(y.max() * 1.08, Y_TOP_MIN * 1.05)

fig, (ax_top, ax_bottom) = plt.subplots(
    2,
    1,
    sharex=True,
    figsize=(10.8, 2),
    gridspec_kw={
        "height_ratios": [1.0, 3.0],
        "hspace": 0.05,
    },
)

ax_top.bar(
    x,
    y,
    width=0.78,
    color=colors,
    edgecolor="black",
    linewidth=0.65,
)

ax_bottom.bar(
    x,
    y,
    width=0.78,
    color=colors,
    edgecolor="black",
    linewidth=0.65,
)

ax_bottom.set_ylim(Y_BOTTOM_MIN, Y_BOTTOM_MAX)
ax_top.set_ylim(Y_TOP_MIN, Y_TOP_MAX)

ax_top.spines["bottom"].set_visible(False)
ax_bottom.spines["top"].set_visible(False)

for ax in [ax_top, ax_bottom]:
    ax.spines["right"].set_visible(False)
    ax.yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    ax.tick_params(axis="y", which="major", labelsize=8, length=3.5, width=1.0)
    ax.yaxis.grid(True, linestyle="-", linewidth=0.35, alpha=0.16)
    ax.xaxis.grid(False)
    ax.spines["left"].set_linewidth(1.1)

ax_top.spines["top"].set_visible(False)
ax_bottom.spines["bottom"].set_linewidth(1.1)

ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
ax_bottom.tick_params(axis="x", which="major", labelsize=8, length=3.5, width=1.0)

ax_bottom.set_xticks(range(1, 44, 2))
ax_bottom.set_xlim(0.4, 43.6)

ax_bottom.set_xlabel("Number of tissues in which an EnhA is present", fontsize=10)
fig.text(
    0.025,
    0.52,
    "Number of EnhA",
    va="center",
    rotation="vertical",
    fontsize=10,
)

ax_top.set_title("Distribution of EnhA across 43 tissues", fontsize=11, pad=5)

d = 0.008
kwargs = dict(color="black", clip_on=False, linewidth=1.1)

ax_top.plot(
    (-d, +d),
    (-d, +d),
    transform=ax_top.transAxes,
    **kwargs,
)

ax_bottom.plot(
    (-d, +d),
    (1 - d, 1 + d),
    transform=ax_bottom.transAxes,
    **kwargs,
)

legend_handles = []
for g in GROUP_ORDER:
    legend_handles.append(
        mpl.patches.Patch(
            facecolor=GROUP_COLORS[g],
            edgecolor="black",
            linewidth=0.5,
            label=GROUP_LABELS[g],
        )
    )

ax_top.legend(
    handles=legend_handles,
    frameon=False,
    fontsize=8,
    ncol=5,
    loc="upper right",
    bbox_to_anchor=(1.0, 1.65),
    handlelength=1.1,
    columnspacing=0.9,
)

plt.tight_layout(rect=[0.045, 0.02, 1.0, 0.88])

pdf_file = os.path.join(PLOT_DIR, "enhancer_tissue_count_distribution.5group_colors.ybreak.half_height.no_numbers.pdf")
png_file = os.path.join(PLOT_DIR, "enhancer_tissue_count_distribution.5group_colors.ybreak.half_height.no_numbers.png")
svg_file = os.path.join(PLOT_DIR, "enhancer_tissue_count_distribution.5group_colors.ybreak.half_height.no_numbers.svg")

plt.savefig(pdf_file, bbox_inches="tight")
plt.savefig(png_file, dpi=600, bbox_inches="tight")
plt.savefig(svg_file, bbox_inches="tight")
plt.close()

print(" .")
print(f"inputfile: {INFILE}")
print(f" : {count_col}")
print(f"EnhA  : {total_enha:,}")
print(f"Summary table: {summary_file}")
print(f"PDF : {pdf_file}")
print(f"PNG : {png_file}")
print(f"SVG : {svg_file}")
