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


# =========================
# =========================

INFILE = "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/all_E5_Gs_one_count.csv"

OUTDIR = "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution"

SPLIT_DIR = os.path.join(OUTDIR, "enhancers_by_tissue_count")
PLOT_DIR = os.path.join(OUTDIR, "plots")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(SPLIT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


# =========================
# 2. readfile
# =========================

df = pd.read_csv(INFILE, sep="\t", header=0, dtype=str)
df.columns = [x.strip() for x in df.columns]

required_cols = ["chr", "start", "end"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"The input file is missing a required column: {col}")

coord_cols = ["chr", "start", "end"]
raw_count_col = df.columns[-1]
tissue_cols = list(df.columns[3:-1])

if len(tissue_cols) != 43:
    raise ValueError(
        f"Expected 43 tissue columns; detected {len(tissue_cols)}  .\n"
        f"Check the input delimiter and column count.\n"
        f"Detected tissue columns:\n{tissue_cols}"
    )


# =========================
# =========================

for col in tissue_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

df[raw_count_col] = pd.to_numeric(df[raw_count_col], errors="coerce")
df["n_tissues_calculated"] = df[tissue_cols].sum(axis=1)

mismatch = df[df[raw_count_col] != df["n_tissues_calculated"]].copy()
mismatch_file = os.path.join(OUTDIR, "check_mismatch_original_count_vs_calculated.tsv")

if mismatch.shape[0] > 0:
    mismatch.to_csv(mismatch_file, sep="\t", index=False)
    print(f"[WARNING]   {mismatch.shape[0]}  .")
    print(f"[WARNING] checkfile: {mismatch_file}")
else:
    print("[OK]  43 .")

df["n_tissues"] = df["n_tissues_calculated"].astype(int)


# =========================
# =========================

summary = (
    df["n_tissues"]
    .value_counts()
    .reindex(range(1, 44), fill_value=0)
    .reset_index()
)
summary.columns = ["n_tissues", "enhancer_count"]
summary["percent"] = summary["enhancer_count"] / summary["enhancer_count"].sum() * 100

summary_file = os.path.join(OUTDIR, "enhancer_tissue_count_distribution.summary.tsv")
summary.to_csv(summary_file, sep="\t", index=False)


# =========================
# =========================

x = summary["n_tissues"]
y = summary["enhancer_count"]

fig, ax = plt.subplots(figsize=(10.5, 4.8))

bars = ax.bar(
    x,
    y,
    width=0.78,
    color="#4D4D4D",
    edgecolor="black",
    linewidth=0.8
)

ax.set_xlabel("Number of tissues in which an enhancer is present", fontsize=12)
ax.set_ylabel("Number of enhancers", fontsize=12)
ax.set_title("Distribution of enhancers across 43 tissues", fontsize=13, pad=10)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

ax.set_xticks(range(1, 44, 2))
ax.set_xlim(0.4, 43.6)
ax.tick_params(axis="both", which="major", labelsize=10, length=4, width=1.1)

ax.yaxis.grid(True, linestyle="-", linewidth=0.4, alpha=0.18)
ax.xaxis.grid(False)

ymax = y.max()
ax.set_ylim(0, ymax * 1.14 if ymax > 0 else 1)

for bar, value in zip(bars, y):
    if value == 0:
        continue
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + ymax * 0.012,
        f"{value}",
        ha="center",
        va="bottom",
        fontsize=8,
        rotation=90
    )

plt.tight_layout()

pdf_file = os.path.join(PLOT_DIR, "enhancer_tissue_count_distribution.nature_style.pdf")
png_file = os.path.join(PLOT_DIR, "enhancer_tissue_count_distribution.nature_style.png")
svg_file = os.path.join(PLOT_DIR, "enhancer_tissue_count_distribution.nature_style.svg")

plt.savefig(pdf_file, bbox_inches="tight")
plt.savefig(png_file, dpi=600, bbox_inches="tight")
plt.savefig(svg_file, bbox_inches="tight")
plt.close()


# =========================
# =========================

out_cols = coord_cols + tissue_cols + ["n_tissues"]

for i in range(1, 44):
    sub = df[df["n_tissues"] == i].copy()
    suffix = "tissue" if i == 1 else "tissues"
    outfile = os.path.join(
        SPLIT_DIR,
        f"enhancers_present_in_{i:02d}_{suffix}.tsv"
    )
    sub[out_cols].to_csv(outfile, sep="\t", index=False)


# =========================
# =========================

print(" .")
print(f"inputfile: {INFILE}")
print(f"  enhancer  : {df.shape[0]}")
print(f" : {len(tissue_cols)}")
print(f"Summary table: {summary_file}")
print(f"PDF : {pdf_file}")
print(f"PNG : {png_file}")
print(f"SVG : {svg_file}")
print(f"43 filedirectory: {SPLIT_DIR}")
