#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
import os
import numpy as np
import pandas as pd
import pyBigWig

import matplotlib
matplotlib.use("Agg")  # HPC/ ( , )

# =========================
# =========================
matplotlib.rcParams["pdf.fonttype"] = 42   # TrueType (Type 42)
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"]  = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False

import matplotlib.pyplot as plt


window = 100          # footprint  :±100bp
tissue_file = "/data/home/sczd644/run/zsw_chrombpnet/tissue.txt"

motifs_of_interest = None

base_finemo = "/data/home/sczd644/run/zsw_chrombpnet/finemo"
base_tobias = "/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/tobias"

fig_root = os.path.join(base_tobias, "fig")
os.makedirs(fig_root, exist_ok=True)

bin_colors = {
    "1": "tab:blue",
    "2": "tab:orange",
    "3": "tab:green",
    "4": "tab:red",
}

tissues = []
with open(tissue_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tissues.append(line.split()[0])

print(" :", tissues)


def plot_one_motif(tissue, motif_name, df_all, bw_path):
    """  +   motif   2x2 footprint  """

    df = df_all[df_all["motif_name"] == motif_name].copy()
    n_raw = df.shape[0]
    print(f"[{tissue}] {motif_name}:   hits = {n_raw}")
    if n_raw == 0:
        return

    score_col = "hit_importance_sq"
    if df[score_col].nunique() < 4:
        df = df.sort_values(score_col).reset_index(drop=True)
        df["score_bin"] = pd.qcut(df.index, q=4, labels=["1", "2", "3", "4"])
    else:
        try:
            df["score_bin"] = pd.qcut(df[score_col], q=4, labels=["1", "2", "3", "4"])
        except ValueError:
            df = df.sort_values(score_col).reset_index(drop=True)
            df["score_bin"] = pd.qcut(df.index, q=4, labels=["1", "2", "3", "4"])

    print("  score_bin  :", df["score_bin"].value_counts().sort_index().to_dict())

    df["center"] = ((df["start"] + df["end"]) / 2).astype(int)

    if not os.path.exists(bw_path):
        print(f"[WARN] {tissue}: bigWig  : {bw_path},  motif")
        return

    bw = pyBigWig.open(bw_path)
    positions = np.arange(-window, window + 1)
    matrix = []

    for _, row in df.iterrows():
        chrom  = row["chr"]
        center = int(row["center"])
        strand = row["strand"]

        start = center - window
        end   = center + window + 1

        try:
            vals = np.array(bw.values(chrom, start, end))
        except RuntimeError:
            continue

        if vals is None or len(vals) != len(positions):
            continue

        vals = np.nan_to_num(vals, nan=0.0)

        if strand == "-":
            vals = vals[::-1]

        matrix.append(vals)

    bw.close()

    if len(matrix) == 0:
        print(f"[WARN] {tissue} {motif_name}:   hits  , ")
        return

    matrix = np.vstack(matrix)
    n_used = matrix.shape[0]
    print(f"    footprint   hits  : n = {n_used}")

    footprint_all = matrix.mean(axis=0)

    bin_profiles = {}
    for b in ["1", "2", "3", "4"]:
        mask = (df["score_bin"] == b).values
        if mask.sum() == 0:
            continue
        fp_bin = matrix[mask].mean(axis=0)
        bin_profiles[b] = (fp_bin, int(mask.sum()))
        print(f"  bin{b} n = {mask.sum()}")

    pdf_path = os.path.join(fig_root, f"{tissue}_{motif_name}_TOBIAS_corrected_bins.pdf")
    png_path = os.path.join(fig_root, f"{tissue}_{motif_name}_TOBIAS_corrected_bins.png")

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    # panel1: All hits
    ax = axes[0, 0]
    ax.plot(positions, footprint_all, linewidth=1.2)
    ax.axvline(0, linestyle="--", linewidth=0.8, color="gray")
    ax.set_title("All hits")
    ax.set_xlabel("Distance from motif center (bp)")
    ax.set_ylabel("Mean corrected signal")
    ax.text(0.02, 0.95, f"n = {n_used}", transform=ax.transAxes,
            va="top", ha="left", fontsize=8)

    # panel2: bin4 only
    ax = axes[0, 1]
    if "4" in bin_profiles:
        fp4, n4 = bin_profiles["4"]
        ax.plot(positions, fp4, color=bin_colors["4"], linewidth=1.2)
        ax.text(0.02, 0.95, f"bin4 (n={n4})", transform=ax.transAxes,
                va="top", ha="left", fontsize=8)
    else:
        ax.text(0.5, 0.5, "bin4  ", transform=ax.transAxes,
                ha="center", va="center")
    ax.axvline(0, linestyle="--", linewidth=0.8, color="gray")
    ax.set_title("Bin4 only")
    ax.set_xlabel("Distance from motif center (bp)")
    ax.set_ylabel("Mean corrected signal")

    # panel3: bin1 vs bin4
    ax = axes[1, 0]
    lines = []
    labels = []
    for b in ["1", "4"]:
        if b in bin_profiles:
            fpb, nb = bin_profiles[b]
            line, = ax.plot(positions, fpb, color=bin_colors[b], linewidth=1.2)
            lines.append(line)
            labels.append(f"bin{b} (n={nb})")
    if lines:
        ax.legend(lines, labels, fontsize=8, loc="best")
    else:
        ax.text(0.5, 0.5, "bin1/bin4  ", transform=ax.transAxes,
                ha="center", va="center")
    ax.axvline(0, linestyle="--", linewidth=0.8, color="gray")
    ax.set_title("Bin1 vs Bin4")
    ax.set_xlabel("Distance from motif center (bp)")
    ax.set_ylabel("Mean corrected signal")

    # panel4: bin1~bin4
    ax = axes[1, 1]
    lines = []
    labels = []
    for b in ["1", "2", "3", "4"]:
        if b in bin_profiles:
            fpb, nb = bin_profiles[b]
            line, = ax.plot(positions, fpb, color=bin_colors[b], linewidth=1.2)
            lines.append(line)
            labels.append(f"bin{b} (n={nb})")
    if lines:
        ax.legend(lines, labels, fontsize=8, loc="best")
    else:
        ax.text(0.5, 0.5, "  bin  ", transform=ax.transAxes,
                ha="center", va="center")
    ax.axvline(0, linestyle="--", linewidth=0.8, color="gray")
    ax.set_title("All quartiles (bin1–4)")
    ax.set_xlabel("Distance from motif center (bp)")
    ax.set_ylabel("Mean corrected signal")

    fig.suptitle(f"{tissue} {motif_name} TOBIAS-corrected footprints", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)

    print(f"  save : {pdf_path}")
    print(f"  save : {png_path}")


for tissue in tissues:
    print(f"\n====  : {tissue} ====")

    hits_tsv = os.path.join(base_finemo, f"{tissue}_finemo", "hits_with_tf_names.tsv")
    bw_path  = os.path.join(base_tobias, f"{tissue}_corrected.bw")

    if not os.path.exists(hits_tsv):
        print(f"[WARN] {tissue}: Not found hits file: {hits_tsv}, ")
        continue

    df_all = pd.read_csv(hits_tsv, sep="\t", header=0)

    if motifs_of_interest is None:
        motif_list = sorted(df_all["motif_name"].unique())
    else:
        motif_list = motifs_of_interest

    print(f"[{tissue}]   motifs: {motif_list}")

    for motif_name in motif_list:
        plot_one_motif(tissue, motif_name, df_all, bw_path)
