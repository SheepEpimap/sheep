#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
information(information，information):
    tissue    assay        state   value

information:
    python as_snp_chrstate_fold_enrichment_raw_plain.py
    python as_snp_chrstate_fold_enrichment_raw_plain.py /your/input/file
    python as_snp_chrstate_fold_enrichment_raw_plain.py /your/input/file /your/output_dir
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_INPUTS = [
    Path("/vol2/wulingyun/R_graph_test/sandian/260418_AS-SNP_ChrState_fold_enrichment"),
    Path(__file__).resolve().parent / "260418_AS-SNP_ChrState_fold_enrichment",
]

ASSAY_ORDER = [
    "AS_ATAC",
    "AS_H3K4me3",
    "AS_H3K4me1",
    "AS_H3K27ac",
    "AS_H3K27me3",
    "AS_RNASeq",
]

ASSAY_LABELS = {
    "AS_ATAC": "ATAC",
    "AS_H3K4me3": "H3K4me3",
    "AS_H3K4me1": "H3K4me1",
    "AS_H3K27ac": "H3K27ac",
    "AS_H3K27me3": "H3K27me3",
    "AS_RNASeq": "RNASeq",
}

ASSAY_COLORS = {
    "AS_ATAC": "#f0679f",
    "AS_H3K4me3": "#2ca25f",
    "AS_H3K4me1": "#e3b207",
    "AS_H3K27ac": "#e41a1c",
    "AS_H3K27me3": "#8c8c8c",
    "AS_RNASeq": "#6a3d9a",
}

STATE_ORDER = [f"E{i}" for i in range(1, 12)]
STATE_TO_X = {s: i for i, s in enumerate(STATE_ORDER, start=1)}

OFFSETS = {
    "AS_ATAC": -0.30,
    "AS_H3K4me3": -0.18,
    "AS_H3K4me1": -0.06,
    "AS_H3K27ac": 0.06,
    "AS_H3K27me3": 0.18,
    "AS_RNASeq": 0.30,
}

FIGSIZE = (8.8, 5.2)
DPI = 600
OUTPUT_BASENAME = "260418_AS-SNP_ChrState_fold_enrichment_raw_meanSEM_plain"


def choose_input() -> Path:
    if len(sys.argv) >= 2:
        p = Path(sys.argv[1]).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {p}")
        return p
    for p in DEFAULT_INPUTS:
        if p.exists():
            return p.resolve()
    tried = "\n  - " + "\n  - ".join(str(p) for p in DEFAULT_INPUTS)
    raise FileNotFoundError(f"Input file not found. Tried:{tried}")


def choose_output_dir(input_file: Path) -> Path:
    if len(sys.argv) >= 3:
        outdir = Path(sys.argv[2]).expanduser().resolve()
    else:
        outdir = input_file.parent
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


INPUT_FILE = choose_input()
OUTPUT_DIR = choose_output_dir(INPUT_FILE)

raw = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    header=None,
    names=["tissue", "assay", "state", "value"],
    dtype={"tissue": "string", "assay": "string", "state": "string", "value": "string"},
)

raw = raw.dropna(how="all").copy()
raw["value"] = pd.to_numeric(raw["value"], errors="coerce")
raw = raw.dropna(subset=["tissue", "assay", "state", "value"]).copy()
raw = raw[raw["assay"].isin(ASSAY_ORDER) & raw["state"].isin(STATE_ORDER)].copy()

if raw.empty:
    raise ValueError("No valid rows remained after filtering assay/state.")

summary = (
    raw.groupby(["assay", "state"], as_index=False)
       .agg(
           n=("value", "size"),
           mean_value=("value", "mean"),
           sd_value=("value", lambda x: np.std(x, ddof=1) if len(x) > 1 else 0.0),
       )
)
summary["sem_value"] = summary["sd_value"] / np.sqrt(summary["n"])
summary["x"] = summary["state"].map(STATE_TO_X).astype(float)
summary["x_plot"] = summary.apply(lambda r: r["x"] + OFFSETS[r["assay"]], axis=1)
summary["assay"] = pd.Categorical(summary["assay"], categories=ASSAY_ORDER, ordered=True)
summary["state"] = pd.Categorical(summary["state"], categories=STATE_ORDER, ordered=True)
summary = summary.sort_values(["state", "assay"]).reset_index(drop=True)
summary["assay_label"] = summary["assay"].map(ASSAY_LABELS)

summary_out = OUTPUT_DIR / f"{OUTPUT_BASENAME}_summary.tsv"
summary.to_csv(summary_out, sep="\t", index=False)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.unicode_minus": False,
})

fig, ax = plt.subplots(figsize=FIGSIZE)
ax.set_facecolor("#f7f7f7")
fig.patch.set_facecolor("white")
ax.axhline(0, color="#9aa0a6", linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)

for assay in ASSAY_ORDER:
    sub = summary[summary["assay"] == assay].copy()
    if sub.empty:
        continue
    ax.errorbar(
        sub["x_plot"],
        sub["mean_value"],
        yerr=sub["sem_value"],
        fmt="none",
        ecolor=ASSAY_COLORS[assay],
        elinewidth=1.4,
        capsize=3,
        alpha=0.95,
        zorder=3,
    )
    ax.scatter(
        sub["x_plot"],
        sub["mean_value"],
        s=40,
        color=ASSAY_COLORS[assay],
        edgecolor="white",
        linewidth=0.6,
        alpha=0.98,
        zorder=4,
        label=ASSAY_LABELS[assay],
    )

ymax = float((summary["mean_value"] + summary["sem_value"]).max())
ax.set_xlim(0.4, len(STATE_ORDER) + 0.65)
ax.set_ylim(-2, ymax * 1.10)
ax.set_xticks(range(1, len(STATE_ORDER) + 1))
ax.set_xticklabels(STATE_ORDER)
ax.set_xlabel("Chromatin state")
ax.set_ylabel("Fold enrichment")

ax.yaxis.grid(True, linestyle="-", linewidth=0.6, color="#e3e3e3")
ax.xaxis.grid(False)
ax.set_axisbelow(True)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_linewidth(1.0)
ax.spines["bottom"].set_linewidth(1.0)

legend = ax.legend(
    title=None,
    ncol=3,
    frameon=True,
    fontsize=9,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.02),
    handletextpad=0.4,
    columnspacing=1.2,
)
legend.get_frame().set_edgecolor("#cccccc")
legend.get_frame().set_linewidth(0.8)
legend.get_frame().set_facecolor("white")

plt.subplots_adjust(left=0.09, right=0.97, bottom=0.13, top=0.82)

pdf_out = OUTPUT_DIR / f"{OUTPUT_BASENAME}.pdf"
png_out = OUTPUT_DIR / f"{OUTPUT_BASENAME}.png"
plt.savefig(pdf_out, dpi=DPI, bbox_inches="tight")
plt.savefig(png_out, dpi=DPI, bbox_inches="tight")
plt.close()

print(f"[OK] Input        : {INPUT_FILE}")
print(f"[OK] Summary TSV  : {summary_out}")
print(f"[OK] Figure PDF   : {pdf_out}")
print(f"[OK] Figure PNG   : {png_out}")
