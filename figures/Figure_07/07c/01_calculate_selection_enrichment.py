#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

"""
 :
1. read file,  BED3,  merge
2. read  tissue   E1-E10   BED
3.   fold enrichment = (C/A) / (B/D)
4. output ,summary,heatmap,combined  ,facet

 :
-  "  overlap bp"
-   MERGE_SELECTION_INTERVALS=True,  merge
-   MERGE_STATE_INTERVALS=True,  merge state
- PLOT_STYLE  :
    "box"         ->   +
    "half_violin" ->   +   (+  )
"""

from pathlib import Path
from collections import defaultdict
import re
import sys

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns


# =========================================================
# =========================================================
BASE = Path("/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature")
STATE_DIR = Path("/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/All_chromatin_state")
TISSUE_FILE = BASE / "tissue.txt"

OUTDIR = BASE / "selection_enrichment_python"
BED_OUTDIR = OUTDIR / "selection_bed"
TABLE_OUTDIR = OUTDIR / "tables"
PLOT_OUTDIR = OUTDIR / "plots"
LOG_OUTDIR = OUTDIR / "logs"

GENOME_SIZE = 2628146905

MERGE_SELECTION_INTERVALS = True
MERGE_STATE_INTERVALS = True

SELECTION_FILES = {
    "fst_CEA_EUR": BASE / "fst_CEA_EUR" / "chrAuto.50000_10000.windowed.weir.fst.filter.t0.01.annotation.chr",
    "fst_ancient_EUR": BASE / "fst.ancient_EUR" / "ancienteurope_4000y.10000_10000.windowed.weir.fst.filter.t0.01.annotation",
    "fst_ancient_CEA": BASE / "fst.ancient_CEA" / "ancientasia_4000y.10000_10000.windowed.weir.fst.filter.t0.01.annotation.chr",
    "Demestication_domestic": BASE / "Demestication_domestic" / "Demestication_domestic_sheep.bed",
}

STATE_ORDER = [f"E{i}" for i in range(1, 11)]

STATE_LABEL_MAP = {
    "E1": "E1 TssA",
    "E2": "E2 TssFlnk",
    "E3": "E3 TSSWk",
    "E4": "E4 TssBiv",
    "E5": "E5 EnhA",
    "E6": "E6 EnhAMe",
    "E7": "E7 EnhAHet",
    "E8": "E8 EnhPois",
    "E9": "E9 Repr",
    "E10": "E10 QuiW",
}
STATE_LABEL_ORDER = [STATE_LABEL_MAP[s] for s in STATE_ORDER]

SIGNAL_ORDER = [
    "fst_CEA_EUR",
    "fst_ancient_EUR",
    "fst_ancient_CEA",
    "Demestication_domestic",
]

SIGNAL_COLOR_LIST = [
    "#DE582B",  # fst_CEA_EUR
    "#1868B2",  # fst_ancient_EUR
    "#018A67",  # fst_ancient_CEA
    "#F3A332",  # Demestication_domestic
]
SIGNAL_PALETTE = dict(zip(SIGNAL_ORDER, SIGNAL_COLOR_LIST))

PLOT_STYLE = "box"

STATE_TOP_TO_BOTTOM = STATE_LABEL_ORDER[:]

COMBINED_FIGSIZE = (9.5, 8.6)
FACET_FIGSIZE = (11.2, 9.6)

BOX_WIDTH_COMBINED = 0.50
BOX_WIDTH_FACET = 0.42
BOX_LINEWIDTH = 1.0
BOX_SATURATION = 1.0

HALF_VIOLIN_SIDE = "upper"      # "upper"   "lower"
HALF_VIOLIN_WIDTH = 0.20
HALF_VIOLIN_ALPHA = 0.60
HALF_VIOLIN_POINT_SIZE = 8
HALF_VIOLIN_JITTER = 0.015
HALF_VIOLIN_INNER_BOX = True
HALF_VIOLIN_INNER_BOX_WIDTH = 0.06

STRIP_SIZE = 1.8
STRIP_ALPHA = 0.28

REFERENCE_LINE_COLOR = "#9E9E9E"
REFERENCE_LINE_STYLE = "--"
REFERENCE_LINE_WIDTH = 1.0

RNG = np.random.default_rng(1234)


# =========================================================
# =========================================================
def ensure_dirs():
    """ Output directory."""
    for d in [OUTDIR, BED_OUTDIR, TABLE_OUTDIR, PLOT_OUTDIR, LOG_OUTDIR]:
        d.mkdir(parents=True, exist_ok=True)


def normalize_chrom(chrom: str) -> str:
    """
     ,  chrchr12  .
     :
    1) chrchr12 -> chr12
    2)   / X / Y / M / MT,  chr*
    3)
    """
    chrom = chrom.strip()

    if chrom.startswith("chrchr"):
        chrom = chrom[3:]   # chrchr12 -> chr12

    if chrom.startswith("chr"):
        return chrom

    if re.fullmatch(r"(\d+|X|Y|M|MT)", chrom):
        return f"chr{chrom}"

    return chrom


def read_tissues(tissue_file: Path):
    """
    read tissue.txt
     : ; , .
    """
    tissues = []
    if not tissue_file.exists():
        return tissues

    with tissue_file.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tissue = line.split()[0]
            if tissue.lower() == "tissue":
                continue
            tissues.append(tissue)

    return tissues


def infer_tissues_from_state_dir(state_dir: Path):
    """
      tissue.txt  ,  *_E*.bed file .
    file :
        abomasum_E1.bed
        lymph-node_E7.bed
    """
    tissues = set()
    for bed in state_dir.glob("*_E*.bed"):
        m = re.match(r"(.+)_E\d+$", bed.stem)
        if m:
            tissues.add(m.group(1))
    return sorted(tissues)


def read_bed3_as_dict(path: Path):
    """
    read BED-like file,  3  , :
        {chrom: [(start, end), ...]}
     :
    -   /   /
    -   chrchr*
    - end <= start
    """
    d = defaultdict(list)

    if not path.exists():
        raise FileNotFoundError(f"file : {path}")

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = re.split(r"\s+", line)
            if len(parts) < 3:
                continue

            chrom = parts[0]

            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            if end <= start:
                continue

            chrom = normalize_chrom(chrom)
            d[chrom].append((start, end))

    return d


def merge_intervals(interval_dict):
    """
      merge.
     :
        {chrom: [(merged_start, merged_end), ...]}
    """
    merged = {}

    for chrom, intervals in interval_dict.items():
        if not intervals:
            merged[chrom] = []
            continue

        intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
        out = [list(intervals[0])]

        for s, e in intervals[1:]:
            last_s, last_e = out[-1]
            if s <= last_e:
                out[-1][1] = max(last_e, e)
            else:
                out.append([s, e])

        merged[chrom] = [(s, e) for s, e in out]

    return merged


def total_length(interval_dict):
    """ ."""
    total = 0
    for chrom, intervals in interval_dict.items():
        for s, e in intervals:
            total += (e - s)
    return total


def overlap_length(a_dict, b_dict):
    """
      overlap bp, .
     :
    -  ;  merge
    """
    total = 0
    common_chroms = set(a_dict.keys()) & set(b_dict.keys())

    for chrom in common_chroms:
        a = a_dict[chrom]
        b = b_dict[chrom]
        i, j = 0, 0

        while i < len(a) and j < len(b):
            a_s, a_e = a[i]
            b_s, b_e = b[j]

            start = max(a_s, b_s)
            end = min(a_e, b_e)

            if end > start:
                total += (end - start)

            if a_e <= b_e:
                i += 1
            else:
                j += 1

    return total


def write_bed(interval_dict, out_path: Path):
    """  BED3 file."""
    with out_path.open("w") as out:
        for chrom in sorted(interval_dict.keys()):
            for s, e in interval_dict[chrom]:
                out.write(f"{chrom}\t{s}\t{e}\n")


def fold_enrichment(C, A, B, D):
    """enrichment = (C/A) / (B/D)"""
    if A == 0 or B == 0 or D == 0:
        return np.nan
    return (C / A) / (B / D)


def get_horizontal_plot_order():
    """
    seaborn  ,order  ,
      E1 -> E10.
    """
    return list(STATE_TOP_TO_BOTTOM)


def add_reference_line(ax):
    """  fold enrichment = 1  ."""
    ax.axvline(
        1,
        color=REFERENCE_LINE_COLOR,
        linestyle=REFERENCE_LINE_STYLE,
        linewidth=REFERENCE_LINE_WIDTH
    )


def add_signal_legend(ax, hue_order):
    """
     :
    1)
    2)
    """
    handles = [
        Patch(
            facecolor=SIGNAL_PALETTE[h],
            edgecolor="black",
            linewidth=1.0,
            alpha=1.0,
            label=h
        )
        for h in hue_order
    ]

    leg = ax.legend(
        handles=handles,
        title="signal",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        borderpad=0.4,
        labelspacing=0.5
    )

    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(1.0)
    leg.get_frame().set_facecolor("white")


def draw_half_violin(ax, values, pos, color):
    """
     .
    values:   enrichment
    pos:    y
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return

    if len(values) >= 2:
        vp = ax.violinplot(
            values,
            positions=[pos],
            vert=False,
            widths=HALF_VIOLIN_WIDTH,
            showmeans=False,
            showmedians=False,
            showextrema=False
        )

        body = vp["bodies"][0]
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(HALF_VIOLIN_ALPHA)
        body.set_linewidth(0.8)

        verts = body.get_paths()[0].vertices
        if HALF_VIOLIN_SIDE == "upper":
            verts[:, 1] = np.maximum(verts[:, 1], pos)
        else:
            verts[:, 1] = np.minimum(verts[:, 1], pos)

    if HALF_VIOLIN_INNER_BOX and len(values) >= 2:
        ax.boxplot(
            values,
            positions=[pos],
            vert=False,
            widths=HALF_VIOLIN_INNER_BOX_WIDTH,
            manage_ticks=False,
            showfliers=False,
            patch_artist=True,
            boxprops=dict(facecolor="white", edgecolor=color, linewidth=0.9),
            medianprops=dict(color=color, linewidth=1.2),
            whiskerprops=dict(color=color, linewidth=0.8),
            capprops=dict(color=color, linewidth=0.8),
        )

    jitter_center = pos + (0.03 if HALF_VIOLIN_SIDE == "upper" else -0.03)
    ys = RNG.normal(loc=jitter_center, scale=HALF_VIOLIN_JITTER, size=len(values))
    ax.scatter(values, ys, s=HALF_VIOLIN_POINT_SIZE, color=color, alpha=0.25, linewidth=0)


def draw_half_violin_panel(ax, data, order, hue_order):
    """
     ,  state × signal  .
      signal   state   dodge.
    """
    if len(hue_order) == 1:
        offsets = np.array([0.0])
    else:
        offsets = np.linspace(-0.30, 0.30, len(hue_order))

    for state_idx, state_label in enumerate(order):
        for hue_idx, signal_name in enumerate(hue_order):
            vals = data.loc[
                (data["state_label"] == state_label) &
                (data["signal"] == signal_name),
                "fold_enrichment"
            ].dropna().to_numpy()

            if len(vals) == 0:
                continue

            pos = state_idx + offsets[hue_idx]
            draw_half_violin(ax, vals, pos, SIGNAL_PALETTE[signal_name])

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.invert_yaxis()


# =========================================================
# =========================================================
def prepare_selection_signals():
    """
    read file, , output  BED.
     :
        signal_intervals: {signal_name: interval_dict}
        signal_lengths:   {signal_name: total_bp}
        signal_meta_df:   DataFrame
    """
    signal_intervals = {}
    meta_rows = []

    for signal_name in SIGNAL_ORDER:
        path = SELECTION_FILES[signal_name]
        raw_dict = read_bed3_as_dict(path)
        raw_n = sum(len(v) for v in raw_dict.values())
        raw_bp = total_length(raw_dict)

        if MERGE_SELECTION_INTERVALS:
            final_dict = merge_intervals(raw_dict)
        else:
            final_dict = {k: sorted(v) for k, v in raw_dict.items()}

        final_n = sum(len(v) for v in final_dict.values())
        final_bp = total_length(final_dict)

        signal_intervals[signal_name] = final_dict

        out_bed = BED_OUTDIR / f"{signal_name}.bed"
        write_bed(final_dict, out_bed)

        meta_rows.append({
            "signal": signal_name,
            "source_file": str(path),
            "raw_interval_n": raw_n,
            "raw_total_bp": raw_bp,
            "final_interval_n": final_n,
            "final_total_bp": final_bp,
            "merged": MERGE_SELECTION_INTERVALS,
            "standard_bed": str(out_bed)
        })

    signal_meta_df = pd.DataFrame(meta_rows)
    signal_lengths = {row["signal"]: row["final_total_bp"] for _, row in signal_meta_df.iterrows()}

    return signal_intervals, signal_lengths, signal_meta_df


# =========================================================
# =========================================================
def compute_all_enrichment(tissues, signal_intervals, signal_lengths):
    """
      tissue   E1-E10   enrichment.
      DataFrame.
    """
    results = []
    missing_files = []

    for tissue in tissues:
        for state in STATE_ORDER:
            state_file = STATE_DIR / f"{tissue}_{state}.bed"

            if not state_file.exists():
                missing_files.append(str(state_file))
                continue

            state_dict = read_bed3_as_dict(state_file)
            raw_n = sum(len(v) for v in state_dict.values())
            raw_bp = total_length(state_dict)

            if MERGE_STATE_INTERVALS:
                state_dict = merge_intervals(state_dict)
            else:
                state_dict = {k: sorted(v) for k, v in state_dict.items()}

            final_n = sum(len(v) for v in state_dict.values())
            A = total_length(state_dict)

            for signal_name in SIGNAL_ORDER:
                sel_dict = signal_intervals[signal_name]
                B = signal_lengths[signal_name]
                C = overlap_length(state_dict, sel_dict)
                FE = fold_enrichment(C=C, A=A, B=B, D=GENOME_SIZE)

                results.append({
                    "tissue": tissue,
                    "state": state,
                    "state_label": STATE_LABEL_MAP[state],
                    "signal": signal_name,
                    "state_file": str(state_file),
                    "A_state_bp": A,
                    "B_signal_bp": B,
                    "C_overlap_bp": C,
                    "D_background_bp": GENOME_SIZE,
                    "state_raw_interval_n": raw_n,
                    "state_final_interval_n": final_n,
                    "state_raw_total_bp": raw_bp,
                    "fold_enrichment": FE
                })

    if missing_files:
        with (LOG_OUTDIR / "missing_state_files.txt").open("w") as out:
            for x in missing_files:
                out.write(x + "\n")

    df = pd.DataFrame(results)
    if df.empty:
        raise RuntimeError("  enrichment results, check tissue.txt,state file .")

    df["state_label"] = pd.Categorical(
        df["state_label"],
        categories=STATE_LABEL_ORDER,
        ordered=True
    )
    df["signal"] = pd.Categorical(
        df["signal"],
        categories=SIGNAL_ORDER,
        ordered=True
    )

    return df


# =========================================================
# =========================================================
def write_tables(df, signal_meta_df):
    """
    output:
    1)
    2) summary
    3) signal x state   median
    4)   signal   tissue x state
    """
    signal_meta_df.to_csv(
        TABLE_OUTDIR / "selection_signal_meta.tsv",
        sep="\t", index=False
    )

    df.to_csv(
        TABLE_OUTDIR / "fold_enrichment.long.tsv",
        sep="\t", index=False
    )

    summary_df = (
        df.groupby(["signal", "state", "state_label"], observed=False)["fold_enrichment"]
        .agg(
            n="count",
            mean="mean",
            median="median",
            sd="std",
            min="min",
            max="max"
        )
        .reset_index()
    )
    summary_df.to_csv(
        TABLE_OUTDIR / "fold_enrichment.summary_by_signal_state.tsv",
        sep="\t", index=False
    )

    median_mat = (
        df.groupby(["state_label", "signal"], observed=False)["fold_enrichment"]
        .median()
        .unstack("signal")
        .reindex(STATE_LABEL_ORDER)
        .reindex(columns=SIGNAL_ORDER)
    )
    median_mat.to_csv(
        TABLE_OUTDIR / "fold_enrichment.median_matrix.state_x_signal.tsv",
        sep="\t"
    )

    for signal_name in SIGNAL_ORDER:
        sub = df[df["signal"] == signal_name].copy()
        mat = sub.pivot(index="tissue", columns="state_label", values="fold_enrichment")
        mat = mat.reindex(columns=STATE_LABEL_ORDER)
        mat.to_csv(
            TABLE_OUTDIR / f"{signal_name}.tissue_x_state.tsv",
            sep="\t"
        )

    return summary_df, median_mat


# =========================================================
# =========================================================
def plot_combined_boxplot(df):
    """
      signal  .
      PLOT_STYLE  :
        - "box"
        - "half_violin"
    """
    order = get_horizontal_plot_order()

    plt.figure(figsize=COMBINED_FIGSIZE)
    ax = plt.gca()

    if PLOT_STYLE == "box":
        sns.boxplot(
            data=df,
            y="state_label",
            x="fold_enrichment",
            hue="signal",
            order=order,
            hue_order=SIGNAL_ORDER,
            palette=SIGNAL_PALETTE,
            showfliers=False,
            width=BOX_WIDTH_COMBINED,
            linewidth=BOX_LINEWIDTH,
            saturation=BOX_SATURATION,
            ax=ax
        )

        sns.stripplot(
            data=df,
            y="state_label",
            x="fold_enrichment",
            hue="signal",
            order=order,
            hue_order=SIGNAL_ORDER,
            dodge=True,
            palette=SIGNAL_PALETTE,
            alpha=STRIP_ALPHA,
            size=STRIP_SIZE,
            linewidth=0,
            ax=ax
        )

        if ax.legend_ is not None:
            ax.legend_.remove()
        add_signal_legend(ax, SIGNAL_ORDER)

    elif PLOT_STYLE == "half_violin":
        draw_half_violin_panel(ax, df, order, SIGNAL_ORDER)
        add_signal_legend(ax, SIGNAL_ORDER)

    else:
        raise ValueError(f"  PLOT_STYLE: {PLOT_STYLE}")

    add_reference_line(ax)
    ax.set_xlabel("Fold enrichment")
    ax.set_ylabel("")
    ax.set_title(f"Selection-signal enrichment across chromatin states ({PLOT_STYLE})")
    ax.margins(y=0.03)

    plt.tight_layout()
    plt.savefig(PLOT_OUTDIR / f"fold_enrichment.combined.{PLOT_STYLE}.pdf")
    plt.savefig(PLOT_OUTDIR / f"fold_enrichment.combined.{PLOT_STYLE}.png", dpi=300)
    plt.close()


def plot_facet_boxplot(df):
    """
      signal  .
      PLOT_STYLE  :
        - "box"
        - "half_violin"
    """
    order = get_horizontal_plot_order()

    fig, axes = plt.subplots(2, 2, figsize=FACET_FIGSIZE, sharey=True)
    axes = axes.ravel()

    for i, signal_name in enumerate(SIGNAL_ORDER):
        ax = axes[i]
        sub = df[df["signal"] == signal_name].copy()

        if PLOT_STYLE == "box":
            sns.boxplot(
                data=sub,
                y="state_label",
                x="fold_enrichment",
                order=order,
                color=SIGNAL_PALETTE[signal_name],
                showfliers=False,
                width=BOX_WIDTH_FACET,
                linewidth=BOX_LINEWIDTH,
                saturation=BOX_SATURATION,
                ax=ax
            )

            sns.stripplot(
                data=sub,
                y="state_label",
                x="fold_enrichment",
                order=order,
                color=SIGNAL_PALETTE[signal_name],
                alpha=STRIP_ALPHA,
                size=STRIP_SIZE,
                linewidth=0,
                ax=ax
            )

        elif PLOT_STYLE == "half_violin":
            draw_half_violin_panel(ax, sub, order, [signal_name])

        else:
            raise ValueError(f"  PLOT_STYLE: {PLOT_STYLE}")

        add_reference_line(ax)
        ax.set_title(signal_name)
        ax.set_xlabel("Fold enrichment")
        ax.set_ylabel("")
        ax.margins(y=0.03)

    for j in range(len(SIGNAL_ORDER), len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle(
        f"Selection-signal enrichment across chromatin states ({PLOT_STYLE})",
        fontsize=13
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    fig.savefig(PLOT_OUTDIR / f"fold_enrichment.facet.{PLOT_STYLE}.pdf")
    fig.savefig(PLOT_OUTDIR / f"fold_enrichment.facet.{PLOT_STYLE}.png", dpi=300)
    plt.close(fig)


def plot_median_heatmap(median_mat):
    """
      signal x state  .
    state   E1 -> E10.
    """
    plt.figure(figsize=(6.8, 4.8))
    ax = sns.heatmap(
        median_mat,
        cmap="RdYlBu_r",
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Median fold enrichment"}
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Median enrichment across tissues")
    plt.tight_layout()

    plt.savefig(PLOT_OUTDIR / "fold_enrichment.median_heatmap.pdf")
    plt.savefig(PLOT_OUTDIR / "fold_enrichment.median_heatmap.png", dpi=300)
    plt.close()


# =========================================================
# =========================================================
def main():
    ensure_dirs()

    tissues = read_tissues(TISSUE_FILE)
    if not tissues:
        tissues = infer_tissues_from_state_dir(STATE_DIR)
        sys.stderr.write(
            "[WARN] tissue.txt  ,  All_chromatin_state   tissue  .\n"
        )

    with (LOG_OUTDIR / "tissues_used.txt").open("w") as out:
        for t in tissues:
            out.write(t + "\n")

    signal_intervals, signal_lengths, signal_meta_df = prepare_selection_signals()
    df = compute_all_enrichment(tissues, signal_intervals, signal_lengths)
    summary_df, median_mat = write_tables(df, signal_meta_df)

    sns.set_theme(style="whitegrid", font_scale=0.95)
    plot_combined_boxplot(df)
    plot_facet_boxplot(df)
    plot_median_heatmap(median_mat)

    print("Done.")
    print(f"Output dir: {OUTDIR}")
    print(f"Long table: {TABLE_OUTDIR / 'fold_enrichment.long.tsv'}")
    print(f"Summary   : {TABLE_OUTDIR / 'fold_enrichment.summary_by_signal_state.tsv'}")
    print(f"Plots dir : {PLOT_OUTDIR}")
    print(f"PLOT_STYLE: {PLOT_STYLE}")
    print(f"Colors    : {SIGNAL_PALETTE}")


if __name__ == "__main__":
    main()
