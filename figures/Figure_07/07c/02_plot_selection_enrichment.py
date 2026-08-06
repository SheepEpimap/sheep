#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

"""
 output :
1) combined  grouping
2) facet
3) median heatmap

input:
- fold_enrichment.long.tsv
- fold_enrichment.median_matrix.state_x_signal.tsv( ;  long  )

 ,  BED.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns


# =========================================================
# =========================================================
BASE = Path("/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature")
TABLE_DIR = BASE / "selection_enrichment_python" / "tables"
PLOT_DIR = BASE / "selection_enrichment_python" / "plots_replot_only"

LONG_TABLE = TABLE_DIR / "fold_enrichment.long.tsv"
MEDIAN_MATRIX_FILE = TABLE_DIR / "fold_enrichment.median_matrix.state_x_signal.tsv"

PLOT_DIR.mkdir(parents=True, exist_ok=True)


# =========================================================
# =========================================================
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


# =========================================================
# =========================================================
COMBINED_FIGSIZE = (10.4, 9.4)
FACET_FIGSIZE = (11.8, 10.2)
HEATMAP_FIGSIZE = (7.2, 5.0)

ROW_STEP = 1.35

GROUP_OFFSETS = np.array([-0.36, -0.12, 0.12, 0.36])

BOX_HEIGHT_COMBINED = 0.18
BOX_HEIGHT_FACET = 0.22

STRIP_SIZE = 10
STRIP_ALPHA = 0.22
JITTER_SD = 0.020

REFERENCE_LINE_X = 1.0
REFERENCE_LINE_COLOR = "#9E9E9E"
REFERENCE_LINE_STYLE = "--"
REFERENCE_LINE_WIDTH = 1.0

AXIS_TEXT_SIZE = 11
TITLE_SIZE = 13

RNG = np.random.default_rng(1234)


# =========================================================
# =========================================================
def add_reference_line(ax):
    """  x=1  ."""
    ax.axvline(
        REFERENCE_LINE_X,
        color=REFERENCE_LINE_COLOR,
        linestyle=REFERENCE_LINE_STYLE,
        linewidth=REFERENCE_LINE_WIDTH,
        zorder=1
    )


def add_signal_legend(ax):
    """ , ."""
    handles = [
        Patch(
            facecolor=SIGNAL_PALETTE[s],
            edgecolor="black",
            linewidth=1.0,
            alpha=1.0,
            label=s
        )
        for s in SIGNAL_ORDER
    ]

    leg = ax.legend(
        handles=handles,
        title="signal",
        bbox_to_anchor=(1.02, 1.0),
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


def format_axes(ax, y_ticks, y_ticklabels, title):
    """ ."""
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticklabels, fontsize=AXIS_TEXT_SIZE)
    ax.set_xlabel("Fold enrichment", fontsize=AXIS_TEXT_SIZE)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=TITLE_SIZE)
    ax.tick_params(axis="x", labelsize=AXIS_TEXT_SIZE - 1)
    ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.set_axisbelow(True)


def draw_horizontal_box_group(ax, values, pos, color, box_height, add_points=True):
    """
      y  , .
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return

    bp = ax.boxplot(
        values,
        positions=[pos],
        vert=False,
        widths=box_height,
        patch_artist=True,
        manage_ticks=False,
        showfliers=False,
        boxprops=dict(facecolor=color, edgecolor="black", linewidth=1.0),
        medianprops=dict(color="black", linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=0.9),
        capprops=dict(color="black", linewidth=0.9),
        zorder=3
    )

    if add_points:
        y = RNG.normal(loc=pos, scale=JITTER_SD, size=len(values))
        ax.scatter(
            values, y,
            s=STRIP_SIZE,
            color=color,
            alpha=STRIP_ALPHA,
            linewidth=0,
            zorder=4
        )


def build_state_base_positions():
    """
      state   y  .
     :
    -   E1   y,  invert_yaxis()
    -   E1,  E10
    """
    return np.arange(len(STATE_LABEL_ORDER)) * ROW_STEP


def prepare_long_df():
    """
    read long  .
    """
    if not LONG_TABLE.exists():
        raise FileNotFoundError(f"Not found long  : {LONG_TABLE}")

    df = pd.read_csv(LONG_TABLE, sep="\t")

    if "state_label" not in df.columns:
        if "state" not in df.columns:
            raise ValueError("long   state_label   state, .")
        df["state_label"] = df["state"].map(STATE_LABEL_MAP)

    df["signal"] = pd.Categorical(df["signal"], categories=SIGNAL_ORDER, ordered=True)

    df["state_label"] = pd.Categorical(
        df["state_label"],
        categories=STATE_LABEL_ORDER,
        ordered=True
    )

    df["fold_enrichment"] = pd.to_numeric(df["fold_enrichment"], errors="coerce")

    return df


def prepare_median_matrix(df):
    """
     read  median matrix;
     ,  long  .
    """
    if MEDIAN_MATRIX_FILE.exists():
        mat = pd.read_csv(MEDIAN_MATRIX_FILE, sep="\t", index_col=0)
        mat = mat.reindex(index=STATE_LABEL_ORDER, columns=SIGNAL_ORDER)
        return mat

    mat = (
        df.groupby(["state_label", "signal"], observed=False)["fold_enrichment"]
        .median()
        .unstack("signal")
        .reindex(index=STATE_LABEL_ORDER, columns=SIGNAL_ORDER)
    )
    return mat


# =========================================================
# =========================================================
def plot_combined_boxplot(df):
    """
      signal  .
      state  ,  4   signal  .
    """
    fig, ax = plt.subplots(figsize=COMBINED_FIGSIZE)

    base_positions = build_state_base_positions()

    for i, state_label in enumerate(STATE_LABEL_ORDER):
        y0 = base_positions[i]

        for j, signal_name in enumerate(SIGNAL_ORDER):
            vals = df.loc[
                (df["state_label"] == state_label) &
                (df["signal"] == signal_name),
                "fold_enrichment"
            ].dropna().to_numpy()

            pos = y0 + GROUP_OFFSETS[j]
            draw_horizontal_box_group(
                ax=ax,
                values=vals,
                pos=pos,
                color=SIGNAL_PALETTE[signal_name],
                box_height=BOX_HEIGHT_COMBINED,
                add_points=True
            )

    add_reference_line(ax)
    format_axes(
        ax=ax,
        y_ticks=base_positions,
        y_ticklabels=STATE_LABEL_ORDER,
        title="Selection-signal enrichment across chromatin states"
    )
    add_signal_legend(ax)

    ax.invert_yaxis()

    ax.set_ylim(base_positions[-1] + 0.75, -0.75)

    plt.tight_layout()
    plt.savefig(PLOT_DIR / "fold_enrichment.combined.boxplot_only.pdf")
    plt.savefig(PLOT_DIR / "fold_enrichment.combined.boxplot_only.png", dpi=300)
    plt.close(fig)


# =========================================================
# =========================================================
def plot_facet_boxplot(df):
    """
      signal  .
    """
    fig, axes = plt.subplots(2, 2, figsize=FACET_FIGSIZE, sharey=True)
    axes = axes.ravel()

    base_positions = build_state_base_positions()

    for k, signal_name in enumerate(SIGNAL_ORDER):
        ax = axes[k]
        sub = df[df["signal"] == signal_name].copy()

        for i, state_label in enumerate(STATE_LABEL_ORDER):
            vals = sub.loc[
                sub["state_label"] == state_label,
                "fold_enrichment"
            ].dropna().to_numpy()

            y0 = base_positions[i]
            draw_horizontal_box_group(
                ax=ax,
                values=vals,
                pos=y0,
                color=SIGNAL_PALETTE[signal_name],
                box_height=BOX_HEIGHT_FACET,
                add_points=True
            )

        add_reference_line(ax)
        format_axes(
            ax=ax,
            y_ticks=base_positions,
            y_ticklabels=STATE_LABEL_ORDER,
            title=signal_name
        )
        ax.invert_yaxis()
        ax.set_ylim(base_positions[-1] + 0.75, -0.75)

    for j in range(len(SIGNAL_ORDER), len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle("Selection-signal enrichment across chromatin states", fontsize=TITLE_SIZE + 1)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    fig.savefig(PLOT_DIR / "fold_enrichment.facet.boxplot_only.pdf")
    fig.savefig(PLOT_DIR / "fold_enrichment.facet.boxplot_only.png", dpi=300)
    plt.close(fig)


# =========================================================
# =========================================================
def plot_median_heatmap(median_mat):
    """
      median heatmap.
      E1 -> E10.
    """
    plt.figure(figsize=HEATMAP_FIGSIZE)
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
    ax.set_title("Median enrichment across tissues", fontsize=TITLE_SIZE)
    ax.tick_params(axis="x", labelsize=AXIS_TEXT_SIZE - 1, rotation=45)
    ax.tick_params(axis="y", labelsize=AXIS_TEXT_SIZE - 1)
    plt.tight_layout()

    plt.savefig(PLOT_DIR / "fold_enrichment.median_heatmap.replot_only.pdf")
    plt.savefig(PLOT_DIR / "fold_enrichment.median_heatmap.replot_only.png", dpi=300)
    plt.close()


# =========================================================
# =========================================================
def main():
    sns.set_theme(style="whitegrid", font_scale=0.95)

    df = prepare_long_df()
    median_mat = prepare_median_matrix(df)

    plot_combined_boxplot(df)
    plot_facet_boxplot(df)
    plot_median_heatmap(median_mat)

    print("Done.")
    print(f"Input long table : {LONG_TABLE}")
    print(f"Output plot dir  : {PLOT_DIR}")


if __name__ == "__main__":
    main()
