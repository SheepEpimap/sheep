#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import itertools
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
mpl.rcParams["axes.linewidth"] = 1.2
mpl.rcParams["xtick.major.width"] = 1.1
mpl.rcParams["ytick.major.width"] = 1.1
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"

import matplotlib.pyplot as plt
from scipy import stats

try:
    import pyBigWig
except ImportError:
    raise ImportError(
        "pyBigWig was not detected. Install it with: conda install -c bioconda pybigwig -y"
    )


# ============================================================
# 1. inputfile
# ============================================================

PAIR_FILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_linked_enhancer_gene_pairs.tsv"
)

EXPR_FILE = "/vol2/mengzhu/snakemake_sheep/expressiondir/rna/all_tisssues_expression_tpm.csv"

TAU_FILE = "/vol2/mengzhu/snakemake_sheep/expressiondir/all_tisssues_expression_tpm.median.tau.csv"

PHYLOP_BW = (
    "/vol2/zhangshiwen/sheep_cor/comparative_analysis_v3/"
    "G2_per_tissue/conservation_bw/phyloP_v2.bw"
)

PHASTCONS_BW = (
    "/vol2/zhangshiwen/sheep_cor/comparative_analysis_v3/"
    "G2_per_tissue/conservation_bw/phastCons_v2.bw"
)

GERP_BW = (
    "/data/home/sczd644/run/zsw_chrombpnet/phylop/v3_v2_grep/"
    "gerp_conservation_scores.ovis_aries.ARS-UI_Ramb_v2.0.bw"
)


# ============================================================
# 2. Output directory
# ============================================================

OUTROOT = "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/zutu"

OUTDIR = os.path.join(
    OUTROOT,
    "E5_5groups_EnhA_gene_expression_tau_enhancer_conservation_remove_outliers"
)

FIG_DIR = os.path.join(OUTDIR, "figures")
TABLE_DIR = os.path.join(OUTDIR, "tables")
CHECK_DIR = os.path.join(OUTDIR, "check")

for d in [OUTDIR, FIG_DIR, TABLE_DIR, CHECK_DIR]:
    os.makedirs(d, exist_ok=True)


# ============================================================
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

GROUP_FULL_NAME = {
    "TS-EnhA": "Tissue-specific enhancer activity",
    "NS-EnhA": "Narrowly shared enhancer activity",
    "MS-EnhA": "Moderately shared enhancer activity",
    "BS-EnhA": "Broadly shared enhancer activity",
    "ES-EnhA": "Extensively shared enhancer activity",
}

GROUP_TISSUE_RANGE = {
    "TS-EnhA": "1 tissue",
    "NS-EnhA": "2-5 tissues",
    "MS-EnhA": "6-10 tissues",
    "BS-EnhA": "11-20 tissues",
    "ES-EnhA": "21-43 tissues",
}

PLOT_GROUPS = [GROUP_LABEL[g] for g in GROUP_ORDER]

GROUP_COLORS = {
    "TS-EnhA": "#0072B2",
    "NS-EnhA": "#E69F00",
    "MS-EnhA": "#009E73",
    "BS-EnhA": "#D55E00",
    "ES-EnhA": "#CC79A7",
}


# ============================================================
# ============================================================

USE_LOG10_EXPR = True
ALPHA = 0.05

GENE_LEVEL_METRICS = [
    "expression_for_plot",
    "tau",
]

ENHANCER_LEVEL_METRICS = [
    "enhancer_phyloP_mean",
    "enhancer_phastCons_mean",
    "enhancer_GERP_mean",
]

METRICS_TO_PLOT = GENE_LEVEL_METRICS + ENHANCER_LEVEL_METRICS

METRIC_LEVEL = {
    "expression_for_plot": "target_gene",
    "tau": "target_gene",
    "enhancer_phyloP_mean": "enhancer_region",
    "enhancer_phastCons_mean": "enhancer_region",
    "enhancer_GERP_mean": "enhancer_region",
}

METRIC_LABELS = {
    "expression_for_plot": "Average target-gene expression\nlog10(TPM + 1)",
    "tau": "Target-gene tau value",
    "enhancer_phyloP_mean": "Mean phyloP score\nenhancer region",
    "enhancer_phastCons_mean": "Mean phastCons score\nenhancer region",
    "enhancer_GERP_mean": "Mean GERP score\nenhancer region",
}

METRIC_PREFIX = {
    "expression_for_plot": "target_gene_expression_EnhA_5groups",
    "tau": "target_gene_tau_EnhA_5groups",
    "enhancer_phyloP_mean": "enhancer_region_phyloP_EnhA_5groups",
    "enhancer_phastCons_mean": "enhancer_region_phastCons_EnhA_5groups",
    "enhancer_GERP_mean": "enhancer_region_GERP_EnhA_5groups",
}


# ============================================================
# ============================================================

REMOVE_OUTLIER_METRICS = {
    "enhancer_phyloP_mean",
    "enhancer_GERP_mean",
}

OUTLIER_METHOD = "global_iqr"
OUTLIER_IQR_MULTIPLIER = 1.5


# ============================================================
# 5. helper functions
# ============================================================

def read_table_auto(file_path, header="infer"):
    """
     read file.
     :
    1. tab
    2. comma
    3.
    """
    seps = ["\t", ",", r"\s+"]
    best_df = None
    best_ncol = 0

    for sep in seps:
        try:
            df = pd.read_csv(
                file_path,
                sep=sep,
                header=header,
                dtype=str,
                engine="python",
            )
        except Exception:
            continue

        if df.shape[1] > best_ncol:
            best_df = df
            best_ncol = df.shape[1]

        if df.shape[1] > 1:
            best_df = df
            break

    if best_df is None:
        raise ValueError(f" readfile: {file_path}")

    df = best_df.dropna(how="all").copy()
    df.columns = [str(c).strip() for c in df.columns]

    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    return df


def bh_adjust(pvalues):
    """
    Benjamini-Hochberg FDR  .
    """
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)

    if n == 0:
        return np.array([])

    order = np.argsort(pvalues)
    ranked_p = pvalues[order]

    adjusted = np.empty(n, dtype=float)
    prev = 1.0

    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked_p[i] * n / rank
        prev = min(prev, val)
        adjusted[i] = prev

    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adjusted, 1.0)

    return out


def compact_letter_display(pmat, alpha=0.05):
    """
      pairwise adjusted p-value  .
    p < alpha: , .
    p >= alpha: , .
    """
    groups = list(pmat.index)
    letters = {g: "" for g in groups}
    letter_pool = list("abcdefghijklmnopqrstuvwxyz")

    for g in groups:
        placed = False

        for L in letter_pool:
            holders = [h for h in groups if L in letters[h]]

            if len(holders) == 0:
                letters[g] += L
                placed = True
                break

            can_share = True

            for h in holders:
                if pmat.loc[g, h] < alpha:
                    can_share = False
                    break

            if can_share:
                letters[g] += L
                placed = True
                break

        if not placed:
            raise RuntimeError(" , grouping .")

    return letters


def pairwise_tests_for_letters(df, value_col):
    """
    Kruskal-Wallis +   Mann-Whitney U + BH   +  .
    """
    sub = df[["group", value_col]].dropna().copy()
    sub["group"] = pd.Categorical(sub["group"], categories=PLOT_GROUPS, ordered=True)

    arrays = []
    valid_groups = []

    for g in PLOT_GROUPS:
        values = sub.loc[sub["group"] == g, value_col].dropna().astype(float).values

        if len(values) > 0:
            arrays.append(values)
            valid_groups.append(g)

    if len(valid_groups) < 2:
        return pd.DataFrame(), {g: "a" for g in PLOT_GROUPS}, np.nan

    try:
        kruskal_p = stats.kruskal(*arrays).pvalue
    except ValueError:
        kruskal_p = 1.0

    raw_records = []

    for g1, g2 in itertools.combinations(PLOT_GROUPS, 2):
        x = sub.loc[sub["group"] == g1, value_col].dropna().astype(float).values
        y = sub.loc[sub["group"] == g2, value_col].dropna().astype(float).values

        if len(x) == 0 or len(y) == 0:
            p = np.nan
        else:
            p = stats.mannwhitneyu(
                x,
                y,
                alternative="two-sided"
            ).pvalue

        raw_records.append([g1, g2, p])

    pairwise = pd.DataFrame(raw_records, columns=["group1", "group2", "p_raw"])
    pairwise["p_adj"] = np.nan

    non_nan_mask = pairwise["p_raw"].notna()

    if non_nan_mask.sum() > 0:
        pairwise.loc[non_nan_mask, "p_adj"] = bh_adjust(
            pairwise.loc[non_nan_mask, "p_raw"].values
        )

    pmat = pd.DataFrame(
        1.0,
        index=PLOT_GROUPS,
        columns=PLOT_GROUPS,
        dtype=float,
    )

    for _, row in pairwise.iterrows():
        if pd.notna(row["p_adj"]):
            pmat.loc[row["group1"], row["group2"]] = row["p_adj"]
            pmat.loc[row["group2"], row["group1"]] = row["p_adj"]

    letters = compact_letter_display(pmat, alpha=ALPHA)

    return pairwise, letters, kruskal_p


def make_summary(df, value_col, metric):
    """
      summary   pairwise  .
    """
    pairwise, letters, kruskal_p = pairwise_tests_for_letters(df, value_col)

    rows = []

    for g in PLOT_GROUPS:
        values = df.loc[df["group"] == g, value_col].dropna().astype(float)

        if len(values) == 0:
            rows.append({
                "metric": metric,
                "level": METRIC_LEVEL[metric],
                "group": g,
                "group_full_name": GROUP_FULL_NAME[g],
                "tissue_range": GROUP_TISSUE_RANGE[g],
                "n": 0,
                "mean": np.nan,
                "median": np.nan,
                "sd": np.nan,
                "q25": np.nan,
                "q75": np.nan,
                "min": np.nan,
                "max": np.nan,
                "kruskal_p": kruskal_p,
                "letter": "NA",
            })
        else:
            rows.append({
                "metric": metric,
                "level": METRIC_LEVEL[metric],
                "group": g,
                "group_full_name": GROUP_FULL_NAME[g],
                "tissue_range": GROUP_TISSUE_RANGE[g],
                "n": int(len(values)),
                "mean": values.mean(),
                "median": values.median(),
                "sd": values.std(ddof=1),
                "q25": values.quantile(0.25),
                "q75": values.quantile(0.75),
                "min": values.min(),
                "max": values.max(),
                "kruskal_p": kruskal_p,
                "letter": letters.get(g, ""),
            })

    summary = pd.DataFrame(rows)

    if pairwise.shape[0] > 0:
        pairwise.insert(0, "metric", metric)
        pairwise.insert(1, "level", METRIC_LEVEL[metric])

    return summary, pairwise


def nature_style_ax(ax):
    """
     .
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)

    ax.tick_params(axis="both", which="major", labelsize=10, length=4, width=1.0)
    ax.grid(False)


def remove_outliers_global_iqr(df, value_col, metric):
    """
      IQR  .

     :
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR

     :
    1.  ,  Q1/Q3;
    2.  , ;
    3.   filtered_df,outlier_df,setting_dict.
    """
    if value_col not in df.columns:
        raise ValueError(f"{value_col}  .")

    work = df.copy()
    values = pd.to_numeric(work[value_col], errors="coerce")

    valid_mask = values.notna()

    if valid_mask.sum() == 0:
        setting = {
            "metric": metric,
            "method": OUTLIER_METHOD,
            "iqr_multiplier": OUTLIER_IQR_MULTIPLIER,
            "q1": np.nan,
            "q3": np.nan,
            "iqr": np.nan,
            "lower_bound": np.nan,
            "upper_bound": np.nan,
            "n_total_rows": work.shape[0],
            "n_valid_values": 0,
            "n_removed_outliers": 0,
            "n_kept_after_filter": 0,
        }

        return work, work.iloc[0:0].copy(), setting

    q1 = values[valid_mask].quantile(0.25)
    q3 = values[valid_mask].quantile(0.75)
    iqr = q3 - q1

    if pd.isna(iqr) or iqr == 0:
        lower = values[valid_mask].min()
        upper = values[valid_mask].max()
    else:
        lower = q1 - OUTLIER_IQR_MULTIPLIER * iqr
        upper = q3 + OUTLIER_IQR_MULTIPLIER * iqr

    outlier_mask = valid_mask & ((values < lower) | (values > upper))
    keep_mask = ~outlier_mask

    filtered_df = work.loc[keep_mask].copy()
    outlier_df = work.loc[outlier_mask].copy()

    if outlier_df.shape[0] > 0:
        outlier_df["outlier_metric"] = metric
        outlier_df["outlier_method"] = OUTLIER_METHOD
        outlier_df["outlier_iqr_multiplier"] = OUTLIER_IQR_MULTIPLIER
        outlier_df["outlier_lower_bound"] = lower
        outlier_df["outlier_upper_bound"] = upper

    setting = {
        "metric": metric,
        "method": OUTLIER_METHOD,
        "iqr_multiplier": OUTLIER_IQR_MULTIPLIER,
        "q1": q1,
        "q3": q3,
        "iqr": iqr,
        "lower_bound": lower,
        "upper_bound": upper,
        "n_total_rows": work.shape[0],
        "n_valid_values": int(valid_mask.sum()),
        "n_removed_outliers": int(outlier_mask.sum()),
        "n_kept_after_filter": int(filtered_df[value_col].notna().sum()),
    }

    return filtered_df, outlier_df, setting


# ============================================================
# ============================================================

def draw_violin_on_ax(ax, df, value_col, y_label, summary_df, panel_label=None):
    """
      ax   violin + boxplot +  .

     :
    1.   median  ;
    2.   y  ;
    3. c/e  ;
    4.  statistics .
    """
    violin_data = []
    violin_positions = []
    violin_groups = []

    box_data = []
    box_positions = []
    box_groups = []

    for i, g in enumerate(PLOT_GROUPS, start=1):
        values = df.loc[df["group"] == g, value_col].dropna().astype(float).values

        if len(values) == 0:
            continue

        box_data.append(values)
        box_positions.append(i)
        box_groups.append(g)

        if len(values) >= 2:
            violin_data.append(values)
            violin_positions.append(i)
            violin_groups.append(g)

    if len(box_data) == 0:
        ax.text(
            0.5,
            0.5,
            f"No valid values\nfor {value_col}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10
        )

        ax.set_xticks(np.arange(1, len(PLOT_GROUPS) + 1))
        ax.set_xticklabels(PLOT_GROUPS, rotation=35, ha="right")
        ax.set_xlabel("")
        ax.set_ylabel(y_label, fontsize=10)

        nature_style_ax(ax)

        if panel_label is not None:
            ax.text(
                -0.18,
                1.05,
                panel_label,
                transform=ax.transAxes,
                fontsize=13,
                fontweight="bold",
                va="top",
                ha="left"
            )

        return

    if len(violin_data) > 0:
        parts = ax.violinplot(
            violin_data,
            positions=violin_positions,
            widths=0.82,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

        for i, body in enumerate(parts["bodies"]):
            g = violin_groups[i]
            body.set_facecolor(GROUP_COLORS[g])
            body.set_edgecolor("none")
            body.set_alpha(0.82)

    ax.boxplot(
        box_data,
        positions=box_positions,
        widths=0.18,
        patch_artist=True,
        showfliers=False,
        showcaps=True,
        boxprops=dict(facecolor="white", edgecolor="black", linewidth=0.7),
        medianprops=dict(color="black", linewidth=0.8),
        whiskerprops=dict(color="black", linewidth=0.6),
        capprops=dict(color="black", linewidth=0.6),
    )

    for pos, values, g in zip(box_positions, box_data, box_groups):
        if len(values) == 1:
            ax.scatter(
                [pos],
                values,
                s=20,
                color=GROUP_COLORS[g],
                edgecolor="black",
                linewidth=0.4,
                zorder=3
            )

    ax.set_xticks(np.arange(1, len(PLOT_GROUPS) + 1))
    ax.set_xticklabels(PLOT_GROUPS, rotation=35, ha="right")

    ax.set_xlabel("")
    ax.set_ylabel(y_label, fontsize=10)

    nature_style_ax(ax)

    all_values = df[value_col].dropna().astype(float)

    ymin = all_values.min()
    ymax = all_values.max()
    y_range = ymax - ymin if ymax > ymin else 1

    letter_y = ymax + 0.07 * y_range
    ax.set_ylim(ymin - 0.04 * y_range, ymax + 0.20 * y_range)

    letters = dict(zip(summary_df["group"], summary_df["letter"]))

    for i, g in enumerate(PLOT_GROUPS, start=1):
        if g in box_groups:
            ax.text(
                i,
                letter_y,
                letters.get(g, ""),
                ha="center",
                va="bottom",
                fontsize=11,
                color="black",
            )
        else:
            ax.text(
                i,
                letter_y,
                "NA",
                ha="center",
                va="bottom",
                fontsize=9,
                color="black",
            )

    if panel_label is not None:
        ax.text(
            -0.18,
            1.05,
            panel_label,
            transform=ax.transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
            ha="left"
        )


def plot_violin(df, value_col, y_label, summary_df, out_prefix):
    """
    output .
    """
    fig, ax = plt.subplots(figsize=(4.9, 3.9))

    draw_violin_on_ax(
        ax=ax,
        df=df,
        value_col=value_col,
        y_label=y_label,
        summary_df=summary_df,
        panel_label=None
    )

    ax.set_xlabel("EnhA tissue-breadth group", fontsize=11)

    plt.tight_layout()

    pdf_file = os.path.join(FIG_DIR, f"{out_prefix}.violin.pdf")
    png_file = os.path.join(FIG_DIR, f"{out_prefix}.violin.png")
    svg_file = os.path.join(FIG_DIR, f"{out_prefix}.violin.svg")

    plt.savefig(pdf_file, bbox_inches="tight")
    plt.savefig(png_file, dpi=600, bbox_inches="tight")
    plt.savefig(svg_file, bbox_inches="tight")
    plt.close()


def plot_combined_violin_figure(data_by_metric, summary_all):
    """
      expression,tau,enhancer phyloP,enhancer phastCons,enhancer GERP
      2 x 3 panel  .
    """
    panel_metrics = [
        ("expression_for_plot", "Average target-gene expression\nlog10(TPM + 1)", "a"),
        ("tau", "Target-gene tau value", "b"),
        ("enhancer_phyloP_mean", "Mean phyloP score\nenhancer region", "c"),
        ("enhancer_phastCons_mean", "Mean phastCons score\nenhancer region", "d"),
        ("enhancer_GERP_mean", "Mean GERP score\nenhancer region", "e"),
    ]

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(11.2, 6.8)
    )

    axes = axes.flatten()

    for i, (metric, y_label, panel_label) in enumerate(panel_metrics):
        metric_summary = summary_all[summary_all["metric"] == metric].copy()
        metric_df = data_by_metric[metric]

        draw_violin_on_ax(
            ax=axes[i],
            df=metric_df,
            value_col=metric,
            y_label=y_label,
            summary_df=metric_summary,
            panel_label=panel_label
        )

    axes[5].axis("off")

    fig.text(
        0.5,
        0.02,
        "EnhA tissue-breadth group",
        ha="center",
        va="center",
        fontsize=12
    )

    plt.tight_layout(rect=[0.02, 0.05, 1, 1])

    combined_pdf = os.path.join(
        FIG_DIR,
        "combined_target_gene_expression_tau_enhancer_phyloP_phastCons_GERP_EnhA_remove_outliers.violin.pdf"
    )

    combined_png = os.path.join(
        FIG_DIR,
        "combined_target_gene_expression_tau_enhancer_phyloP_phastCons_GERP_EnhA_remove_outliers.violin.png"
    )

    combined_svg = os.path.join(
        FIG_DIR,
        "combined_target_gene_expression_tau_enhancer_phyloP_phastCons_GERP_EnhA_remove_outliers.violin.svg"
    )

    plt.savefig(combined_pdf, bbox_inches="tight")
    plt.savefig(combined_png, dpi=600, bbox_inches="tight")
    plt.savefig(combined_svg, bbox_inches="tight")
    plt.close()

    print("[INFO] Combined figure:")
    print(combined_pdf)
    print(combined_png)
    print(combined_svg)


# ============================================================
# ============================================================

def resolve_chrom_name(chrom, bw_chroms):
    """
      chr1 vs 1  .
    """
    chrom = str(chrom)

    if chrom in bw_chroms:
        return chrom

    if chrom.startswith("chr"):
        alt = chrom.replace("chr", "", 1)
        if alt in bw_chroms:
            return alt
    else:
        alt = "chr" + chrom
        if alt in bw_chroms:
            return alt

    return None


def interval_bw_stats(bw, bw_chroms, chrom, start, end):
    """
      enhancer   mean,coverage,length.
    """
    bw_chrom = resolve_chrom_name(chrom, bw_chroms)

    if bw_chrom is None:
        return np.nan, 0.0, 0

    chrom_len = bw_chroms[bw_chrom]

    start = max(0, int(start))
    end = min(int(end), int(chrom_len))

    if end <= start:
        return np.nan, 0.0, 0

    try:
        try:
            mean_val = bw.stats(
                bw_chrom,
                start,
                end,
                type="mean",
                nBins=1,
                exact=True
            )[0]

            cov_val = bw.stats(
                bw_chrom,
                start,
                end,
                type="coverage",
                nBins=1,
                exact=True
            )[0]
        except TypeError:
            mean_val = bw.stats(
                bw_chrom,
                start,
                end,
                type="mean",
                nBins=1
            )[0]

            cov_val = bw.stats(
                bw_chrom,
                start,
                end,
                type="coverage",
                nBins=1
            )[0]

    except RuntimeError:
        return np.nan, 0.0, 0

    if mean_val is None:
        mean_val = np.nan

    if cov_val is None:
        cov_val = 0.0

    length = end - start

    return mean_val, cov_val, length


def enhancer_bw_mean(enhancer_regions, bw_path, score_name):
    """
      enhancer region   bigWig  .
    statistics  enhancer  ,  gene body.
    """
    print(f"[INFO] Reading bigWig for enhancer regions: {bw_path}")

    bw = pyBigWig.open(bw_path)
    bw_chroms = bw.chroms()

    print(f"[INFO] First 10 chroms in {score_name}: {list(bw_chroms.keys())[:10]}")

    rows = []

    for row in enhancer_regions.itertuples(index=False):
        tissue_count_group = getattr(row, "tissue_count_group")
        group = getattr(row, "group")
        enhancer_id = getattr(row, "enhancer_id")
        chrom = getattr(row, "chr")
        start = getattr(row, "start")
        end = getattr(row, "end")

        mean_val, cov_frac, length = interval_bw_stats(
            bw=bw,
            bw_chroms=bw_chroms,
            chrom=chrom,
            start=start,
            end=end
        )

        if length > 0 and cov_frac > 0:
            covered_bp = cov_frac * length
        else:
            covered_bp = 0.0

        rows.append({
            "tissue_count_group": tissue_count_group,
            "group": group,
            "group_full_name": GROUP_FULL_NAME[group],
            "tissue_range": GROUP_TISSUE_RANGE[group],
            "enhancer_id": enhancer_id,
            "chr": chrom,
            "start": int(start),
            "end": int(end),
            f"{score_name}_mean": mean_val,
            f"{score_name}_coverage": cov_frac,
            f"{score_name}_covered_bp": covered_bp,
            f"{score_name}_region_bp": length,
        })

    bw.close()

    return pd.DataFrame(rows)


# ============================================================
# ============================================================

pair = read_table_auto(PAIR_FILE, header=0)

required_pair_cols = {
    "tissue_count_group",
    "tissue_count_range",
    "n_tissues",
    "enhancer_id",
    "chr",
    "start",
    "end",
    "gene",
}

missing_cols = required_pair_cols - set(pair.columns)

if missing_cols:
    raise ValueError(f"PAIR_FILE Missing required columns: {missing_cols}")

pair["gene"] = pair["gene"].astype(str).str.strip()
pair["enhancer_id"] = pair["enhancer_id"].astype(str).str.strip()
pair["chr"] = pair["chr"].astype(str).str.strip()

pair["start"] = pd.to_numeric(pair["start"], errors="coerce")
pair["end"] = pd.to_numeric(pair["end"], errors="coerce")

pair = pair.dropna(subset=["chr", "start", "end"]).copy()
pair["start"] = pair["start"].astype(int)
pair["end"] = pair["end"].astype(int)

pair = pair[pair["end"] > pair["start"]].copy()

pair = pair[
    (~pair["gene"].isin(["", "NA", "nan", "None"]))
    & pair["gene"].notna()
].copy()

pair = pair[pair["tissue_count_group"].isin(GROUP_ORDER)].copy()
pair["group"] = pair["tissue_count_group"].map(GROUP_LABEL)

bad_enhancer_id = pair["enhancer_id"].isin(["", "NA", "nan", "None"])

pair.loc[bad_enhancer_id, "enhancer_id"] = (
    pair.loc[bad_enhancer_id, "chr"].astype(str)
    + ":"
    + pair.loc[bad_enhancer_id, "start"].astype(str)
    + "-"
    + pair.loc[bad_enhancer_id, "end"].astype(str)
)

pair["group"] = pd.Categorical(pair["group"], categories=PLOT_GROUPS, ordered=True)
pair["group_full_name"] = pair["group"].astype(str).map(GROUP_FULL_NAME)
pair["tissue_range"] = pair["group"].astype(str).map(GROUP_TISSUE_RANGE)


# ============================================================
# ============================================================

group_label_table = pd.DataFrame({
    "original_group": GROUP_ORDER,
    "EnhA_group": [GROUP_LABEL[g] for g in GROUP_ORDER],
    "full_name": [GROUP_FULL_NAME[GROUP_LABEL[g]] for g in GROUP_ORDER],
    "tissue_range": [GROUP_TISSUE_RANGE[GROUP_LABEL[g]] for g in GROUP_ORDER],
})

group_label_file = os.path.join(
    TABLE_DIR,
    "EnhA_group_label_mapping.tsv"
)

group_label_table.to_csv(
    group_label_file,
    sep="\t",
    index=False
)


# ============================================================
# ============================================================

group_gene = (
    pair[[
        "tissue_count_group",
        "group",
        "group_full_name",
        "tissue_range",
        "gene"
    ]]
    .drop_duplicates()
    .copy()
)

group_gene["group"] = pd.Categorical(
    group_gene["group"],
    categories=PLOT_GROUPS,
    ordered=True
)


# ============================================================
# ============================================================

group_enhancer = (
    pair[[
        "tissue_count_group",
        "group",
        "group_full_name",
        "tissue_range",
        "enhancer_id",
        "chr",
        "start",
        "end"
    ]]
    .drop_duplicates()
    .copy()
)

group_enhancer["group"] = pd.Categorical(
    group_enhancer["group"],
    categories=PLOT_GROUPS,
    ordered=True
)

group_enhancer = group_enhancer.sort_values(
    ["group", "chr", "start", "end", "enhancer_id"]
).copy()

print(f"[INFO] group-gene records: {group_gene.shape[0]}")
print(f"[INFO] unique genes in 5 groups: {group_gene['gene'].nunique()}")
print(f"[INFO] group-enhancer records: {group_enhancer.shape[0]}")
print(f"[INFO] unique enhancer regions in 5 groups: {group_enhancer['enhancer_id'].nunique()}")


# ============================================================
# ============================================================

expr = read_table_auto(EXPR_FILE, header=0)
expr.columns = [str(c).strip() for c in expr.columns]

if "Gene" in expr.columns:
    expr_gene_col = "Gene"
elif "ID" in expr.columns:
    expr_gene_col = "ID"
elif "gene" in expr.columns:
    expr_gene_col = "gene"
elif "id" in expr.columns:
    expr_gene_col = "id"
else:
    expr_gene_col = expr.columns[0]

expr["gene"] = expr[expr_gene_col].astype(str).str.strip()

lower_to_col = {str(c).lower(): c for c in expr.columns}

if "end" in lower_to_col:
    end_col = lower_to_col["end"]
    end_idx = list(expr.columns).index(end_col)
    sample_cols = list(expr.columns[(end_idx + 1):])
else:
    first6_lower = [str(c).lower() for c in expr.columns[:6]]
    metadata_words = {"chr", "chrom", "chromosome", "start", "end", "gene", "id", "strand"}

    if len(expr.columns) > 6 and len(set(first6_lower) & metadata_words) >= 2:
        sample_cols = list(expr.columns[6:])
    else:
        candidate_cols = [c for c in expr.columns if c not in [expr_gene_col, "gene"]]
        sample_cols = []

        for c in candidate_cols:
            numeric_values = pd.to_numeric(expr[c], errors="coerce")
            if numeric_values.notna().sum() > 0:
                sample_cols.append(c)

if len(sample_cols) == 0:
    raise ValueError(" file , check EXPR_FILE  .")

expr_values = expr[sample_cols].apply(pd.to_numeric, errors="coerce")
expr["average_expression"] = expr_values.mean(axis=1, skipna=True)

expr_keep = (
    expr[["gene", "average_expression"]]
    .dropna(subset=["gene"])
    .sort_values(["gene", "average_expression"], ascending=[True, False])
    .drop_duplicates("gene", keep="first")
    .copy()
)

print(f"[INFO] genes in expression file: {expr_keep['gene'].nunique()}")


# ============================================================
# 13. read tau file
# ============================================================

tau_raw = read_table_auto(TAU_FILE, header=None)

if tau_raw.shape[1] < 2:
    raise ValueError("TAU_FILE  :gene   tau")

tau = tau_raw.iloc[:, :2].copy()
tau.columns = ["gene", "tau"]

tau["gene"] = tau["gene"].astype(str).str.strip()
tau["tau"] = pd.to_numeric(tau["tau"], errors="coerce")

tau = tau[~tau["gene"].str.lower().isin(["id", "gene", "genes", "symbol"])].copy()

tau_keep = (
    tau[["gene", "tau"]]
    .dropna(subset=["gene"])
    .sort_values(["gene", "tau"], ascending=[True, False])
    .drop_duplicates("gene", keep="first")
    .copy()
)

print(f"[INFO] genes in tau file: {tau_keep['gene'].nunique()}")


# ============================================================
# 14. enhancer region bigWig conservation
# ============================================================

phyloP_df = enhancer_bw_mean(
    enhancer_regions=group_enhancer,
    bw_path=PHYLOP_BW,
    score_name="enhancer_phyloP",
)

phastCons_df = enhancer_bw_mean(
    enhancer_regions=group_enhancer,
    bw_path=PHASTCONS_BW,
    score_name="enhancer_phastCons",
)

gerp_df = enhancer_bw_mean(
    enhancer_regions=group_enhancer,
    bw_path=GERP_BW,
    score_name="enhancer_GERP",
)

merge_keys = [
    "tissue_count_group",
    "group",
    "group_full_name",
    "tissue_range",
    "enhancer_id",
    "chr",
    "start",
    "end",
]

enhancer_bw_score = (
    phyloP_df
    .merge(phastCons_df, on=merge_keys, how="outer")
    .merge(gerp_df, on=merge_keys, how="outer")
)

enhancer_bw_score_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_enhancer_region_bigwig_conservation_scores.raw.tsv",
)

enhancer_bw_score.to_csv(
    enhancer_bw_score_file,
    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# ============================================================

gene_data = (
    group_gene
    .merge(expr_keep, on="gene", how="left")
    .merge(tau_keep, on="gene", how="left")
)

if USE_LOG10_EXPR:
    gene_data["expression_for_plot"] = np.log10(gene_data["average_expression"] + 1)
else:
    gene_data["expression_for_plot"] = gene_data["average_expression"]

gene_data["group"] = pd.Categorical(
    gene_data["group"],
    categories=PLOT_GROUPS,
    ordered=True
)


# ============================================================
# 16. enhancer-level conservation data
# ============================================================

enhancer_data_raw = enhancer_bw_score.copy()

enhancer_data_raw["group"] = pd.Categorical(
    enhancer_data_raw["group"],
    categories=PLOT_GROUPS,
    ordered=True
)


# ============================================================
# ============================================================

data_by_metric = {
    "expression_for_plot": gene_data.copy(),
    "tau": gene_data.copy(),
    "enhancer_phyloP_mean": enhancer_data_raw.copy(),
    "enhancer_phastCons_mean": enhancer_data_raw.copy(),
    "enhancer_GERP_mean": enhancer_data_raw.copy(),
}

outlier_setting_rows = []
removed_outlier_tables = []

for metric in METRICS_TO_PLOT:
    if metric in REMOVE_OUTLIER_METRICS:
        filtered_df, outlier_df, setting = remove_outliers_global_iqr(
            df=data_by_metric[metric],
            value_col=metric,
            metric=metric
        )

        data_by_metric[metric] = filtered_df
        outlier_setting_rows.append(setting)

        if outlier_df.shape[0] > 0:
            removed_outlier_tables.append(outlier_df)

    else:
        valid_n = data_by_metric[metric][metric].notna().sum()

        outlier_setting_rows.append({
            "metric": metric,
            "method": "none",
            "iqr_multiplier": "NA",
            "q1": "NA",
            "q3": "NA",
            "iqr": "NA",
            "lower_bound": "NA",
            "upper_bound": "NA",
            "n_total_rows": data_by_metric[metric].shape[0],
            "n_valid_values": int(valid_n),
            "n_removed_outliers": 0,
            "n_kept_after_filter": int(valid_n),
        })

outlier_setting = pd.DataFrame(outlier_setting_rows)

outlier_setting_file = os.path.join(
    CHECK_DIR,
    "outlier_filter_settings_EnhA.tsv"
)

outlier_setting.to_csv(
    outlier_setting_file,
    sep="\t",
    index=False,
    na_rep="NA"
)

if len(removed_outlier_tables) > 0:
    removed_outliers = pd.concat(removed_outlier_tables, axis=0, ignore_index=True)
else:
    removed_outliers = pd.DataFrame()

removed_outlier_file = os.path.join(
    CHECK_DIR,
    "removed_outliers_phyloP_GERP_EnhA.tsv"
)

removed_outliers.to_csv(
    removed_outlier_file,
    sep="\t",
    index=False,
    na_rep="NA"
)


enhancer_values_filtered_long = []

for metric in ENHANCER_LEVEL_METRICS:
    tmp = data_by_metric[metric].copy()
    tmp["metric_for_plot_or_test"] = metric
    tmp["outlier_filter_applied"] = metric in REMOVE_OUTLIER_METRICS
    enhancer_values_filtered_long.append(tmp)

enhancer_values_filtered_long = pd.concat(
    enhancer_values_filtered_long,
    axis=0,
    ignore_index=True
)

enhancer_values_filtered_long_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_enhancer_region_conservation.filtered_by_metric.long.tsv"
)

enhancer_values_filtered_long.to_csv(
    enhancer_values_filtered_long_file,
    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# ============================================================

summary_list = []
pairwise_list = []

for metric in METRICS_TO_PLOT:
    metric_df = data_by_metric[metric]

    metric_summary, metric_pairwise = make_summary(
        metric_df,
        metric,
        metric
    )

    summary_list.append(metric_summary)

    if metric_pairwise is not None and metric_pairwise.shape[0] > 0:
        pairwise_list.append(metric_pairwise)

summary_all = pd.concat(summary_list, axis=0, ignore_index=True)

if len(pairwise_list) > 0:
    pairwise_all = pd.concat(pairwise_list, axis=0, ignore_index=True)
else:
    pairwise_all = pd.DataFrame(
        columns=["metric", "level", "group1", "group2", "p_raw", "p_adj"]
    )


# ============================================================
# ============================================================

plot_check_rows = []

for metric in METRICS_TO_PLOT:
    metric_df = data_by_metric[metric]

    for g in PLOT_GROUPS:
        n_valid = metric_df.loc[metric_df["group"] == g, metric].dropna().shape[0]

        plot_check_rows.append({
            "metric": metric,
            "level": METRIC_LEVEL[metric],
            "group": g,
            "group_full_name": GROUP_FULL_NAME[g],
            "tissue_range": GROUP_TISSUE_RANGE[g],
            "n_valid_after_outlier_filter": n_valid,
            "outlier_filter_applied": metric in REMOVE_OUTLIER_METRICS,
        })

plot_check = pd.DataFrame(plot_check_rows)

plot_check_file = os.path.join(
    CHECK_DIR,
    "check_valid_values_before_plot_after_outlier_filter_EnhA.tsv"
)

plot_check.to_csv(plot_check_file, sep="\t", index=False)

print("[INFO] valid values before plot after outlier filtering:")
print(plot_check)


# ============================================================
# ============================================================

for metric in METRICS_TO_PLOT:
    metric_summary = summary_all[summary_all["metric"] == metric].copy()
    metric_df = data_by_metric[metric]

    plot_violin(
        df=metric_df,
        value_col=metric,
        y_label=METRIC_LABELS[metric],
        summary_df=metric_summary,
        out_prefix=METRIC_PREFIX[metric],
    )

plot_combined_violin_figure(
    data_by_metric=data_by_metric,
    summary_all=summary_all
)


# ============================================================
# ============================================================

gene_values_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_target_gene_expression_tau.values.tsv",
)

summary_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_gene_expression_tau_enhancer_region_conservation.remove_outliers.summary_with_letters.tsv",
)

pairwise_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_gene_expression_tau_enhancer_region_conservation.remove_outliers.pairwise_mannwhitney.tsv",
)

match_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_feature_match_summary.remove_outliers.tsv",
)

unmatched_gene_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_unmatched_target_genes_expression_tau.tsv",
)

unmatched_enhancer_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_unmatched_enhancer_regions_conservation.raw.tsv",
)

gene_group_match_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_gene_feature_match_by_group.tsv",
)

enhancer_group_match_file = os.path.join(
    TABLE_DIR,
    "five_groups_EnhA_enhancer_conservation_match_by_group.after_outlier_filter.tsv",
)

gene_data.to_csv(
    gene_values_file,
    sep="\t",
    index=False,
    na_rep="NA"
)

summary_all.to_csv(
    summary_file,
    sep="\t",
    index=False,
    na_rep="NA"
)

pairwise_all.to_csv(
    pairwise_file,
    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# ============================================================

match_summary = pd.DataFrame({
    "item": [
        "group_gene_records",
        "unique_genes_in_5groups",
        "genes_in_expression_file",
        "genes_in_tau_file",
        "matched_expression_gene_records",
        "matched_tau_gene_records",
        "unmatched_expression_gene_records",
        "unmatched_tau_gene_records",

        "group_enhancer_records_raw",
        "unique_enhancer_ids_in_5groups_raw",
        "matched_phyloP_enhancer_records_raw",
        "matched_phastCons_enhancer_records_raw",
        "matched_GERP_enhancer_records_raw",

        "matched_phyloP_enhancer_records_after_outlier_filter",
        "matched_phastCons_enhancer_records_after_outlier_filter",
        "matched_GERP_enhancer_records_after_outlier_filter",

        "removed_phyloP_outliers",
        "removed_GERP_outliers",
    ],
    "n": [
        group_gene.shape[0],
        group_gene["gene"].nunique(),
        expr_keep["gene"].nunique(),
        tau_keep["gene"].nunique(),
        gene_data["average_expression"].notna().sum(),
        gene_data["tau"].notna().sum(),
        gene_data["average_expression"].isna().sum(),
        gene_data["tau"].isna().sum(),

        group_enhancer.shape[0],
        group_enhancer["enhancer_id"].nunique(),
        enhancer_data_raw["enhancer_phyloP_mean"].notna().sum(),
        enhancer_data_raw["enhancer_phastCons_mean"].notna().sum(),
        enhancer_data_raw["enhancer_GERP_mean"].notna().sum(),

        data_by_metric["enhancer_phyloP_mean"]["enhancer_phyloP_mean"].notna().sum(),
        data_by_metric["enhancer_phastCons_mean"]["enhancer_phastCons_mean"].notna().sum(),
        data_by_metric["enhancer_GERP_mean"]["enhancer_GERP_mean"].notna().sum(),

        int(
            outlier_setting.loc[
                outlier_setting["metric"] == "enhancer_phyloP_mean",
                "n_removed_outliers"
            ].iloc[0]
        ),
        int(
            outlier_setting.loc[
                outlier_setting["metric"] == "enhancer_GERP_mean",
                "n_removed_outliers"
            ].iloc[0]
        ),
    ],
})

match_summary.to_csv(match_file, sep="\t", index=False)


# ============================================================
# ============================================================

unmatched_gene = gene_data[
    gene_data["average_expression"].isna()
    | gene_data["tau"].isna()
].copy()

unmatched_gene.to_csv(
    unmatched_gene_file,
    sep="\t",
    index=False,
    na_rep="NA"
)

unmatched_enhancer = enhancer_data_raw[
    enhancer_data_raw["enhancer_phyloP_mean"].isna()
    | enhancer_data_raw["enhancer_phastCons_mean"].isna()
    | enhancer_data_raw["enhancer_GERP_mean"].isna()
].copy()

unmatched_enhancer.to_csv(
    unmatched_enhancer_file,
    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# ============================================================

gene_group_match = (
    gene_data
    .groupby("group", observed=False)
    .agg(
        total_group_gene_records=("gene", "count"),
        matched_expression=("average_expression", lambda x: x.notna().sum()),
        matched_tau=("tau", lambda x: x.notna().sum()),
        mean_expression=("average_expression", "mean"),
        mean_tau=("tau", "mean"),
    )
    .reset_index()
)

gene_group_match["group_full_name"] = gene_group_match["group"].astype(str).map(GROUP_FULL_NAME)
gene_group_match["tissue_range"] = gene_group_match["group"].astype(str).map(GROUP_TISSUE_RANGE)

gene_group_match["matched_expression_percent"] = (
    gene_group_match["matched_expression"]
    / gene_group_match["total_group_gene_records"]
    * 100
)

gene_group_match["matched_tau_percent"] = (
    gene_group_match["matched_tau"]
    / gene_group_match["total_group_gene_records"]
    * 100
)

gene_group_match = gene_group_match[
    [
        "group",
        "group_full_name",
        "tissue_range",
        "total_group_gene_records",
        "matched_expression",
        "matched_tau",
        "matched_expression_percent",
        "matched_tau_percent",
        "mean_expression",
        "mean_tau",
    ]
].copy()

gene_group_match.to_csv(
    gene_group_match_file,
    sep="\t",
    index=False,
    na_rep="NA"
)


enhancer_group_rows = []

for metric in ENHANCER_LEVEL_METRICS:
    tmp = data_by_metric[metric].copy()

    one = (
        tmp
        .groupby("group", observed=False)
        .agg(
            total_records_after_filter=("enhancer_id", "count"),
            matched_values=(metric, lambda x: x.notna().sum()),
        )
        .reset_index()
    )

    one["metric"] = metric
    one["group_full_name"] = one["group"].astype(str).map(GROUP_FULL_NAME)
    one["tissue_range"] = one["group"].astype(str).map(GROUP_TISSUE_RANGE)
    one["outlier_filter_applied"] = metric in REMOVE_OUTLIER_METRICS
    one["matched_percent_after_filter"] = (
        one["matched_values"] / one["total_records_after_filter"] * 100
    )

    enhancer_group_rows.append(one)

enhancer_group_match = pd.concat(enhancer_group_rows, axis=0, ignore_index=True)

enhancer_group_match = enhancer_group_match[
    [
        "metric",
        "group",
        "group_full_name",
        "tissue_range",
        "outlier_filter_applied",
        "total_records_after_filter",
        "matched_values",
        "matched_percent_after_filter",
    ]
].copy()

enhancer_group_match.to_csv(
    enhancer_group_match_file,
    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# ============================================================

print("")
print(" .")
print(f"Main output directory: {OUTDIR}")

print("")
print("EnhA grouping :")
print(group_label_file)

print("")
print(" :")
print(outlier_setting_file)

print("")
print(" :")
print(removed_outlier_file)

print("")
print(" filedirectory:")
print(FIG_DIR)

print("")
print(" file:")
for metric in METRICS_TO_PLOT:
    print(os.path.join(FIG_DIR, f"{METRIC_PREFIX[metric]}.violin.pdf"))

print("")
print("Combined figure:")
print(os.path.join(
    FIG_DIR,
    "combined_target_gene_expression_tau_enhancer_phyloP_phastCons_GERP_EnhA_remove_outliers.violin.pdf"
))
print(os.path.join(
    FIG_DIR,
    "combined_target_gene_expression_tau_enhancer_phyloP_phastCons_GERP_EnhA_remove_outliers.violin.png"
))
print(os.path.join(
    FIG_DIR,
    "combined_target_gene_expression_tau_enhancer_phyloP_phastCons_GERP_EnhA_remove_outliers.violin.svg"
))

print("")
print("results directory:")
print(TABLE_DIR)

print("")
print("results :")
print(gene_values_file)
print(enhancer_bw_score_file)
print(enhancer_values_filtered_long_file)
print(summary_file)
print(pairwise_file)
print(match_file)
print(unmatched_gene_file)
print(unmatched_enhancer_file)
print(gene_group_match_file)
print(enhancer_group_match_file)
print(plot_check_file)
