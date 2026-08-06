#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import re
import glob
import io
import shutil
import tempfile
import subprocess

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

from scipy.stats import fisher_exact, kruskal, mannwhitneyu


# ============================================================
# ============================================================

CHROMBPNET_ANNOT_DIR = "/data/home/sczd644/run/zsw_chrombpnet/snpscore/02_annotation"

ENHANCER_FILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_enhancers.tsv"
)

OUTDIR = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "zutu/ATAC_ChromBPNet_logfc_pval_5groups_final_combined_v3"
)

os.makedirs(OUTDIR, exist_ok=True)

OUTPDF = os.path.join(
    OUTDIR,
    "combined_ATAC_ChromBPNet_logfc_pval_5groups.final_v3.pdf"
)


# ============================================================
# ============================================================

LOGFC_PVAL_CUTOFF = 0.05
ALPHA = 0.05

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

GROUP_COLORS = {
    "G1_1_tissue": "#0072B2",
    "G2_2_5_tissues": "#E69F00",
    "G3_6_10_tissues": "#009E73",
    "G4_11_20_tissues": "#D55E00",
    "G5_21_43_tissues": "#CC79A7",
}

ATAC_COLOR = "#0072B2"

OUTLIER_IQR_MULTIPLIER = 1.5

LEFT_YLABEL_X = -0.14


# ============================================================
# 3. helper functions
# ============================================================

def check_bedtools():
    bedtools = shutil.which("bedtools")

    if bedtools is None:
        raise RuntimeError(
            "  bedtools.  bedtools, :conda install -c bioconda bedtools"
        )

    return bedtools


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
    ranked_p = p[order]

    adjusted = np.empty(n, dtype=float)
    prev = 1.0

    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked_p[i] * n / rank
        prev = min(prev, val)
        adjusted[i] = prev

    tmp = np.empty(n, dtype=float)
    tmp[order] = np.minimum(adjusted, 1.0)

    out[valid] = tmp

    return out


def compact_letter_display_ordered_by_mean(pmat, group_order_by_mean, alpha=0.05):
    """
      mean  .
    mean   a.
    """
    all_groups = list(pmat.index)
    letters = {g: "" for g in all_groups}
    letter_pool = list("abcdefghijklmnopqrstuvwxyz")

    for g in group_order_by_mean:
        placed = False

        for L in letter_pool:
            holders = [h for h in all_groups if L in letters[h]]

            if len(holders) == 0:
                letters[g] += L
                placed = True
                break

            can_share = True

            for h in holders:
                p = pmat.loc[g, h]

                if pd.isna(p) or p < alpha:
                    can_share = False
                    break

            if can_share:
                letters[g] += L
                placed = True
                break

        if not placed:
            raise RuntimeError(" ,  letter_pool.")

    return letters


def normalize_chr(x):
    """
      chr1 / chr2 / chrX.
    """
    x = str(x).strip()

    if x == "" or x.lower() == "nan":
        return x

    if x.startswith("chr"):
        return x

    return "chr" + x


def chr_sort_key(chrom):
    chrom = str(chrom).replace("chr", "")

    if chrom == "X":
        return 27

    if chrom == "Y":
        return 28

    if chrom in ["M", "MT"]:
        return 29

    try:
        return int(chrom)
    except ValueError:
        return 1000


def sort_bed_df(df):
    df = df.copy()
    df["_chr_order"] = df["chr"].map(chr_sort_key)
    df = df.sort_values(["_chr_order", "start", "end"])
    df = df.drop(columns=["_chr_order"])
    return df


def q_to_star(q):
    if pd.isna(q):
        return ""

    if q < 0.001:
        return "***"

    if q < 0.01:
        return "**"

    if q < 0.05:
        return "*"

    return ""


def nature_style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        labelsize=10,
        length=4,
        width=1.0
    )

    ax.grid(False)


def remove_outliers_global_iqr(df, value_col):
    """
      IQR  .
     .

    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    """
    work = df.copy()
    values = pd.to_numeric(work[value_col], errors="coerce")

    valid = values.notna()

    if valid.sum() == 0:
        return work.iloc[0:0].copy()

    q1 = values[valid].quantile(0.25)
    q3 = values[valid].quantile(0.75)
    iqr = q3 - q1

    if pd.isna(iqr) or iqr == 0:
        lower = values[valid].min()
        upper = values[valid].max()
    else:
        lower = q1 - OUTLIER_IQR_MULTIPLIER * iqr
        upper = q3 + OUTLIER_IQR_MULTIPLIER * iqr

    keep = valid & (values >= lower) & (values <= upper)

    return work.loc[keep].copy()


# ============================================================
# 4. read ChromBPNet ATAC annotation
# ============================================================

def parse_annotation_filename(filepath):
    """
    file :
    abomasum.annotations.tsv
    cerebral-cortex.annotations.tsv

      tissue.
    """
    name = os.path.basename(filepath)

    if not name.endswith(".annotations.tsv"):
        raise ValueError("  annotation file : {}".format(name))

    tissue = name.replace(".annotations.tsv", "")

    return tissue


def read_one_chrombpnet_annotation(filepath):
    """
    read  ChromBPNet annotation file.
      tissue,start,end,chrombpnet_record_id.
    """
    tissue = parse_annotation_filename(filepath)

    if (not os.path.exists(filepath)) or os.path.getsize(filepath) == 0:
        print("[WARNING] empty annotation skipped: {}".format(filepath))
        return pd.DataFrame()

    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            header=0,
            dtype=str,
            engine="python",
        )
    except pd.errors.EmptyDataError:
        print("[WARNING] EmptyDataError annotation skipped: {}".format(filepath))
        return pd.DataFrame()
    except Exception as e:
        print("[WARNING] failed to read annotation {}, error={}".format(filepath, e))
        return pd.DataFrame()

    if df.shape[0] == 0:
        return pd.DataFrame()

    required = {
        "chr",
        "pos",
        "allele1",
        "allele2",
        "variant_id",
        "logfc.pval",
    }

    missing = required - set(df.columns)

    if missing:
        print("[WARNING] annotation file missing columns {}, skipped: {}".format(missing, filepath))
        return pd.DataFrame()

    df["tissue"] = tissue

    df["chr"] = df["chr"].map(normalize_chr)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")

    df = df.dropna(subset=["chr", "pos", "variant_id"]).copy()

    if df.shape[0] == 0:
        return pd.DataFrame()

    df["pos"] = df["pos"].astype(int)
    df = df[df["pos"] > 0].copy()

    if df.shape[0] == 0:
        return pd.DataFrame()

    df["start"] = df["pos"] - 1
    df["end"] = df["pos"]

    df["variant_id"] = df["variant_id"].astype(str).str.strip()
    df["allele1"] = df["allele1"].astype(str).str.strip()
    df["allele2"] = df["allele2"].astype(str).str.strip()

    df["chrombpnet_record_id"] = (
        df["tissue"].astype(str) + "|" +
        df["variant_id"].astype(str)
    )

    numeric_cols = [
        "allele1_pred_counts",
        "allele2_pred_counts",
        "logfc",
        "abs_logfc",
        "jsd",
        "original_jsd",
        "logfc_x_jsd",
        "abs_logfc_x_jsd",
        "logfc.pval",
        "abs_logfc.pval",
        "jsd.pval",
        "logfc_x_jsd.pval",
        "abs_logfc_x_jsd.pval",
        "gene_distance_1",
        "gene_distance_2",
        "gene_distance_3",
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "abs_logfc" not in df.columns and "logfc" in df.columns:
        df["abs_logfc"] = df["logfc"].abs()

    df = (
        df.sort_values(
            ["tissue", "variant_id", "logfc.pval"],
            ascending=[True, True, True]
        )
        .drop_duplicates(subset=["tissue", "variant_id"], keep="first")
        .copy()
    )

    return df


def read_all_chrombpnet_annotations(annotation_dir):
    """
    read 02_annotation   *.annotations.tsv.
    """
    files = sorted(glob.glob(os.path.join(annotation_dir, "*.annotations.tsv")))

    if len(files) == 0:
        raise FileNotFoundError("  ChromBPNet annotation file: {}".format(annotation_dir))

    all_list = []

    for f in files:
        df = read_one_chrombpnet_annotation(f)

        if df.shape[0] > 0:
            all_list.append(df)
            print("[INFO] read annotation {}: {} rows".format(os.path.basename(f), df.shape[0]))

    if len(all_list) == 0:
        raise ValueError("  ChromBPNet annotation file .")

    anno = pd.concat(all_list, axis=0, ignore_index=True)
    anno = sort_bed_df(anno)

    return anno


# ============================================================
# 5. read enhancer grouping
# ============================================================

def read_enhancer_groups(enhancer_file):
    """
    read  enhancer file.
     :
    tissue_count_group, tissue_count_range, n_tissues, enhancer_id, chr, start, end
    """
    if not os.path.exists(enhancer_file):
        raise FileNotFoundError("Not found enhancer groupingfile: {}".format(enhancer_file))

    enh = pd.read_csv(enhancer_file, sep="\t", header=0, dtype=str)

    required_cols = {
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_id",
        "chr",
        "start",
        "end",
    }

    missing = required_cols - set(enh.columns)

    if missing:
        raise ValueError("enhancer fileMissing required columns: {}".format(missing))

    enh["chr"] = enh["chr"].map(normalize_chr)
    enh["start"] = pd.to_numeric(enh["start"], errors="coerce")
    enh["end"] = pd.to_numeric(enh["end"], errors="coerce")
    enh["n_tissues"] = pd.to_numeric(enh["n_tissues"], errors="coerce")

    enh = enh.dropna(subset=["chr", "start", "end", "tissue_count_group"]).copy()

    enh["start"] = enh["start"].astype(int)
    enh["end"] = enh["end"].astype(int)

    enh = enh[enh["end"] > enh["start"]].copy()
    enh = enh[enh["tissue_count_group"].isin(GROUP_ORDER)].copy()

    enh = enh.drop_duplicates(
        subset=[
            "tissue_count_group",
            "enhancer_id",
            "chr",
            "start",
            "end",
        ],
        keep="first"
    ).copy()

    enh["enhancer_length"] = enh["end"] - enh["start"]
    enh["group_label"] = enh["tissue_count_group"].map(GROUP_LABEL)

    enh = sort_bed_df(enh)

    return enh


# ============================================================
# 6. BED + bedtools overlap
# ============================================================

def write_chrombpnet_bed(anno, tmpdir):
    """
      ChromBPNet ATAC AS   BED.

    BED columns:
    chr start end chrombpnet_record_id tissue variant_id pos allele1 allele2
    """
    bed = anno[
        [
            "chr",
            "start",
            "end",
            "chrombpnet_record_id",
            "tissue",
            "variant_id",
            "pos",
            "allele1",
            "allele2",
        ]
    ].copy()

    bed = sort_bed_df(bed)

    outfile = os.path.join(tmpdir, "ATAC_ChromBPNet_AS_sites.bed")
    bed.to_csv(outfile, sep="\t", header=False, index=False)

    return outfile


def write_enhancer_bed(enh, tmpdir):
    """
      enhancer BED.

    BED columns:
    chr start end enhancer_id tissue_count_group tissue_count_range n_tissues enhancer_length
    """
    bed = enh[
        [
            "chr",
            "start",
            "end",
            "enhancer_id",
            "tissue_count_group",
            "tissue_count_range",
            "n_tissues",
            "enhancer_length",
        ]
    ].copy()

    bed = sort_bed_df(bed)

    outfile = os.path.join(tmpdir, "all_5groups.enhancers.bed")
    bed.to_csv(outfile, sep="\t", header=False, index=False)

    return outfile


def run_bedtools_intersect(chrombpnet_bed, enhancer_bed):
    """
    ChromBPNet ATAC AS   enhancer   overlap.
     results, output file.
    """
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        chrombpnet_bed,
        "-b",
        enhancer_bed,
        "-wa",
        "-wb",
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )

    names = [
        "as_chr",
        "as_start",
        "as_end",
        "chrombpnet_record_id",
        "tissue",
        "variant_id",
        "as_pos",
        "allele1",
        "allele2",
        "enh_chr",
        "enh_start",
        "enh_end",
        "enhancer_id",
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_length",
    ]

    if result.stdout.strip() == "":
        pair = pd.DataFrame(columns=names)
    else:
        pair = pd.read_csv(
            io.StringIO(result.stdout),
            sep="\t",
            header=None,
            names=names,
            dtype=str,
            engine="python",
        )

    numeric_cols = [
        "as_start",
        "as_end",
        "as_pos",
        "enh_start",
        "enh_end",
        "n_tissues",
        "enhancer_length",
    ]

    for col in numeric_cols:
        if col in pair.columns:
            pair[col] = pd.to_numeric(pair[col], errors="coerce")

    return pair


def merge_overlap_with_full_annotation(overlap_pairs, anno):
    """
    overlap results  chrombpnet_record_id   ChromBPNet annotation.
    """
    if overlap_pairs.shape[0] == 0:
        return overlap_pairs.copy()

    anno = anno.copy()

    merged = overlap_pairs.merge(
        anno,
        on="chrombpnet_record_id",
        how="left",
        suffixes=("", "_chrombpnet")
    )

    merged["has_chrombpnet_annotation"] = merged["logfc.pval"].notna()
    merged["logfc.pval"] = pd.to_numeric(merged["logfc.pval"], errors="coerce")

    if "logfc" in merged.columns:
        merged["logfc"] = pd.to_numeric(merged["logfc"], errors="coerce")

    if "abs_logfc" in merged.columns:
        merged["abs_logfc"] = pd.to_numeric(merged["abs_logfc"], errors="coerce")
    elif "logfc" in merged.columns:
        merged["abs_logfc"] = merged["logfc"].abs()

    merged["is_logfc_pval_lt_0p05"] = merged["logfc.pval"] < LOGFC_PVAL_CUTOFF

    return merged


# ============================================================
# 7. statistics summary
# ============================================================

def summarize_by_group(merged, enh):
    """
      EnhA breadth group statistics:
    - overlap   enhancer   ATAC AS
    - logfc.pval < 0.05
    -
    -   enhancer
    - abs_logfc
    """
    rows = []

    for group in GROUP_ORDER:
        sub_all = merged[merged["tissue_count_group"] == group].copy()
        group_enh = enh[enh["tissue_count_group"] == group].copy()

        enhancer_total = group_enh["enhancer_id"].nunique()

        enhancer_total_bp = (
            group_enh[["enhancer_id", "enhancer_length"]]
            .drop_duplicates()
            ["enhancer_length"]
            .sum()
        )

        sub_scored = sub_all[sub_all["has_chrombpnet_annotation"] == True].copy()

        event_scored = sub_scored.drop_duplicates(
            subset=["chrombpnet_record_id", "tissue_count_group"],
            keep="first"
        ).copy()

        event_sig = event_scored[event_scored["is_logfc_pval_lt_0p05"] == True].copy()

        overlap_events_with_chrombpnet = event_scored["chrombpnet_record_id"].nunique()
        sig_events = event_sig["chrombpnet_record_id"].nunique()

        sig_percent = (
            sig_events / overlap_events_with_chrombpnet * 100
            if overlap_events_with_chrombpnet > 0 else np.nan
        )

        sig_density_per_mb = (
            sig_events / enhancer_total_bp * 1e6
            if enhancer_total_bp > 0 else np.nan
        )

        enhancers_with_overlap = sub_all["enhancer_id"].nunique()

        enhancers_with_sig = (
            sub_all[sub_all["is_logfc_pval_lt_0p05"] == True]["enhancer_id"].nunique()
            if sub_all.shape[0] > 0 else 0
        )

        enhancers_with_sig_percent = (
            enhancers_with_sig / enhancer_total * 100
            if enhancer_total > 0 else np.nan
        )

        if "abs_logfc" in event_scored.columns and event_scored.shape[0] > 0:
            mean_abs_logfc_all = event_scored["abs_logfc"].mean()
            median_abs_logfc_all = event_scored["abs_logfc"].median()
        else:
            mean_abs_logfc_all = np.nan
            median_abs_logfc_all = np.nan

        if "abs_logfc" in event_sig.columns and event_sig.shape[0] > 0:
            mean_abs_logfc_sig = event_sig["abs_logfc"].mean()
            median_abs_logfc_sig = event_sig["abs_logfc"].median()
        else:
            mean_abs_logfc_sig = np.nan
            median_abs_logfc_sig = np.nan

        rows.append({
            "tissue_count_group": group,
            "tissue_count_range": GROUP_LABEL[group],
            "enhancer_total": enhancer_total,
            "enhancer_total_bp": enhancer_total_bp,
            "enhancers_with_overlap": enhancers_with_overlap,
            "overlap_events_with_ChromBPNet": overlap_events_with_chrombpnet,
            "logfc_pval_lt_0p05_events": sig_events,
            "logfc_pval_lt_0p05_percent": sig_percent,
            "logfc_pval_lt_0p05_density_per_Mb": sig_density_per_mb,
            "enhancers_with_logfc_pval_lt_0p05": enhancers_with_sig,
            "enhancers_with_logfc_pval_lt_0p05_percent": enhancers_with_sig_percent,
            "mean_abs_logfc_all_overlap_scored": mean_abs_logfc_all,
            "median_abs_logfc_all_overlap_scored": median_abs_logfc_all,
            "mean_abs_logfc_sig": mean_abs_logfc_sig,
            "median_abs_logfc_sig": median_abs_logfc_sig,
        })

    summary = pd.DataFrame(rows)

    return summary


def fisher_sig_enrichment(summary):
    """
      group   Fisher  .

      group   logfc.pval < 0.05   groups.

    2×2  :
                    significant     non-significant
    this group            a                b
    other groups          c                d

     :
    fisher_exact([[a, b], [c, d]], alternative="greater")
    """
    rows = []

    for group in GROUP_ORDER:
        this = summary[summary["tissue_count_group"] == group].copy()
        rest = summary[summary["tissue_count_group"] != group].copy()

        if this.shape[0] == 0 or rest.shape[0] == 0:
            rows.append({
                "tissue_count_group": group,
                "fisher_or": np.nan,
                "log2_fisher_or": np.nan,
                "fisher_p_enrichment": np.nan,
            })
            continue

        a = int(this["logfc_pval_lt_0p05_events"].sum())
        b = int(this["overlap_events_with_ChromBPNet"].sum() - a)
        c = int(rest["logfc_pval_lt_0p05_events"].sum())
        d = int(rest["overlap_events_with_ChromBPNet"].sum() - c)

        try:
            fisher_or, fisher_p = fisher_exact(
                [[a, b], [c, d]],
                alternative="greater"
            )
        except Exception:
            fisher_or, fisher_p = np.nan, np.nan

        log2_or = np.log2(
            ((a + 0.5) / (b + 0.5)) /
            ((c + 0.5) / (d + 0.5))
        )

        rows.append({
            "tissue_count_group": group,
            "fisher_or": fisher_or,
            "log2_fisher_or": log2_or,
            "fisher_p_enrichment": fisher_p,
            "a_this_group_significant": a,
            "b_this_group_non_significant": b,
            "c_other_groups_significant": c,
            "d_other_groups_non_significant": d,
        })

    fisher_df = pd.DataFrame(rows)

    fisher_df["fisher_q_BH_enrichment"] = bh_adjust(
        fisher_df["fisher_p_enrichment"].values
    )

    fisher_df["significance_enrichment"] = fisher_df[
        "fisher_q_BH_enrichment"
    ].apply(q_to_star)

    return fisher_df


# ============================================================
# ============================================================

def get_significant_abs_logfc_plot_df(merged):
    """
     :
    1. overlap   enhancer   ATAC-AS
    2.   ChromBPNet annotation
    3. logfc.pval < 0.05
    4. abs_logfc   NA

      tissue + variant_id + tissue_count_group  .
    """
    plot_df = merged[
        (merged["has_chrombpnet_annotation"] == True) &
        (merged["is_logfc_pval_lt_0p05"] == True) &
        (merged["abs_logfc"].notna())
    ].copy()

    plot_df = plot_df.drop_duplicates(
        subset=["tissue", "variant_id", "tissue_count_group"],
        keep="first"
    ).copy()

    return plot_df


def test_abs_logfc_letters_for_significant_sites(merged):
    """
      ATAC-AS   abs_logfc  :
    Kruskal-Wallis   +   Mann-Whitney U   + BH  +  .

     :
     .
     .
    """
    plot_df = get_significant_abs_logfc_plot_df(merged)

    arrays = []
    valid_groups = []

    for group in GROUP_ORDER:
        values = plot_df.loc[
            plot_df["tissue_count_group"] == group,
            "abs_logfc"
        ].dropna().astype(float).values

        if len(values) > 0:
            arrays.append(values)
            valid_groups.append(group)

    if len(valid_groups) >= 2:
        kruskal_p = kruskal(*arrays).pvalue
    else:
        kruskal_p = np.nan

    pairwise_rows = []

    for i in range(len(GROUP_ORDER)):
        for j in range(i + 1, len(GROUP_ORDER)):
            g1 = GROUP_ORDER[i]
            g2 = GROUP_ORDER[j]

            x = plot_df.loc[
                plot_df["tissue_count_group"] == g1,
                "abs_logfc"
            ].dropna().astype(float).values

            y = plot_df.loc[
                plot_df["tissue_count_group"] == g2,
                "abs_logfc"
            ].dropna().astype(float).values

            if len(x) == 0 or len(y) == 0:
                p = np.nan
            else:
                p = mannwhitneyu(
                    x,
                    y,
                    alternative="two-sided"
                ).pvalue

            pairwise_rows.append({
                "group1": g1,
                "group2": g2,
                "group1_label": GROUP_LABEL[g1],
                "group2_label": GROUP_LABEL[g2],
                "p_raw": p,
            })

    pairwise = pd.DataFrame(pairwise_rows)
    pairwise["p_adj_BH"] = bh_adjust(pairwise["p_raw"].values)

    pmat = pd.DataFrame(
        1.0,
        index=GROUP_ORDER,
        columns=GROUP_ORDER,
        dtype=float,
    )

    for _, row in pairwise.iterrows():
        if pd.notna(row["p_adj_BH"]):
            pmat.loc[row["group1"], row["group2"]] = row["p_adj_BH"]
            pmat.loc[row["group2"], row["group1"]] = row["p_adj_BH"]

    group_mean = {}

    for group in GROUP_ORDER:
        values = plot_df.loc[
            plot_df["tissue_count_group"] == group,
            "abs_logfc"
        ].dropna().astype(float)

        if values.shape[0] > 0:
            group_mean[group] = values.mean()
        else:
            group_mean[group] = -np.inf

    group_order_by_mean = sorted(
        GROUP_ORDER,
        key=lambda g: group_mean[g],
        reverse=True
    )

    letters_raw = compact_letter_display_ordered_by_mean(
        pmat=pmat,
        group_order_by_mean=group_order_by_mean,
        alpha=ALPHA
    )

    summary_rows = []

    for group in GROUP_ORDER:
        values = plot_df.loc[
            plot_df["tissue_count_group"] == group,
            "abs_logfc"
        ].dropna().astype(float)

        if values.shape[0] == 0:
            letter = "NA"
        else:
            letter = letters_raw.get(group, "")

        summary_rows.append({
            "tissue_count_group": group,
            "tissue_count_range": GROUP_LABEL[group],
            "n": int(values.shape[0]),
            "mean": values.mean() if values.shape[0] > 0 else np.nan,
            "median": values.median() if values.shape[0] > 0 else np.nan,
            "q25": values.quantile(0.25) if values.shape[0] > 0 else np.nan,
            "q75": values.quantile(0.75) if values.shape[0] > 0 else np.nan,
            "min": values.min() if values.shape[0] > 0 else np.nan,
            "max": values.max() if values.shape[0] > 0 else np.nan,
            "kruskal_p": kruskal_p,
            "letter": letter,
            "letter_order_basis": "mean_descending",
        })

    summary = pd.DataFrame(summary_rows)

    return summary, pairwise, plot_df


# ============================================================
# ============================================================

def draw_violin_panel(ax, plot_df, stat_summary):
    """
     .
      stat_summary  .
     .
    """
    violin_data = []
    violin_positions = []
    violin_groups = []

    box_data = []
    box_positions = []

    for i, group in enumerate(GROUP_ORDER, start=1):
        values = plot_df.loc[
            plot_df["tissue_count_group"] == group,
            "abs_logfc"
        ].dropna().astype(float).values

        if len(values) > 0:
            box_data.append(values)
            box_positions.append(i)

        if len(values) >= 2:
            violin_data.append(values)
            violin_positions.append(i)
            violin_groups.append(group)

    if len(violin_data) > 0:
        parts = ax.violinplot(
            violin_data,
            positions=violin_positions,
            widths=0.75,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

        for body, group in zip(parts["bodies"], violin_groups):
            body.set_facecolor(GROUP_COLORS[group])
            body.set_edgecolor("none")
            body.set_alpha(0.78)

    if len(box_data) > 0:
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

    ax.set_xticks(np.arange(1, len(GROUP_ORDER) + 1))
    ax.set_xticklabels(
        [GROUP_LABEL[g] for g in GROUP_ORDER],
        rotation=35,
        ha="right"
    )

    ax.set_xlabel("EnhA tissue-breadth group", fontsize=10)
    ax.set_ylabel("abs(logFC), significant ATAC AS sites", fontsize=10)

    nature_style_ax(ax)

    all_values = plot_df["abs_logfc"].dropna().astype(float)

    if all_values.shape[0] > 0:
        ymin = all_values.min()
        ymax = all_values.max()
        y_range = ymax - ymin if ymax > ymin else 1

        letter_y = ymax + 0.06 * y_range
        ax.set_ylim(ymin - 0.04 * y_range, ymax + 0.18 * y_range)

        letter_dict = dict(zip(stat_summary["tissue_count_group"], stat_summary["letter"]))

        for i, group in enumerate(GROUP_ORDER, start=1):
            letter = letter_dict.get(group, "NA")
            ax.text(
                i,
                letter_y,
                letter,
                ha="center",
                va="bottom",
                fontsize=11,
                color="black",
            )


def draw_vector_log2or_heatmap(ax_heat, mat, qmat, group_labels, x, norm, cmap, vmax):
    """
      heatmap,  imshow   PDF/AI  .
      heatmap   Rectangle,  Illustrator  .
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
        ["ATAC"],
        fontsize=9
    )

    ax_heat.tick_params(
        axis="y",
        left=False,
        right=False,
        labelleft=True,
        pad=7
    )

    ax_heat.set_ylabel(
        "AS marker",
        fontsize=10
    )

    ax_heat.yaxis.set_label_coords(LEFT_YLABEL_X, 0.5)

    for spine in ax_heat.spines.values():
        spine.set_visible(False)


def draw_vector_horizontal_colorbar(ax_cbar, cmap, norm, vmin, vmax, n_steps=256):
    """
      colorbar.
      colorbar   imshow  ,  Illustrator  .
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


def plot_combined_final(summary, merged):
    """
     Combined figure:
     :
         : ,  significant ATAC AS sites per Mb enhancer
         :  heatmap,  log2OR
         :  colorbar

     :
          heatmap   imshow;
          heatmap   Rectangle  ;
        colorbar   Rectangle  ;
        OUTPDF output file .
    """

    stat_summary_raw, pairwise_raw, plot_df_raw = test_abs_logfc_letters_for_significant_sites(merged)
    plot_df_no_outlier = remove_outliers_global_iqr(plot_df_raw, "abs_logfc")

    fig = plt.figure(figsize=(14.8, 5.2))

    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.35, 1.0, 1.0],
        wspace=0.40
    )

    group_labels = [GROUP_LABEL[g] for g in GROUP_ORDER]

    x = np.arange(1, len(GROUP_ORDER) + 1)

    sub = summary.set_index("tissue_count_group")

    # ========================================================
    # ========================================================

    left_bbox = gs[0, 0].get_position(fig)

    fig_w, fig_h = fig.get_size_inches()

    left_x0 = left_bbox.x0
    left_y0 = left_bbox.y0
    left_w = left_bbox.width
    left_y1 = left_bbox.y1

    heat_h = (left_w * fig_w / len(GROUP_ORDER)) / fig_h

    cbar_h = 0.035

    gap_heat_cbar = 0.145

    gap_line_heat = 0.018

    cbar_y0 = left_y0
    heat_y0 = cbar_y0 + cbar_h + gap_heat_cbar
    line_y0 = heat_y0 + heat_h + gap_line_heat
    line_h = left_y1 - line_y0

    if line_h <= 0:
        raise ValueError(
            " .  figsize  ,  gap_heat_cbar / heat_h."
        )

    ax_line = fig.add_axes(
        [left_x0, line_y0, left_w, line_h]
    )

    ax_heat = fig.add_axes(
        [left_x0, heat_y0, left_w, heat_h],
        sharex=ax_line
    )

    ax_cbar = fig.add_axes(
        [left_x0, cbar_y0, left_w, cbar_h]
    )

    # ========================================================
    # ========================================================

    y_density = [
        sub.loc[g, "logfc_pval_lt_0p05_density_per_Mb"] if g in sub.index else np.nan
        for g in GROUP_ORDER
    ]

    ax_line.plot(
        x,
        y_density,
        marker="o",
        linewidth=2.2,
        markersize=5,
        color=ATAC_COLOR,
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
        "Significant ATAC AS sites per Mb enhancer",
        fontsize=10
    )

    ax_line.yaxis.set_label_coords(LEFT_YLABEL_X, 0.5)

    nature_style_ax(ax_line)

    # ========================================================
    # ========================================================

    mat = pd.DataFrame(
        index=["ATAC"],
        columns=group_labels,
        dtype=float
    )

    qmat = pd.DataFrame(
        index=["ATAC"],
        columns=group_labels,
        dtype=float
    )

    for _, row in summary.iterrows():
        group_label = row["tissue_count_range"]
        mat.loc["ATAC", group_label] = row["log2_fisher_or"]
        qmat.loc["ATAC", group_label] = row["fisher_q_BH_enrichment"]

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

    # ========================================================
    # ========================================================

    ax_v1 = fig.add_subplot(gs[0, 1])

    draw_violin_panel(
        ax=ax_v1,
        plot_df=plot_df_raw,
        stat_summary=stat_summary_raw
    )

    # ========================================================
    # ========================================================

    ax_v2 = fig.add_subplot(gs[0, 2])

    draw_violin_panel(
        ax=ax_v2,
        plot_df=plot_df_no_outlier,
        stat_summary=stat_summary_raw
    )

    plt.savefig(
        OUTPDF,
        bbox_inches="tight",
        dpi=600,
        facecolor="white",
        transparent=False
    )

    plt.close()


# ============================================================
# 10. main workflow
# ============================================================

def main():
    check_bedtools()

    print("[INFO] Reading ChromBPNet ATAC annotation files...")
    anno = read_all_chrombpnet_annotations(CHROMBPNET_ANNOT_DIR)

    print("[INFO] ChromBPNet ATAC annotation rows:", anno.shape[0])
    print("[INFO] tissues:", anno["tissue"].nunique())

    print("[INFO] Reading enhancer groups...")
    enh = read_enhancer_groups(ENHANCER_FILE)

    print("[INFO] enhancer rows:", enh.shape[0])

    with tempfile.TemporaryDirectory(prefix="ATAC_ChromBPNet_final_v3_") as tmpdir:
        print("[INFO] Writing temporary BED files...")
        chrombpnet_bed = write_chrombpnet_bed(anno, tmpdir)
        enhancer_bed = write_enhancer_bed(enh, tmpdir)

        print("[INFO] Running bedtools intersect...")
        overlap_pairs = run_bedtools_intersect(chrombpnet_bed, enhancer_bed)

    print("[INFO] overlap pairs:", overlap_pairs.shape[0])

    print("[INFO] Merging overlap pairs with full ChromBPNet annotation...")
    merged = merge_overlap_with_full_annotation(overlap_pairs, anno)

    print("[INFO] merged rows:", merged.shape[0])

    print("[INFO] Summarizing...")
    summary = summarize_by_group(merged, enh)

    fisher_df = fisher_sig_enrichment(summary)

    summary = summary.merge(
        fisher_df,
        on="tissue_count_group",
        how="left"
    )

    print("[INFO] Plotting final combined figure only...")
    plot_combined_final(summary, merged)

    print("")
    print(" . output .")
    print("Output directory:", OUTDIR)
    print(" :", OUTPDF)


if __name__ == "__main__":
    main()
