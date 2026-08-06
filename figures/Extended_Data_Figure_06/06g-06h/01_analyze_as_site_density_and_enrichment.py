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
from scipy.stats import fisher_exact


# ============================================================
# ============================================================

AS_DIR = "/vol2/zhangshiwen/GWAS/as_gwas/merge_as"

ENHANCER_FILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_enhancers.tsv"
)

OUTDIR = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "zutu/AS_overlap_5groups_EnhA_enrichment_only_combined"
)

os.makedirs(OUTDIR, exist_ok=True)

OUTPDF = os.path.join(
    OUTDIR,
    "combined_AS_overlap_5groups_EnhA.enrichment_only.square_heatmap.star_only.pdf"
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

GROUP_LABEL = {
    "G1_1_tissue": "TS-EnhA",
    "G2_2_5_tissues": "NS-EnhA",
    "G3_6_10_tissues": "MS-EnhA",
    "G4_11_20_tissues": "BS-EnhA",
    "G5_21_43_tissues": "ES-EnhA",
}

MARKER_ORDER = [
    "ATAC",
    "H3K27ac",
    "H3K4me1",
    "H3K4me3",
    "H3K27me3",
]

MARKER_COLORS = {
    "ATAC": "#FF83FA",
    "H3K27ac": "#FF0000",
    "H3K4me1": "#FFD700",
    "H3K4me3": "#32CD32",
    "H3K27me3": "#363636",
}


# ============================================================
# ============================================================

def check_bedtools():
    """
    check bedtools  .
    """
    bedtools = shutil.which("bedtools")

    if bedtools is None:
        raise RuntimeError(
            "  bedtools.  bedtools,"
            " :conda install -c bioconda bedtools"
        )

    return bedtools


def bh_adjust(pvalues):
    """
    Benjamini-Hochberg FDR correction.

    input:
        pvalues:   p

    output:
        BH   q
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
    """
     .
    """
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
    """
      chr/start/end  .
    """
    df = df.copy()
    df["_chr_order"] = df["chr"].map(chr_sort_key)
    df = df.sort_values(["_chr_order", "start", "end"])
    df = df.drop(columns=["_chr_order"])

    return df


def q_to_star(q):
    """
    FDR  .
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


def nature_style_ax(ax):
    """
    Nature-like axis style.
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
# 4. read AS file
# ============================================================

def parse_as_filename(filepath):
    """
    file :
    ATAC_abomasum_AS.tsv
    H3K27ac_abomasum_AS.tsv

     :
        marker, tissue
    """
    name = os.path.basename(filepath)

    m = re.match(r"^(.+?)_(.+)_AS\.tsv$", name)

    if m is None:
        raise ValueError(f"  AS file : {name}")

    marker = m.group(1)
    tissue = m.group(2)

    return marker, tissue


def empty_as_dataframe():
    """
      AS DataFrame.
    """
    return pd.DataFrame(
        columns=[
            "chr",
            "start",
            "end",
            "as_site_id",
            "marker",
            "tissue",
            "pos",
            "ref",
            "alt",
            "variant_id",
        ]
    )


def read_as_file(filepath):
    """
    read  AS file.

    input :
        chr pos ref alt variant_id

      1-bp BED:
        chr start end as_site_id marker tissue pos ref alt variant_id

     :
        start = pos - 1
        end   = pos
    """
    marker, tissue = parse_as_filename(filepath)
    empty_df = empty_as_dataframe()

    if (not os.path.exists(filepath)) or os.path.getsize(filepath) == 0:
        print(f"[WARNING] empty AS file skipped: {filepath}")
        return empty_df

    try:
        has_content = False

        with open(filepath, "r") as f:
            for line in f:
                if line.strip():
                    has_content = True
                    break

        if not has_content:
            print(f"[WARNING] blank AS file skipped: {filepath}")
            return empty_df

    except Exception as e:
        print(f"[WARNING] cannot check AS file, skipped: {filepath}; error={e}")
        return empty_df

    try:
        df = pd.read_csv(
            filepath,
            sep=r"\s+",
            header=None,
            dtype=str,
            engine="python",
        )
    except pd.errors.EmptyDataError:
        print(f"[WARNING] EmptyDataError AS file skipped: {filepath}")
        return empty_df
    except Exception as e:
        print(f"[WARNING] failed to read AS file, skipped: {filepath}; error={e}")
        return empty_df

    if df.shape[0] == 0:
        print(f"[WARNING] no rows in AS file skipped: {filepath}")
        return empty_df

    if df.shape[1] < 5:
        print(f"[WARNING] AS file has <5 columns, skipped: {filepath}")
        return empty_df

    df = df.iloc[:, :5].copy()
    df.columns = ["chr", "pos", "ref", "alt", "variant_id"]

    df["chr"] = df["chr"].map(normalize_chr)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")

    df = df.dropna(subset=["chr", "pos"]).copy()

    if df.shape[0] == 0:
        print(f"[WARNING] no valid AS positions after parsing, skipped: {filepath}")
        return empty_df

    df["pos"] = df["pos"].astype(int)
    df = df[df["pos"] > 0].copy()

    if df.shape[0] == 0:
        print(f"[WARNING] no positive AS positions, skipped: {filepath}")
        return empty_df

    df["start"] = df["pos"] - 1
    df["end"] = df["pos"]

    df["marker"] = marker
    df["tissue"] = tissue

    df["ref"] = df["ref"].astype(str)
    df["alt"] = df["alt"].astype(str)
    df["variant_id"] = df["variant_id"].astype(str)

    df["as_site_id"] = (
        df["marker"] + "|"
        + df["chr"] + ":"
        + df["pos"].astype(str) + ":"
        + df["ref"] + ">"
        + df["alt"]
    )

    df = df[
        [
            "chr",
            "start",
            "end",
            "as_site_id",
            "marker",
            "tissue",
            "pos",
            "ref",
            "alt",
            "variant_id",
        ]
    ].copy()

    return df


def read_all_as_files(as_dir):
    """
    read  *_AS.tsv file.

      marker  :
      marker   AS  , ,
      AS  .
    """
    files = sorted(glob.glob(os.path.join(as_dir, "*_AS.tsv")))

    if len(files) == 0:
        raise FileNotFoundError(f" directory  *_AS.tsv file: {as_dir}")

    all_list = []
    skipped_n = 0

    for f in files:
        marker, tissue = parse_as_filename(f)

        if marker not in MARKER_ORDER:
            print(f"[WARNING] skip unknown marker file: {f}")
            skipped_n += 1
            continue

        df = read_as_file(f)

        if df.shape[0] == 0:
            skipped_n += 1
            continue

        all_list.append(df)

        print(f"[INFO] read {os.path.basename(f)}: {df.shape[0]} AS sites")

    if len(all_list) == 0:
        raise ValueError(
            "  AS file file file,  overlap   AS  ."
        )

    all_as = pd.concat(all_list, axis=0, ignore_index=True)

    as_unique = (
        all_as
        .groupby(
            [
                "marker",
                "chr",
                "start",
                "end",
                "as_site_id",
                "pos",
                "ref",
                "alt",
                "variant_id",
            ],
            as_index=False
        )
        .agg(
            n_AS_tissues=("tissue", "nunique"),
            AS_tissue_list=("tissue", lambda x: ",".join(sorted(set(x.astype(str)))))
        )
    )

    as_unique = sort_bed_df(as_unique)

    print(f"[INFO] skipped empty/invalid/unknown AS files: {skipped_n}")

    return all_as, as_unique


# ============================================================
# 5. read enhancer groupingfile
# ============================================================

def read_enhancer_groups(enhancer_file):
    """
    read  enhancer file.

     :
        tissue_count_group, tissue_count_range, n_tissues,
        enhancer_id, chr, start, end
    """
    if not os.path.exists(enhancer_file):
        raise FileNotFoundError(f"Not found enhancer groupingfile: {enhancer_file}")

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
        raise ValueError(f"enhancer fileMissing required columns: {missing}")

    enh["chr"] = enh["chr"].map(normalize_chr)
    enh["start"] = pd.to_numeric(enh["start"], errors="coerce")
    enh["end"] = pd.to_numeric(enh["end"], errors="coerce")
    enh["n_tissues"] = pd.to_numeric(enh["n_tissues"], errors="coerce")

    enh = enh.dropna(
        subset=[
            "chr",
            "start",
            "end",
            "tissue_count_group",
        ]
    ).copy()

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
# ============================================================

def write_marker_as_beds(as_unique, tmpdir):
    """
    output  marker   AS BED  directory.

    BED  4 :
        chr start end as_site_id
    """
    marker_bed = {}

    for marker in MARKER_ORDER:
        sub = as_unique[as_unique["marker"] == marker].copy()

        if sub.shape[0] == 0:
            print(f"[WARNING] marker {marker}   AS  ")
            continue

        bed = sub[["chr", "start", "end", "as_site_id"]].copy()
        bed = sort_bed_df(bed)

        outfile = os.path.join(tmpdir, f"{marker}.AS_sites.bed")
        bed.to_csv(outfile, sep="\t", header=False, index=False)

        marker_bed[marker] = outfile

    return marker_bed


def write_group_enhancer_beds(enh, tmpdir):
    """
    output  enhancer group   BED  directory.

    BED  4 :
        chr start end enhancer_id
    """
    group_bed = {}

    for group in GROUP_ORDER:
        sub = enh[enh["tissue_count_group"] == group].copy()

        if sub.shape[0] == 0:
            print(f"[WARNING] group {group}   enhancer")
            continue

        bed = sub[["chr", "start", "end", "enhancer_id"]].copy()
        bed = sort_bed_df(bed)

        outfile = os.path.join(tmpdir, f"{group}.enhancers.bed")
        bed.to_csv(outfile, sep="\t", header=False, index=False)

        group_bed[group] = outfile

    return group_bed


# ============================================================
# 7. bedtools intersect
# ============================================================

def run_bedtools_intersect(marker, group, as_bed, enh_bed):
    """
      bedtools intersect -wa -wb   AS site   enhancer   overlap pair.

     output pair file,  DataFrame.
    """
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        as_bed,
        "-b",
        enh_bed,
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

    if result.stdout.strip() == "":
        pair = pd.DataFrame(
            columns=[
                "as_chr",
                "as_start",
                "as_end",
                "as_site_id",
                "enh_chr",
                "enh_start",
                "enh_end",
                "enhancer_id",
                "marker",
                "tissue_count_group",
                "tissue_count_range",
            ]
        )
    else:
        pair = pd.read_csv(
            io.StringIO(result.stdout),
            sep="\t",
            header=None,
            names=[
                "as_chr",
                "as_start",
                "as_end",
                "as_site_id",
                "enh_chr",
                "enh_start",
                "enh_end",
                "enhancer_id",
            ],
            dtype=str,
        )

        for c in ["as_start", "as_end", "enh_start", "enh_end"]:
            pair[c] = pd.to_numeric(pair[c], errors="coerce").astype("Int64")

        pair["marker"] = marker
        pair["tissue_count_group"] = group
        pair["tissue_count_range"] = GROUP_LABEL[group]

    return pair


# ============================================================
# ============================================================

def fisher_enhancer_enrichment(summary_df):
    """
      marker,  group   Fisher  .

     :
          group   AS-overlapped enhancer
          groups.

    contingency table:
                        with_AS     without_AS
        this group         a            b
        other groups       c            d

     :
        alternative="greater"

     :
          group   odds   other groups.
    """
    rows = []

    for marker in MARKER_ORDER:
        marker_df = summary_df[summary_df["marker"] == marker].copy()

        for group in GROUP_ORDER:
            this = marker_df[marker_df["tissue_count_group"] == group].copy()
            rest = marker_df[marker_df["tissue_count_group"] != group].copy()

            if this.shape[0] == 0 or rest.shape[0] == 0:
                rows.append({
                    "marker": marker,
                    "tissue_count_group": group,
                    "fisher_or": np.nan,
                    "log2_fisher_or": np.nan,
                    "fisher_p_enrichment": np.nan,
                })
                continue

            a = int(this["enhancers_with_AS"].sum())
            b = int(this["enhancer_total"].sum() - a)
            c = int(rest["enhancers_with_AS"].sum())
            d = int(rest["enhancer_total"].sum() - c)

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
                "marker": marker,
                "tissue_count_group": group,
                "fisher_or": fisher_or,
                "log2_fisher_or": log2_or,
                "fisher_p_enrichment": fisher_p,
                "a_this_group_with_AS": a,
                "b_this_group_without_AS": b,
                "c_other_groups_with_AS": c,
                "d_other_groups_without_AS": d,
            })

    fisher_df = pd.DataFrame(rows)

    fisher_df["fisher_q_BH_enrichment"] = bh_adjust(
        fisher_df["fisher_p_enrichment"].values
    )

    return fisher_df


# ============================================================
# 9. Combined figure
# ============================================================

def plot_combined(summary_df):
    """
     outputCombined figure.

    panel a:
        AS density per Mb enhancer

    panel b:
        enhancers with AS percent

    panel c:
        one-sided Fisher enrichment heatmap

    c  :
        1. heatmap   log2OR;
        2. log2OR > 0  ;
        3.   log2OR > 0   enrichment q < 0.05  ;
        4.  ;
        5. heatmap  .
    """
    fig = plt.figure(figsize=(13.8, 4.5))

    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.25, 1.05, 1.18],
        wspace=0.42
    )

    x = np.arange(len(GROUP_ORDER))

    # -------------------------
    # Panel a: density bar
    # -------------------------
    ax1 = fig.add_subplot(gs[0, 0])

    width = 0.15

    for i, marker in enumerate(MARKER_ORDER):
        sub = summary_df[summary_df["marker"] == marker].set_index("tissue_count_group")

        values = [
            sub.loc[g, "AS_site_density_per_Mb"] if g in sub.index else 0
            for g in GROUP_ORDER
        ]

        ax1.bar(
            x + (i - 2) * width,
            values,
            width=width,
            color=MARKER_COLORS[marker],
            label=marker,
            linewidth=0,
        )

    ax1.set_xticks(x)
    ax1.set_xticklabels(
        [GROUP_LABEL[g] for g in GROUP_ORDER],
        rotation=35,
        ha="right"
    )

    ax1.set_xlabel("EnhA tissue-breadth group", fontsize=10)
    ax1.set_ylabel("AS sites per Mb enhancer", fontsize=10)

    nature_style_ax(ax1)

    ax1.legend(
        frameon=False,
        fontsize=8,
        ncol=1
    )

    ax1.text(
        -0.15,
        1.05,
        "a",
        transform=ax1.transAxes,
        fontsize=13,
        fontweight="bold"
    )

    # -------------------------
    # Panel b: percent line
    # -------------------------
    ax2 = fig.add_subplot(gs[0, 1])

    for marker in MARKER_ORDER:
        sub = summary_df[summary_df["marker"] == marker].set_index("tissue_count_group")

        y = [
            sub.loc[g, "enhancers_with_AS_percent"] if g in sub.index else np.nan
            for g in GROUP_ORDER
        ]

        ax2.plot(
            x,
            y,
            marker="o",
            linewidth=2.0,
            markersize=5,
            color=MARKER_COLORS[marker],
            label=marker,
        )

    ax2.set_xticks(x)
    ax2.set_xticklabels(
        [GROUP_LABEL[g] for g in GROUP_ORDER],
        rotation=35,
        ha="right"
    )

    ax2.set_xlabel("EnhA tissue-breadth group", fontsize=10)
    ax2.set_ylabel("Enhancers with AS sites (%)", fontsize=10)

    nature_style_ax(ax2)

    ax2.text(
        -0.15,
        1.05,
        "b",
        transform=ax2.transAxes,
        fontsize=13,
        fontweight="bold"
    )

    # -------------------------
    # Panel c: one-sided enrichment heatmap
    # -------------------------
    ax3 = fig.add_subplot(gs[0, 2])

    mat = pd.DataFrame(
        index=MARKER_ORDER,
        columns=[GROUP_LABEL[g] for g in GROUP_ORDER],
        dtype=float
    )

    qmat = pd.DataFrame(
        index=MARKER_ORDER,
        columns=[GROUP_LABEL[g] for g in GROUP_ORDER],
        dtype=float
    )

    for _, row in summary_df.iterrows():
        marker = row["marker"]
        group_label = row["tissue_count_range"]

        mat.loc[marker, group_label] = row["log2_fisher_or"]
        qmat.loc[marker, group_label] = row["fisher_q_BH_enrichment"]

    mat_values = mat.values.astype(float)

    vmax = np.nanmax(np.abs(mat_values))

    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1

    norm = TwoSlopeNorm(
        vmin=-vmax,
        vcenter=0,
        vmax=vmax
    )

    im = ax3.imshow(
        mat_values,
        cmap="RdBu_r",
        norm=norm,
        aspect="equal",
        interpolation="nearest",
    )

    ax3.set_aspect("equal", adjustable="box")

    ax3.set_xticks(np.arange(len(mat.columns)))
    ax3.set_xticklabels(
        mat.columns,
        fontsize=9,
        rotation=35,
        ha="right"
    )

    ax3.set_yticks(np.arange(len(mat.index)))
    ax3.set_yticklabels(
        mat.index,
        fontsize=9
    )

    ax3.set_xlabel("EnhA tissue-breadth group", fontsize=10)
    ax3.set_ylabel("AS marker", fontsize=10)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            value = mat.iloc[i, j]
            q = qmat.iloc[i, j]

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

            ax3.text(
                j,
                i,
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

    for spine in ax3.spines.values():
        spine.set_visible(False)

    ax3.tick_params(length=0)

    cbar = fig.colorbar(
        im,
        ax=ax3,
        fraction=0.046,
        pad=0.04
    )

    cbar.set_label("log2 OR", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax3.text(
        -0.15,
        1.05,
        "c",
        transform=ax3.transAxes,
        fontsize=13,
        fontweight="bold"
    )

    plt.savefig(OUTPDF, bbox_inches="tight")
    plt.close()


# ============================================================
# 10. main workflow
# ============================================================

def main():
    check_bedtools()

    print("[INFO] Reading AS files...")
    all_as, as_unique = read_all_as_files(AS_DIR)

    print(f"[INFO] all AS per tissue rows: {all_as.shape[0]}")
    print(f"[INFO] unique AS by marker rows: {as_unique.shape[0]}")

    print("[INFO] Reading enhancer groups...")
    enh = read_enhancer_groups(ENHANCER_FILE)

    print(f"[INFO] enhancer rows: {enh.shape[0]}")

    summary_rows = []

    with tempfile.TemporaryDirectory(prefix="AS_overlap_EnhA_enrichment_only_") as tmpdir:
        print("[INFO] Writing temporary BED files...")
        marker_bed = write_marker_as_beds(as_unique, tmpdir)
        group_bed = write_group_enhancer_beds(enh, tmpdir)

        print("[INFO] Running bedtools intersect...")

        for marker in MARKER_ORDER:
            if marker not in marker_bed:
                print(f"[WARNING] No AS BED for marker {marker}, skipped.")
                continue

            marker_as = as_unique[as_unique["marker"] == marker].copy()
            marker_total_as = marker_as["as_site_id"].nunique()

            for group in GROUP_ORDER:
                if group not in group_bed:
                    print(f"[WARNING] No enhancer BED for group {group}, skipped.")
                    continue

                group_enh = enh[enh["tissue_count_group"] == group].copy()

                enhancer_total = group_enh["enhancer_id"].nunique()
                enhancer_bp = group_enh["enhancer_length"].sum()

                pair = run_bedtools_intersect(
                    marker=marker,
                    group=group,
                    as_bed=marker_bed[marker],
                    enh_bed=group_bed[group],
                )

                overlap_pair_count = pair.shape[0]

                as_overlap_sites = (
                    pair["as_site_id"].nunique()
                    if pair.shape[0] > 0 else 0
                )

                enhancers_with_as = (
                    pair["enhancer_id"].nunique()
                    if pair.shape[0] > 0 else 0
                )

                as_overlap_percent = (
                    as_overlap_sites / marker_total_as * 100
                    if marker_total_as > 0 else 0
                )

                enhancers_with_as_percent = (
                    enhancers_with_as / enhancer_total * 100
                    if enhancer_total > 0 else 0
                )

                as_density_per_mb = (
                    as_overlap_sites / enhancer_bp * 1e6
                    if enhancer_bp > 0 else np.nan
                )

                summary_rows.append({
                    "marker": marker,
                    "tissue_count_group": group,
                    "tissue_count_range": GROUP_LABEL[group],
                    "marker_total_AS_sites": marker_total_as,
                    "enhancer_total": enhancer_total,
                    "enhancer_total_bp": enhancer_bp,
                    "AS_overlap_sites": as_overlap_sites,
                    "AS_overlap_sites_percent_of_marker_AS": as_overlap_percent,
                    "AS_site_density_per_Mb": as_density_per_mb,
                    "enhancers_with_AS": enhancers_with_as,
                    "enhancers_with_AS_percent": enhancers_with_as_percent,
                    "overlap_pair_count": overlap_pair_count,
                })

                print(
                    f"[INFO] {marker} {GROUP_LABEL[group]}: "
                    f"AS_overlap_sites={as_overlap_sites}, "
                    f"enhancers_with_AS={enhancers_with_as}, "
                    f"density_per_Mb={as_density_per_mb:.3f}"
                )

    summary = pd.DataFrame(summary_rows)

    if summary.shape[0] == 0:
        raise ValueError("No summary results were generated; check the AS and enhancer files.")

    fisher_df = fisher_enhancer_enrichment(summary)

    summary = summary.merge(
        fisher_df,
        on=["marker", "tissue_count_group"],
        how="left"
    )

    summary["significance_enrichment"] = summary["fisher_q_BH_enrichment"].apply(q_to_star)

    print("[INFO] Plotting enrichment-only combined figure...")
    plot_combined(summary)

    print(" . output Combined figure.")
    print(f"Output directory: {OUTDIR}")
    print(f"Combined figure: {OUTPDF}")


if __name__ == "__main__":
    main()
