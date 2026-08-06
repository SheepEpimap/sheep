#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import re
import glob
import shutil
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
from scipy.stats import kruskal, mannwhitneyu


# ============================================================
# ============================================================

EPI_AS_DIR = "/vol2/wulingyun/AS/AS_final_outcome/atac_cuttag_AS_fdr0.1"

# RNASeq ASE filedirectory
RNASEQ_ASE_DIR = "/vol2/wulingyun/AS/AS_final_outcome/RNASeq_ASE_fdr0.1"

ENHANCER_FILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_enhancers.tsv"
)

LINKED_PAIR_FILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_linked_enhancer_gene_pairs.tsv"
)

GENE_BED = "/vol2/mengzhu/genome/part_change_esemb100/Gene_esemble100_colin_ncbi.bed"

OUT_DIR = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "zutu/AS_ASE_effect_6panel_enhancer_and_target_gene"
)

COUNT_DIR = os.path.join(OUT_DIR, "01_merged_counts")
BED_DIR = os.path.join(OUT_DIR, "02_bed_for_overlap")
REGION_DIR = os.path.join(OUT_DIR, "03_regions_used")

TARGET_GENE_DIR = os.path.join(REGION_DIR, "target_gene_regions")

OVERLAP_DIR = os.path.join(OUT_DIR, "04_overlap")
VALUE_DIR = os.path.join(OUT_DIR, "05_values_for_violin")
SUMMARY_DIR = os.path.join(OUT_DIR, "06_statistics")
FIG_DIR = os.path.join(OUT_DIR, "07_figures")
TMP_DIR = os.path.join(OUT_DIR, "tmp")

for d in [
    OUT_DIR,
    COUNT_DIR,
    BED_DIR,
    REGION_DIR,
    TARGET_GENE_DIR,
    OVERLAP_DIR,
    VALUE_DIR,
    SUMMARY_DIR,
    FIG_DIR,
    TMP_DIR,
]:
    os.makedirs(d, exist_ok=True)


# ============================================================
# ============================================================

EPI_MARKERS_TO_USE = [
    "ATAC",
    "H3K27ac",
    "H3K4me1",
    "H3K4me3",
    "H3K27me3",
]

PANEL_ORDER = [
    "ATAC",
    "H3K27ac",
    "H3K4me1",
    "H3K4me3",
    "H3K27me3",
    "RNASeq_ASE",
]

PANEL_TITLE = {
    "ATAC": "ATAC",
    "H3K27ac": "H3K27ac",
    "H3K4me1": "H3K4me1",
    "H3K4me3": "H3K4me3",
    "H3K27me3": "H3K27me3",
    "RNASeq_ASE": "RNASeq ASE target genes",
}

REPS_TO_USE = ["39", "40"]

AS_PADJ_BH_CUTOFF = None

MIN_TOTAL_COUNT_MERGED = 0

#   abs_log2FC_alt_vs_ref_merged
#   log2FC_alt_vs_ref_merged
#   abs_alt_ratio_minus_0.5
EFFECT_COL = "abs_log2FC_alt_vs_ref_merged"

# event:
#
# unique:
ANALYSIS_LEVEL = "event"

ALPHA = 0.05

PROTEIN_CODING_ONLY = False

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


# ============================================================
# 3. helper functions
# ============================================================

def check_tools():
    if shutil.which("bedtools") is None:
        raise RuntimeError(
            "  bedtools.  bedtools, :conda install -c bioconda bedtools"
        )


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
    ranked = p[order]

    adj = np.empty(n, dtype=float)
    prev = 1.0

    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[i] = prev

    tmp = np.empty(n, dtype=float)
    tmp[order] = np.minimum(adj, 1.0)

    out[valid] = tmp

    return out


def normalize_chr(x):
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


def sort_bed_df(df, chr_col="chr", start_col="start", end_col="end"):
    df = df.copy()
    df["_chr_order"] = df[chr_col].map(chr_sort_key)
    df = df.sort_values(["_chr_order", start_col, end_col])
    df = df.drop(columns=["_chr_order"])
    return df


def safe_name(x):
    return str(x).replace("/", "_").replace(" ", "_").replace(":", "_")


def make_variant_id(contig, position, ref, alt):
    """
      1_1338992_C_G  .
    """
    chrom = str(contig)
    chrom = re.sub(r"^chr", "", chrom)
    return "{}_{}_{}_{}".format(chrom, int(position), str(ref), str(alt))


def nature_style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        labelsize=9,
        length=4,
        width=1.0
    )

    ax.grid(False)


def format_pvalue(p):
    if pd.isna(p):
        return "NA"

    if p < 1e-300:
        return "<1e-300"

    if p < 0.001:
        return "{:.1e}".format(p)

    return "{:.3f}".format(p)


# ============================================================
# ============================================================

def parse_epi_as_filename(path):
    """
     :

    ATAC:
        0.1_ATAC_abomasum_39.0.1_binomial_results.txt
        filtered_ATAC_abomasum_39_binomial_results.txt

    CUT&Tag:
        filtered_H3K27ac_abomasum_39_binomial_results.txt
        filtered_H3K4me1_abomasum_39_binomial_results.txt
        filtered_H3K4me3_abomasum_39_binomial_results.txt
        filtered_H3K27me3_abomasum_39_binomial_results.txt
    """
    name = os.path.basename(path)

    m1 = re.match(
        r"^0\.1_ATAC_(?P<tissue>.+)_(?P<rep>39|40)\.0\.1_binomial_results\.txt$",
        name
    )

    if m1:
        return {
            "marker": "ATAC",
            "tissue": m1.group("tissue"),
            "rep": m1.group("rep"),
        }

    m2 = re.match(
        r"^filtered_ATAC_(?P<tissue>.+)_(?P<rep>39|40)_binomial_results\.txt$",
        name
    )

    if m2:
        return {
            "marker": "ATAC",
            "tissue": m2.group("tissue"),
            "rep": m2.group("rep"),
        }

    m3 = re.match(
        r"^filtered_(?P<marker>H3K27ac|H3K27me3|H3K4me1|H3K4me3)_(?P<tissue>.+)_(?P<rep>39|40)_binomial_results\.txt$",
        name
    )

    if m3:
        return {
            "marker": m3.group("marker"),
            "tissue": m3.group("tissue"),
            "rep": m3.group("rep"),
        }

    return None


def parse_rnaseq_ase_filename(path):
    """
    RNASeq ASE file :
        filtered_RNASeq_abomasum_39_binomial_results.txt
    """
    name = os.path.basename(path)

    m = re.match(
        r"^filtered_RNASeq_(?P<tissue>.+)_(?P<rep>39|40)_binomial_results\.txt$",
        name
    )

    if m is None:
        return None

    return {
        "marker": "RNASeq_ASE",
        "tissue": m.group("tissue"),
        "rep": m.group("rep"),
    }


def discover_as_files(as_dir, parser_func, markers_to_use=None, out_prefix="AS"):
    files = sorted(glob.glob(os.path.join(as_dir, "*binomial_results.txt")))

    file_map = {}
    skipped = []

    for f in files:
        info = parser_func(f)

        if info is None:
            skipped.append(f)
            continue

        marker = info["marker"]
        tissue = info["tissue"]
        rep = info["rep"]

        if markers_to_use is not None and marker not in markers_to_use:
            continue

        file_map.setdefault(marker, {})
        file_map[marker].setdefault(tissue, {})
        file_map[marker][tissue][rep] = f

    skipped_file = os.path.join(
        OUT_DIR,
        "{}.skipped_unrecognized_files.txt".format(out_prefix)
    )

    with open(skipped_file, "w") as out:
        for f in skipped:
            out.write(f + "\n")

    rows = []

    for marker in sorted(file_map.keys()):
        for tissue in sorted(file_map[marker].keys()):
            reps = file_map[marker][tissue]

            rows.append({
                "marker": marker,
                "tissue": tissue,
                "has_39": "39" in reps,
                "has_40": "40" in reps,
                "file_39": reps.get("39", "NA"),
                "file_40": reps.get("40", "NA"),
            })

    discovered_df = pd.DataFrame(rows)

    discovered_file = os.path.join(
        OUT_DIR,
        "{}.discovered_marker_tissue_rep_files.tsv".format(out_prefix)
    )

    discovered_df.to_csv(discovered_file, sep="\t", index=False, na_rep="NA")

    return file_map, discovered_file


# ============================================================
# ============================================================

def read_one_as_file(path, source_type):
    """
    source_type:
        epi:
              AS file  end  .
        rnaseq:
            RNASeq ASE file  position.
              BED   start=position-1, end=position.
    """
    if path is None or (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t", header=0, dtype=str)
    df.columns = [str(c).strip() for c in df.columns]

    required = [
        "contig",
        "position",
        "ref",
        "alt",
        "refCount",
        "altCount",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError("fileMissing required columns {}: {}".format(missing, path))

    df["contig"] = df["contig"].map(normalize_chr)
    df["position"] = pd.to_numeric(df["position"], errors="coerce")

    if "end" in df.columns:
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
    else:
        if source_type == "rnaseq":
            df["end"] = df["position"]
        else:
            df["end"] = df["position"] + 1

    df["refCount"] = pd.to_numeric(df["refCount"], errors="coerce").fillna(0)
    df["altCount"] = pd.to_numeric(df["altCount"], errors="coerce").fillna(0)

    if "totalCount" in df.columns:
        df["totalCount"] = pd.to_numeric(df["totalCount"], errors="coerce")
        df["totalCount"] = df["totalCount"].fillna(df["refCount"] + df["altCount"])
    else:
        df["totalCount"] = df["refCount"] + df["altCount"]

    if "ratio" in df.columns:
        df["ratio"] = pd.to_numeric(df["ratio"], errors="coerce")
    else:
        df["ratio"] = np.nan

    if "p_value" in df.columns:
        df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    else:
        df["p_value"] = np.nan

    if "padj_bh" in df.columns:
        df["padj_bh"] = pd.to_numeric(df["padj_bh"], errors="coerce")
    else:
        df["padj_bh"] = np.nan

    df = df.dropna(subset=["contig", "position", "end", "ref", "alt"]).copy()

    df["position"] = df["position"].astype(int)
    df["end"] = df["end"].astype(int)
    df["ref"] = df["ref"].astype(str)
    df["alt"] = df["alt"].astype(str)

    df = df[df["position"] > 0].copy()

    if AS_PADJ_BH_CUTOFF is not None:
        df = df[df["padj_bh"] <= AS_PADJ_BH_CUTOFF].copy()

    if df.shape[0] == 0:
        return pd.DataFrame()

    group_cols = [
        "contig",
        "position",
        "end",
        "ref",
        "alt",
    ]

    out = (
        df.groupby(group_cols, as_index=False)
        .agg({
            "refCount": "sum",
            "altCount": "sum",
            "totalCount": "sum",
            "ratio": "mean",
            "p_value": "min",
            "padj_bh": "min",
        })
    )

    return out


def merge_rep_counts(marker, tissue, rep_files, source_type):
    """
      marker + tissue   39/40 replicate.
    """
    df39 = read_one_as_file(rep_files.get("39"), source_type=source_type)
    df40 = read_one_as_file(rep_files.get("40"), source_type=source_type)

    keys = [
        "contig",
        "position",
        "end",
        "ref",
        "alt",
    ]

    if df39.shape[0] == 0 and df40.shape[0] == 0:
        return pd.DataFrame()

    if df39.shape[0] > 0:
        df39 = df39.rename(columns={
            "refCount": "refCount_39",
            "altCount": "altCount_39",
            "totalCount": "totalCount_39",
            "ratio": "ratio_39",
            "p_value": "p_value_39",
            "padj_bh": "padj_bh_39",
        })
    else:
        df39 = pd.DataFrame(columns=keys + [
            "refCount_39",
            "altCount_39",
            "totalCount_39",
            "ratio_39",
            "p_value_39",
            "padj_bh_39",
        ])

    if df40.shape[0] > 0:
        df40 = df40.rename(columns={
            "refCount": "refCount_40",
            "altCount": "altCount_40",
            "totalCount": "totalCount_40",
            "ratio": "ratio_40",
            "p_value": "p_value_40",
            "padj_bh": "padj_bh_40",
        })
    else:
        df40 = pd.DataFrame(columns=keys + [
            "refCount_40",
            "altCount_40",
            "totalCount_40",
            "ratio_40",
            "p_value_40",
            "padj_bh_40",
        ])

    merged = pd.merge(df39, df40, on=keys, how="outer")

    count_cols = [
        "refCount_39",
        "altCount_39",
        "totalCount_39",
        "refCount_40",
        "altCount_40",
        "totalCount_40",
    ]

    for col in count_cols:
        if col not in merged.columns:
            merged[col] = 0

        merged[col] = pd.to_numeric(
            merged[col],
            errors="coerce"
        ).fillna(0)

    for col in [
        "ratio_39",
        "ratio_40",
        "p_value_39",
        "padj_bh_39",
        "p_value_40",
        "padj_bh_40",
    ]:
        if col not in merged.columns:
            merged[col] = np.nan

        merged[col] = pd.to_numeric(merged[col], errors="coerce")

    merged["marker"] = marker
    merged["tissue"] = tissue

    merged["refCount_merged"] = merged["refCount_39"] + merged["refCount_40"]
    merged["altCount_merged"] = merged["altCount_39"] + merged["altCount_40"]
    merged["totalCount_merged"] = merged["refCount_merged"] + merged["altCount_merged"]

    merged = merged[
        merged["totalCount_merged"] >= MIN_TOTAL_COUNT_MERGED
    ].copy()

    if merged.shape[0] == 0:
        return pd.DataFrame()

    merged["alt_ratio_merged"] = np.where(
        merged["totalCount_merged"] > 0,
        merged["altCount_merged"] / merged["totalCount_merged"],
        np.nan
    )

    merged["log2FC_alt_vs_ref_merged"] = np.log2(
        (merged["altCount_merged"] + 0.5) /
        (merged["refCount_merged"] + 0.5)
    )

    merged["abs_log2FC_alt_vs_ref_merged"] = (
        merged["log2FC_alt_vs_ref_merged"].abs()
    )

    merged["abs_alt_ratio_minus_0.5"] = (
        merged["alt_ratio_merged"] - 0.5
    ).abs()

    merged["p_value_min_reps"] = merged[
        ["p_value_39", "p_value_40"]
    ].min(axis=1, skipna=True)

    merged["padj_bh_min_reps"] = merged[
        ["padj_bh_39", "padj_bh_40"]
    ].min(axis=1, skipna=True)

    merged["variant_id"] = [
        make_variant_id(c, p, r, a)
        for c, p, r, a in zip(
            merged["contig"],
            merged["position"],
            merged["ref"],
            merged["alt"],
        )
    ]

    merged["variant_id_chr"] = (
        merged["contig"].astype(str) + "_" +
        merged["position"].astype(str) + "_" +
        merged["ref"].astype(str) + "_" +
        merged["alt"].astype(str)
    )

    merged["_chr_order"] = merged["contig"].map(chr_sort_key)
    merged = merged.sort_values(["_chr_order", "position", "ref", "alt"])
    merged = merged.drop(columns=["_chr_order"])

    out_cols = [
        "marker",
        "tissue",
        "contig",
        "position",
        "end",
        "ref",
        "alt",
        "variant_id",
        "variant_id_chr",
        "refCount_39",
        "altCount_39",
        "totalCount_39",
        "ratio_39",
        "refCount_40",
        "altCount_40",
        "totalCount_40",
        "ratio_40",
        "refCount_merged",
        "altCount_merged",
        "totalCount_merged",
        "alt_ratio_merged",
        "log2FC_alt_vs_ref_merged",
        "abs_log2FC_alt_vs_ref_merged",
        "abs_alt_ratio_minus_0.5",
        "p_value_39",
        "padj_bh_39",
        "p_value_40",
        "padj_bh_40",
        "p_value_min_reps",
        "padj_bh_min_reps",
    ]

    return merged[out_cols].copy()


# ============================================================
# ============================================================

def read_enhancer_groups(enhancer_file):
    enh = pd.read_csv(enhancer_file, sep="\t", header=0, dtype=str)

    required = [
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_id",
        "chr",
        "start",
        "end",
    ]

    missing = [c for c in required if c not in enh.columns]

    if missing:
        raise ValueError("enhancer fileMissing required columns: {}".format(missing))

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
            "chr",
            "start",
            "end",
        ],
        keep="first"
    ).copy()

    enh["enhancer_length"] = enh["end"] - enh["start"]
    enh["group_label"] = enh["tissue_count_group"].map(GROUP_LABEL)

    enh = sort_bed_df(enh, chr_col="chr", start_col="start", end_col="end")

    return enh


def write_enhancer_bed(enh):
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

    outfile = os.path.join(TMP_DIR, "all_5groups.enhancers.bed")
    bed.to_csv(outfile, sep="\t", header=False, index=False)

    return outfile


# ============================================================
# ============================================================

def read_gene_bed(gene_bed):
    """
    Gene_esemble100_colin_ncbi.bed  :
    chr1 start end LOC114110831 LOC114110831 - pseudogene

     4  gene_id  5  gene_name   gene.
    """
    gene = pd.read_csv(
        gene_bed,
        sep="\t",
        header=None,
        dtype=str,
        comment="#"
    )

    if gene.shape[1] < 4:
        raise ValueError("GENE_BED   4  : chr start end gene")

    cols = [
        "chr",
        "start",
        "end",
        "gene_id",
        "gene_name",
        "strand",
        "biotype",
    ]

    gene = gene.iloc[:, :min(gene.shape[1], len(cols))].copy()
    gene.columns = cols[:gene.shape[1]]

    if "gene_name" not in gene.columns:
        gene["gene_name"] = gene["gene_id"]

    if "strand" not in gene.columns:
        gene["strand"] = "NA"

    if "biotype" not in gene.columns:
        gene["biotype"] = "NA"

    gene["chr"] = gene["chr"].map(normalize_chr)
    gene["start"] = pd.to_numeric(gene["start"], errors="coerce")
    gene["end"] = pd.to_numeric(gene["end"], errors="coerce")

    gene = gene.dropna(
        subset=[
            "chr",
            "start",
            "end",
            "gene_id",
        ]
    ).copy()

    gene["start"] = gene["start"].astype(int)
    gene["end"] = gene["end"].astype(int)

    gene = gene[gene["end"] > gene["start"]].copy()

    if PROTEIN_CODING_ONLY:
        gene = gene[gene["biotype"] == "protein_coding"].copy()

    gene = gene.drop_duplicates(
        subset=[
            "chr",
            "start",
            "end",
            "gene_id",
            "gene_name",
        ],
        keep="first"
    ).copy()

    gene = sort_bed_df(gene, chr_col="chr", start_col="start", end_col="end")

    return gene


def read_target_gene_groups(linked_pair_file, gene_bed):
    """
      enhancer-gene pair file  EnhA group  ,
      Gene BED.
    """
    pair = pd.read_csv(linked_pair_file, sep="\t", header=0, dtype=str)

    required = [
        "tissue_count_group",
        "gene",
    ]

    missing = [c for c in required if c not in pair.columns]

    if missing:
        raise ValueError("linked enhancer-gene pair fileMissing required columns: {}".format(missing))

    pair = pair[pair["tissue_count_group"].isin(GROUP_ORDER)].copy()

    pair["gene"] = pair["gene"].replace(
        ["", "NA", "nan", "None", "null"],
        np.nan
    )

    pair = pair.dropna(subset=["gene"]).copy()

    target = (
        pair[["tissue_count_group", "gene"]]
        .drop_duplicates()
        .copy()
    )

    target["group_label"] = target["tissue_count_group"].map(GROUP_LABEL)

    gene = read_gene_bed(gene_bed)

    m1 = target.merge(
        gene,
        left_on="gene",
        right_on="gene_id",
        how="inner"
    )

    m2 = target.merge(
        gene,
        left_on="gene",
        right_on="gene_name",
        how="inner"
    )

    matched = pd.concat([m1, m2], axis=0, ignore_index=True)

    matched = matched.drop_duplicates(
        subset=[
            "tissue_count_group",
            "gene",
            "chr",
            "start",
            "end",
            "gene_id",
            "gene_name",
        ],
        keep="first"
    ).copy()

    matched["gene_length"] = matched["end"] - matched["start"]

    matched = sort_bed_df(
        matched,
        chr_col="chr",
        start_col="start",
        end_col="end"
    )

    unmatched = target[
        ~target["gene"].isin(matched["gene"].unique())
    ].copy()

    matched_file = os.path.join(
        TARGET_GENE_DIR,
        "target_genes_5groups.matched_gene_regions.tsv"
    )

    unmatched_file = os.path.join(
        TARGET_GENE_DIR,
        "target_genes_5groups.unmatched_genes.tsv"
    )

    target_file = os.path.join(
        TARGET_GENE_DIR,
        "target_genes_5groups.unique_gene_list.tsv"
    )

    target.to_csv(target_file, sep="\t", index=False, na_rep="NA")
    matched.to_csv(matched_file, sep="\t", index=False, na_rep="NA")
    unmatched.to_csv(unmatched_file, sep="\t", index=False, na_rep="NA")

    if matched.shape[0] == 0:
        raise ValueError(
            "  target gene   GENE_BED. check linked pair file  gene   Gene_esemble100_colin_ncbi.bed  4/ 5 ."
        )

    return matched, target_file, matched_file, unmatched_file


def write_target_gene_bed(target_gene_regions):
    """
    BED columns:
    chr start end gene gene_id gene_name tissue_count_group group_label strand biotype gene_length
    """
    bed = target_gene_regions[
        [
            "chr",
            "start",
            "end",
            "gene",
            "gene_id",
            "gene_name",
            "tissue_count_group",
            "group_label",
            "strand",
            "biotype",
            "gene_length",
        ]
    ].copy()

    outfile = os.path.join(TMP_DIR, "all_5groups.target_gene_regions.bed")
    bed.to_csv(outfile, sep="\t", header=False, index=False)

    return outfile


# ============================================================
# ============================================================

def write_as_bed(marker, marker_df, source_type):
    """
      AS:
          position/end.
    RNASeq ASE:
        position   1-based SNP position;
        BED   start = position - 1, end = position.
    """
    marker_df = marker_df.copy()

    if source_type == "rnaseq":
        as_start = marker_df["position"].astype(int) - 1
        as_end = marker_df["position"].astype(int)
    else:
        as_start = marker_df["position"].astype(int)
        as_end = marker_df["end"].astype(int)

    bed = pd.DataFrame({
        "as_chr": marker_df["contig"],
        "as_start": as_start,
        "as_end": as_end,
        "marker": marker_df["marker"],
        "tissue": marker_df["tissue"],
        "variant_id": marker_df["variant_id"],
        "variant_id_chr": marker_df["variant_id_chr"],
        "ref": marker_df["ref"],
        "alt": marker_df["alt"],
        "refCount_39": marker_df["refCount_39"],
        "altCount_39": marker_df["altCount_39"],
        "refCount_40": marker_df["refCount_40"],
        "altCount_40": marker_df["altCount_40"],
        "refCount_merged": marker_df["refCount_merged"],
        "altCount_merged": marker_df["altCount_merged"],
        "totalCount_merged": marker_df["totalCount_merged"],
        "alt_ratio_merged": marker_df["alt_ratio_merged"],
        "log2FC_alt_vs_ref_merged": marker_df["log2FC_alt_vs_ref_merged"],
        "abs_log2FC_alt_vs_ref_merged": marker_df["abs_log2FC_alt_vs_ref_merged"],
        "abs_alt_ratio_minus_0.5": marker_df["abs_alt_ratio_minus_0.5"],
        "p_value_min_reps": marker_df["p_value_min_reps"],
        "padj_bh_min_reps": marker_df["padj_bh_min_reps"],
    })

    bed = bed.dropna(subset=["as_chr", "as_start", "as_end"]).copy()

    bed["as_start"] = bed["as_start"].astype(int)
    bed["as_end"] = bed["as_end"].astype(int)

    bed = bed[bed["as_end"] > bed["as_start"]].copy()

    bed["_chr_order"] = bed["as_chr"].map(chr_sort_key)
    bed = bed.sort_values(["_chr_order", "as_start", "as_end"])
    bed = bed.drop(columns=["_chr_order"])

    outfile = os.path.join(
        BED_DIR,
        "{}.all_tissues.merged_counts.bed".format(safe_name(marker))
    )

    bed.to_csv(outfile, sep="\t", header=False, index=False)

    return outfile, bed


# ============================================================
# 9. bedtools overlap
# ============================================================

def run_epi_as_enhancer_overlap(marker, marker_bed, enhancer_bed):
    """
      AS   overlap   enhancer grouping.
    """
    raw_out = os.path.join(
        OVERLAP_DIR,
        "{}.AS_overlap_5groups_enhancer.raw.tsv".format(safe_name(marker))
    )

    cmd = [
        "bedtools",
        "intersect",
        "-a",
        marker_bed,
        "-b",
        enhancer_bed,
        "-wa",
        "-wb",
    ]

    with open(raw_out, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=True)

    a_cols = [
        "as_chr",
        "as_start",
        "as_end",
        "marker",
        "tissue",
        "variant_id",
        "variant_id_chr",
        "ref",
        "alt",
        "refCount_39",
        "altCount_39",
        "refCount_40",
        "altCount_40",
        "refCount_merged",
        "altCount_merged",
        "totalCount_merged",
        "alt_ratio_merged",
        "log2FC_alt_vs_ref_merged",
        "abs_log2FC_alt_vs_ref_merged",
        "abs_alt_ratio_minus_0.5",
        "p_value_min_reps",
        "padj_bh_min_reps",
    ]

    b_cols = [
        "enh_chr",
        "enh_start",
        "enh_end",
        "enhancer_id",
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_length",
    ]

    names = a_cols + b_cols

    if (not os.path.exists(raw_out)) or os.path.getsize(raw_out) == 0:
        overlap = pd.DataFrame(columns=names)
    else:
        overlap = pd.read_csv(raw_out, sep="\t", header=None, names=names, dtype=str)

    numeric_cols = [
        "as_start",
        "as_end",
        "refCount_39",
        "altCount_39",
        "refCount_40",
        "altCount_40",
        "refCount_merged",
        "altCount_merged",
        "totalCount_merged",
        "alt_ratio_merged",
        "log2FC_alt_vs_ref_merged",
        "abs_log2FC_alt_vs_ref_merged",
        "abs_alt_ratio_minus_0.5",
        "p_value_min_reps",
        "padj_bh_min_reps",
        "enh_start",
        "enh_end",
        "n_tissues",
        "enhancer_length",
    ]

    for col in numeric_cols:
        if col in overlap.columns:
            overlap[col] = pd.to_numeric(overlap[col], errors="coerce")

    overlap = overlap[overlap["tissue_count_group"].isin(GROUP_ORDER)].copy()

    full_out = os.path.join(
        OVERLAP_DIR,
        "{}.AS_overlap_5groups_enhancer.with_effect.tsv".format(safe_name(marker))
    )

    overlap.to_csv(full_out, sep="\t", index=False, na_rep="NA")

    return overlap, full_out, raw_out


def run_rnaseq_ase_target_gene_overlap(marker_bed, target_gene_bed):
    """
    RNASeq ASE   overlap   target gene regions.
    """
    marker = "RNASeq_ASE"

    raw_out = os.path.join(
        OVERLAP_DIR,
        "{}.ASE_overlap_5groups_target_genes.raw.tsv".format(marker)
    )

    cmd = [
        "bedtools",
        "intersect",
        "-a",
        marker_bed,
        "-b",
        target_gene_bed,
        "-wa",
        "-wb",
    ]

    with open(raw_out, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=True)

    a_cols = [
        "as_chr",
        "as_start",
        "as_end",
        "marker",
        "tissue",
        "variant_id",
        "variant_id_chr",
        "ref",
        "alt",
        "refCount_39",
        "altCount_39",
        "refCount_40",
        "altCount_40",
        "refCount_merged",
        "altCount_merged",
        "totalCount_merged",
        "alt_ratio_merged",
        "log2FC_alt_vs_ref_merged",
        "abs_log2FC_alt_vs_ref_merged",
        "abs_alt_ratio_minus_0.5",
        "p_value_min_reps",
        "padj_bh_min_reps",
    ]

    b_cols = [
        "gene_chr",
        "gene_start",
        "gene_end",
        "target_gene",
        "gene_id",
        "gene_name",
        "tissue_count_group",
        "group_label",
        "strand",
        "biotype",
        "gene_length",
    ]

    names = a_cols + b_cols

    if (not os.path.exists(raw_out)) or os.path.getsize(raw_out) == 0:
        overlap = pd.DataFrame(columns=names)
    else:
        overlap = pd.read_csv(raw_out, sep="\t", header=None, names=names, dtype=str)

    numeric_cols = [
        "as_start",
        "as_end",
        "refCount_39",
        "altCount_39",
        "refCount_40",
        "altCount_40",
        "refCount_merged",
        "altCount_merged",
        "totalCount_merged",
        "alt_ratio_merged",
        "log2FC_alt_vs_ref_merged",
        "abs_log2FC_alt_vs_ref_merged",
        "abs_alt_ratio_minus_0.5",
        "p_value_min_reps",
        "padj_bh_min_reps",
        "gene_start",
        "gene_end",
        "gene_length",
    ]

    for col in numeric_cols:
        if col in overlap.columns:
            overlap[col] = pd.to_numeric(overlap[col], errors="coerce")

    overlap = overlap[overlap["tissue_count_group"].isin(GROUP_ORDER)].copy()

    full_out = os.path.join(
        OVERLAP_DIR,
        "{}.ASE_overlap_5groups_target_genes.with_effect.tsv".format(marker)
    )

    overlap.to_csv(full_out, sep="\t", index=False, na_rep="NA")

    return overlap, full_out, raw_out


# ============================================================
# ============================================================

def build_epi_analysis_values(overlap, marker):
    """
      AS:overlap enhancer  ,  group   AS effect.
    """
    if overlap.shape[0] == 0:
        return pd.DataFrame()

    df = overlap.copy()
    df = df[df["tissue_count_group"].isin(GROUP_ORDER)].copy()
    df = df[df[EFFECT_COL].notna()].copy()

    if df.shape[0] == 0:
        return pd.DataFrame()

    if ANALYSIS_LEVEL == "event":
        df = (
            df.sort_values(
                ["marker", "tissue", "variant_id", "tissue_count_group", EFFECT_COL],
                ascending=[True, True, True, True, False]
            )
            .drop_duplicates(
                subset=["marker", "tissue", "variant_id", "tissue_count_group"],
                keep="first"
            )
            .copy()
        )

    elif ANALYSIS_LEVEL == "unique":
        df = (
            df.sort_values(
                ["marker", "variant_id", "tissue_count_group", EFFECT_COL],
                ascending=[True, True, True, False]
            )
            .drop_duplicates(
                subset=["marker", "variant_id", "tissue_count_group"],
                keep="first"
            )
            .copy()
        )

    else:
        raise ValueError("ANALYSIS_LEVEL   event   unique.")

    df["group_label"] = df["tissue_count_group"].map(GROUP_LABEL)
    df["effect_value"] = df[EFFECT_COL]
    df["region_type"] = "enhancer"
    df["target_gene"] = "NA"

    out_cols = [
        "marker",
        "region_type",
        "tissue",
        "variant_id",
        "variant_id_chr",
        "as_chr",
        "as_start",
        "as_end",
        "ref",
        "alt",
        "refCount_39",
        "altCount_39",
        "refCount_40",
        "altCount_40",
        "refCount_merged",
        "altCount_merged",
        "totalCount_merged",
        "alt_ratio_merged",
        "log2FC_alt_vs_ref_merged",
        "abs_log2FC_alt_vs_ref_merged",
        "abs_alt_ratio_minus_0.5",
        "effect_value",
        "p_value_min_reps",
        "padj_bh_min_reps",
        "enhancer_id",
        "enh_chr",
        "enh_start",
        "enh_end",
        "target_gene",
        "tissue_count_group",
        "tissue_count_range",
        "group_label",
        "n_tissues",
    ]

    out_cols = [c for c in out_cols if c in df.columns]

    df = df[out_cols].copy()

    value_file = os.path.join(
        VALUE_DIR,
        "{}.{}.{}.enhancer_overlap.values_for_violin.tsv".format(
            safe_name(marker),
            ANALYSIS_LEVEL,
            EFFECT_COL
        )
    )

    df.to_csv(value_file, sep="\t", index=False, na_rep="NA")

    return df


def build_rnaseq_analysis_values(overlap):
    """
    RNASeq ASE:overlap target gene  ,  group   ASE effect.
    """
    marker = "RNASeq_ASE"

    if overlap.shape[0] == 0:
        return pd.DataFrame()

    df = overlap.copy()
    df = df[df["tissue_count_group"].isin(GROUP_ORDER)].copy()
    df = df[df[EFFECT_COL].notna()].copy()

    if df.shape[0] == 0:
        return pd.DataFrame()

    if ANALYSIS_LEVEL == "event":
        df = (
            df.sort_values(
                ["marker", "tissue", "variant_id", "target_gene", "tissue_count_group", EFFECT_COL],
                ascending=[True, True, True, True, True, False]
            )
            .drop_duplicates(
                subset=["marker", "tissue", "variant_id", "target_gene", "tissue_count_group"],
                keep="first"
            )
            .copy()
        )

    elif ANALYSIS_LEVEL == "unique":
        df = (
            df.sort_values(
                ["marker", "variant_id", "target_gene", "tissue_count_group", EFFECT_COL],
                ascending=[True, True, True, True, False]
            )
            .drop_duplicates(
                subset=["marker", "variant_id", "target_gene", "tissue_count_group"],
                keep="first"
            )
            .copy()
        )

    else:
        raise ValueError("ANALYSIS_LEVEL   event   unique.")

    df["group_label"] = df["tissue_count_group"].map(GROUP_LABEL)
    df["effect_value"] = df[EFFECT_COL]
    df["region_type"] = "target_gene"
    df["tissue_count_range"] = df["group_label"]

    out_cols = [
        "marker",
        "region_type",
        "tissue",
        "variant_id",
        "variant_id_chr",
        "as_chr",
        "as_start",
        "as_end",
        "ref",
        "alt",
        "refCount_39",
        "altCount_39",
        "refCount_40",
        "altCount_40",
        "refCount_merged",
        "altCount_merged",
        "totalCount_merged",
        "alt_ratio_merged",
        "log2FC_alt_vs_ref_merged",
        "abs_log2FC_alt_vs_ref_merged",
        "abs_alt_ratio_minus_0.5",
        "effect_value",
        "p_value_min_reps",
        "padj_bh_min_reps",
        "gene_chr",
        "gene_start",
        "gene_end",
        "target_gene",
        "gene_id",
        "gene_name",
        "strand",
        "biotype",
        "tissue_count_group",
        "tissue_count_range",
        "group_label",
    ]

    out_cols = [c for c in out_cols if c in df.columns]

    df = df[out_cols].copy()

    value_file = os.path.join(
        VALUE_DIR,
        "{}.{}.{}.target_gene_overlap.values_for_violin.tsv".format(
            marker,
            ANALYSIS_LEVEL,
            EFFECT_COL
        )
    )

    df.to_csv(value_file, sep="\t", index=False, na_rep="NA")

    return df


# ============================================================
# ============================================================

def compact_letters_from_pairwise(pmat, means, alpha=0.05):
    groups = list(pmat.index)

    order = sorted(
        groups,
        key=lambda g: means.get(g, -np.inf) if pd.notna(means.get(g, np.nan)) else -np.inf,
        reverse=True
    )

    letters = {g: "" for g in groups}
    letter_pool = list("abcdefghijklmnopqrstuvwxyz")

    for g in order:
        placed = False

        for L in letter_pool:
            holders = [h for h in groups if L in letters[h]]

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
            letters[g] += "z"

    return letters


def summarize_and_test(values_df, marker):
    if values_df.shape[0] == 0:
        return pd.DataFrame(), pd.DataFrame()

    arrays = []
    valid_groups = []

    for group in GROUP_ORDER:
        vals = values_df.loc[
            values_df["tissue_count_group"] == group,
            "effect_value"
        ].dropna().astype(float).values

        if len(vals) > 0:
            arrays.append(vals)
            valid_groups.append(group)

    if len(valid_groups) >= 2:
        try:
            kruskal_p = kruskal(*arrays).pvalue
        except Exception:
            kruskal_p = np.nan
    else:
        kruskal_p = np.nan

    pair_rows = []

    for i in range(len(GROUP_ORDER)):
        for j in range(i + 1, len(GROUP_ORDER)):
            g1 = GROUP_ORDER[i]
            g2 = GROUP_ORDER[j]

            x = values_df.loc[
                values_df["tissue_count_group"] == g1,
                "effect_value"
            ].dropna().astype(float).values

            y = values_df.loc[
                values_df["tissue_count_group"] == g2,
                "effect_value"
            ].dropna().astype(float).values

            if len(x) == 0 or len(y) == 0:
                p = np.nan
            else:
                try:
                    p = mannwhitneyu(
                        x,
                        y,
                        alternative="two-sided"
                    ).pvalue
                except Exception:
                    p = np.nan

            pair_rows.append({
                "marker": marker,
                "group1": g1,
                "group2": g2,
                "group1_label": GROUP_LABEL[g1],
                "group2_label": GROUP_LABEL[g2],
                "p_raw": p,
            })

    pairwise = pd.DataFrame(pair_rows)
    pairwise["p_adj_BH"] = bh_adjust(pairwise["p_raw"].values)

    pmat = pd.DataFrame(
        1.0,
        index=GROUP_ORDER,
        columns=GROUP_ORDER,
        dtype=float
    )

    for _, row in pairwise.iterrows():
        if pd.notna(row["p_adj_BH"]):
            pmat.loc[row["group1"], row["group2"]] = row["p_adj_BH"]
            pmat.loc[row["group2"], row["group1"]] = row["p_adj_BH"]

    means = {}
    summary_rows = []

    for group in GROUP_ORDER:
        vals = values_df.loc[
            values_df["tissue_count_group"] == group,
            "effect_value"
        ].dropna().astype(float)

        means[group] = vals.mean() if vals.shape[0] > 0 else np.nan

    letters = compact_letters_from_pairwise(
        pmat,
        means,
        alpha=ALPHA
    )

    for group in GROUP_ORDER:
        vals = values_df.loc[
            values_df["tissue_count_group"] == group,
            "effect_value"
        ].dropna().astype(float)

        summary_rows.append({
            "marker": marker,
            "analysis_level": ANALYSIS_LEVEL,
            "effect_col": EFFECT_COL,
            "tissue_count_group": group,
            "tissue_count_range": GROUP_LABEL[group],
            "n": int(vals.shape[0]),
            "mean": vals.mean() if vals.shape[0] > 0 else np.nan,
            "median": vals.median() if vals.shape[0] > 0 else np.nan,
            "q25": vals.quantile(0.25) if vals.shape[0] > 0 else np.nan,
            "q75": vals.quantile(0.75) if vals.shape[0] > 0 else np.nan,
            "min": vals.min() if vals.shape[0] > 0 else np.nan,
            "max": vals.max() if vals.shape[0] > 0 else np.nan,
            "kruskal_p": kruskal_p,
            "letter": letters.get(group, "NA") if vals.shape[0] > 0 else "NA",
            "letter_order_basis": "mean_descending",
        })

    summary = pd.DataFrame(summary_rows)

    summary_file = os.path.join(
        SUMMARY_DIR,
        "{}.{}.{}.group_summary_with_letters.tsv".format(
            safe_name(marker),
            ANALYSIS_LEVEL,
            EFFECT_COL
        )
    )

    pair_file = os.path.join(
        SUMMARY_DIR,
        "{}.{}.{}.pairwise_mannwhitney_BH.tsv".format(
            safe_name(marker),
            ANALYSIS_LEVEL,
            EFFECT_COL
        )
    )

    summary.to_csv(summary_file, sep="\t", index=False, na_rep="NA")
    pairwise.to_csv(pair_file, sep="\t", index=False, na_rep="NA")

    return summary, pairwise


# ============================================================
# ============================================================

def draw_violin_on_ax(ax, values_df, summary_df, marker, panel_label=None):
    sub = values_df.copy()
    summ = summary_df.copy()

    violin_data = []
    violin_pos = []
    violin_groups = []

    box_data = []
    box_pos = []
    box_groups = []

    for i, group in enumerate(GROUP_ORDER, start=1):
        vals = sub.loc[
            sub["tissue_count_group"] == group,
            "effect_value"
        ].dropna().astype(float).values

        if len(vals) >= 2:
            violin_data.append(vals)
            violin_pos.append(i)
            violin_groups.append(group)

        if len(vals) > 0:
            box_data.append(vals)
            box_pos.append(i)
            box_groups.append(group)

    if len(violin_data) > 0:
        parts = ax.violinplot(
            violin_data,
            positions=violin_pos,
            widths=0.78,
            showmeans=False,
            showmedians=False,
            showextrema=False
        )

        for body, group in zip(parts["bodies"], violin_groups):
            body.set_facecolor(GROUP_COLORS[group])
            body.set_edgecolor("none")
            body.set_alpha(0.78)

    if len(box_data) > 0:
        ax.boxplot(
            box_data,
            positions=box_pos,
            widths=0.18,
            patch_artist=True,
            showfliers=False,
            showcaps=True,
            boxprops=dict(facecolor="white", edgecolor="black", linewidth=0.7),
            medianprops=dict(color="black", linewidth=0.8),
            whiskerprops=dict(color="black", linewidth=0.6),
            capprops=dict(color="black", linewidth=0.6),
        )

    for pos, vals, group in zip(box_pos, box_data, box_groups):
        if len(vals) == 1:
            ax.scatter(
                [pos],
                vals,
                s=18,
                color=GROUP_COLORS[group],
                edgecolor="black",
                linewidth=0.4,
                zorder=3
            )

    ax.set_xticks(np.arange(1, len(GROUP_ORDER) + 1))
    ax.set_xticklabels(
        [GROUP_LABEL[g] for g in GROUP_ORDER],
        rotation=35,
        ha="right"
    )

    ax.set_xlabel("EnhA tissue-breadth group", fontsize=9)

    if EFFECT_COL == "abs_log2FC_alt_vs_ref_merged":
        ax.set_ylabel("abs(log2FC alt/ref)", fontsize=9)
    elif EFFECT_COL == "log2FC_alt_vs_ref_merged":
        ax.set_ylabel("log2FC alt/ref", fontsize=9)
    elif EFFECT_COL == "abs_alt_ratio_minus_0.5":
        ax.set_ylabel("|alt ratio - 0.5|", fontsize=9)
    else:
        ax.set_ylabel(EFFECT_COL, fontsize=9)

    ax.set_title(PANEL_TITLE.get(marker, marker), fontsize=11)

    nature_style_ax(ax)

    vals_all = sub["effect_value"].dropna().astype(float)

    if vals_all.shape[0] > 0:
        ymin = vals_all.min()
        ymax = vals_all.max()
        yr = ymax - ymin if ymax > ymin else 1.0

        letter_y = ymax + 0.06 * yr
        ax.set_ylim(ymin - 0.04 * yr, ymax + 0.24 * yr)

        letter_dict = dict(zip(summ["tissue_count_group"], summ["letter"]))

        for i, group in enumerate(GROUP_ORDER, start=1):
            letter = letter_dict.get(group, "NA")

            if pd.notna(letter) and letter != "NA":
                ax.text(
                    i,
                    letter_y,
                    str(letter),
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="black"
                )

    if summ.shape[0] > 0 and "kruskal_p" in summ.columns:
        kp = summ["kruskal_p"].iloc[0]

        ax.text(
            0.03,
            0.97,
            "K-W P = {}".format(format_pvalue(kp)),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5
        )

    if panel_label is not None:
        ax.text(
            -0.18,
            1.08,
            panel_label,
            transform=ax.transAxes,
            fontsize=13,
            fontweight="bold"
        )


def plot_six_panel_figure(all_values, all_summary):
    """
    output 6 panel  :
    a-e:   AS-overlap enhancer effect
    f: RNASeq ASE-overlap target gene effect
    """
    markers_found = set(all_values["marker"].dropna().astype(str).unique())

    markers = [m for m in PANEL_ORDER if m in markers_found]

    if len(markers) == 0:
        print("[WARNING] no markers for combined figure.")
        return

    ncol = 3
    nrow = 2

    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(15.0, 8.2),
        squeeze=False
    )

    panel_letters = list("abcdef")

    for idx, marker in enumerate(markers):
        r = idx // ncol
        c = idx % ncol
        ax = axes[r][c]

        sub_values = all_values[all_values["marker"] == marker].copy()
        sub_summary = all_summary[all_summary["marker"] == marker].copy()

        draw_violin_on_ax(
            ax=ax,
            values_df=sub_values,
            summary_df=sub_summary,
            marker=marker,
            panel_label=panel_letters[idx]
        )

    for idx in range(len(markers), nrow * ncol):
        r = idx // ncol
        c = idx % ncol
        axes[r][c].axis("off")

    plt.tight_layout(w_pad=2.0, h_pad=2.2)

    prefix = "six_panels.epigenomic_AS_enhancer_and_RNASeq_ASE_target_gene.{}.{}".format(
        ANALYSIS_LEVEL,
        EFFECT_COL
    )

    pdf = os.path.join(FIG_DIR, prefix + ".pdf")
    png = os.path.join(FIG_DIR, prefix + ".png")
    svg = os.path.join(FIG_DIR, prefix + ".svg")

    plt.savefig(pdf, bbox_inches="tight")
    plt.savefig(png, dpi=600, bbox_inches="tight")
    plt.savefig(svg, bbox_inches="tight")
    plt.close()

    print("[INFO] six-panel figure:")
    print(pdf)
    print(png)
    print(svg)


# ============================================================
# 13. main workflow
# ============================================================

def process_epigenomic_as(enhancer_bed):
    """
      marker:
    AS   overlap enhancer  .
    """
    print("[INFO] Discovering epigenomic AS files...")

    file_map, discovered_file = discover_as_files(
        as_dir=EPI_AS_DIR,
        parser_func=parse_epi_as_filename,
        markers_to_use=EPI_MARKERS_TO_USE,
        out_prefix="epigenomic_AS"
    )

    print("[INFO] discovered epigenomic AS file table:", discovered_file)

    all_summaries = []
    all_pairwise = []
    all_values = []
    match_rows = []

    for marker in EPI_MARKERS_TO_USE:
        if marker not in file_map:
            print("[WARNING] marker not found:", marker)
            continue

        print("\n[INFO] Processing epigenomic marker:", marker)

        marker_count_list = []
        marker_count_dir = os.path.join(COUNT_DIR, safe_name(marker))
        os.makedirs(marker_count_dir, exist_ok=True)

        for tissue in sorted(file_map[marker].keys()):
            rep_files = file_map[marker][tissue]

            merged_counts = merge_rep_counts(
                marker=marker,
                tissue=tissue,
                rep_files=rep_files,
                source_type="epi"
            )

            if merged_counts.shape[0] == 0:
                print("[WARNING] empty merged count:", marker, tissue)
                continue

            out_count = os.path.join(
                marker_count_dir,
                "{}_{}_AS_merged_counts.tsv".format(
                    safe_name(marker),
                    tissue
                )
            )

            merged_counts.to_csv(out_count, sep="\t", index=False, na_rep="NA")
            marker_count_list.append(merged_counts)

        if len(marker_count_list) == 0:
            print("[WARNING] no count rows for marker:", marker)
            continue

        marker_counts = pd.concat(marker_count_list, axis=0, ignore_index=True)

        marker_counts_file = os.path.join(
            COUNT_DIR,
            "{}.all_tissues.AS_merged_counts.tsv".format(safe_name(marker))
        )

        marker_counts.to_csv(
            marker_counts_file,
            sep="\t",
            index=False,
            na_rep="NA"
        )

        marker_bed, marker_bed_df = write_as_bed(
            marker=marker,
            marker_df=marker_counts,
            source_type="epi"
        )

        overlap, overlap_file, raw_overlap_file = run_epi_as_enhancer_overlap(
            marker=marker,
            marker_bed=marker_bed,
            enhancer_bed=enhancer_bed
        )

        values_df = build_epi_analysis_values(overlap, marker)

        if values_df.shape[0] == 0:
            print("[WARNING] {} has no analysis values after enhancer overlap.".format(marker))

            match_rows.append({
                "marker": marker,
                "region_type": "enhancer",
                "n_tissue_pairs": len(file_map[marker]),
                "merged_count_rows": marker_counts.shape[0],
                "marker_bed_rows": marker_bed_df.shape[0],
                "raw_overlap_rows": overlap.shape[0],
                "analysis_values_rows_after_dedup": 0,
                "analysis_level": ANALYSIS_LEVEL,
                "effect_col": EFFECT_COL,
                "count_file": marker_counts_file,
                "overlap_file": overlap_file,
                "raw_overlap_file": raw_overlap_file,
            })

            continue

        summary, pairwise = summarize_and_test(values_df, marker)

        all_summaries.append(summary)
        all_pairwise.append(pairwise)
        all_values.append(values_df)

        match_rows.append({
            "marker": marker,
            "region_type": "enhancer",
            "n_tissue_pairs": len(file_map[marker]),
            "merged_count_rows": marker_counts.shape[0],
            "merged_count_unique_marker_tissue_variant": marker_counts[
                ["marker", "tissue", "variant_id"]
            ].drop_duplicates().shape[0],
            "marker_bed_rows": marker_bed_df.shape[0],
            "raw_overlap_rows": overlap.shape[0],
            "analysis_values_rows_after_dedup": values_df.shape[0],
            "analysis_level": ANALYSIS_LEVEL,
            "effect_col": EFFECT_COL,
            "count_file": marker_counts_file,
            "overlap_file": overlap_file,
            "raw_overlap_file": raw_overlap_file,
        })

    return all_summaries, all_pairwise, all_values, match_rows


def process_rnaseq_ase(target_gene_bed):
    """
      RNASeq ASE:
    ASE   overlap   target gene regions.
    """
    marker = "RNASeq_ASE"

    print("\n[INFO] Discovering RNASeq ASE files...")

    file_map, discovered_file = discover_as_files(
        as_dir=RNASEQ_ASE_DIR,
        parser_func=parse_rnaseq_ase_filename,
        markers_to_use=["RNASeq_ASE"],
        out_prefix="RNASeq_ASE"
    )

    print("[INFO] discovered RNASeq ASE file table:", discovered_file)

    if marker not in file_map:
        print("[WARNING] RNASeq_ASE files not found.")
        return [], [], [], [{
            "marker": marker,
            "region_type": "target_gene",
            "n_tissue_pairs": 0,
            "merged_count_rows": 0,
            "marker_bed_rows": 0,
            "raw_overlap_rows": 0,
            "analysis_values_rows_after_dedup": 0,
            "analysis_level": ANALYSIS_LEVEL,
            "effect_col": EFFECT_COL,
            "count_file": "NA",
            "overlap_file": "NA",
            "raw_overlap_file": "NA",
        }]

    count_list = []
    count_dir = os.path.join(COUNT_DIR, marker)
    os.makedirs(count_dir, exist_ok=True)

    for tissue in sorted(file_map[marker].keys()):
        rep_files = file_map[marker][tissue]

        merged_counts = merge_rep_counts(
            marker=marker,
            tissue=tissue,
            rep_files=rep_files,
            source_type="rnaseq"
        )

        if merged_counts.shape[0] == 0:
            print("[WARNING] empty RNASeq ASE merged count:", tissue)
            continue

        out_count = os.path.join(
            count_dir,
            "{}_{}_ASE_merged_counts.tsv".format(
                marker,
                tissue
            )
        )

        merged_counts.to_csv(out_count, sep="\t", index=False, na_rep="NA")
        count_list.append(merged_counts)

    if len(count_list) == 0:
        print("[WARNING] no RNASeq ASE count rows.")
        return [], [], [], [{
            "marker": marker,
            "region_type": "target_gene",
            "n_tissue_pairs": len(file_map[marker]),
            "merged_count_rows": 0,
            "marker_bed_rows": 0,
            "raw_overlap_rows": 0,
            "analysis_values_rows_after_dedup": 0,
            "analysis_level": ANALYSIS_LEVEL,
            "effect_col": EFFECT_COL,
            "count_file": "NA",
            "overlap_file": "NA",
            "raw_overlap_file": "NA",
        }]

    counts = pd.concat(count_list, axis=0, ignore_index=True)

    counts_file = os.path.join(
        COUNT_DIR,
        "{}.all_tissues.ASE_merged_counts.tsv".format(marker)
    )

    counts.to_csv(counts_file, sep="\t", index=False, na_rep="NA")

    marker_bed, marker_bed_df = write_as_bed(
        marker=marker,
        marker_df=counts,
        source_type="rnaseq"
    )

    overlap, overlap_file, raw_overlap_file = run_rnaseq_ase_target_gene_overlap(
        marker_bed=marker_bed,
        target_gene_bed=target_gene_bed
    )

    values_df = build_rnaseq_analysis_values(overlap)

    if values_df.shape[0] == 0:
        print("[WARNING] RNASeq_ASE has no analysis values after target-gene overlap.")

        return [], [], [], [{
            "marker": marker,
            "region_type": "target_gene",
            "n_tissue_pairs": len(file_map[marker]),
            "merged_count_rows": counts.shape[0],
            "marker_bed_rows": marker_bed_df.shape[0],
            "raw_overlap_rows": overlap.shape[0],
            "analysis_values_rows_after_dedup": 0,
            "analysis_level": ANALYSIS_LEVEL,
            "effect_col": EFFECT_COL,
            "count_file": counts_file,
            "overlap_file": overlap_file,
            "raw_overlap_file": raw_overlap_file,
        }]

    summary, pairwise = summarize_and_test(values_df, marker)

    match_rows = [{
        "marker": marker,
        "region_type": "target_gene",
        "n_tissue_pairs": len(file_map[marker]),
        "merged_count_rows": counts.shape[0],
        "merged_count_unique_marker_tissue_variant": counts[
            ["marker", "tissue", "variant_id"]
        ].drop_duplicates().shape[0],
        "marker_bed_rows": marker_bed_df.shape[0],
        "raw_overlap_rows": overlap.shape[0],
        "analysis_values_rows_after_dedup": values_df.shape[0],
        "analysis_level": ANALYSIS_LEVEL,
        "effect_col": EFFECT_COL,
        "count_file": counts_file,
        "overlap_file": overlap_file,
        "raw_overlap_file": raw_overlap_file,
    }]

    return [summary], [pairwise], [values_df], match_rows


def main():
    check_tools()

    print("[INFO] Reading enhancer groups...")
    enh = read_enhancer_groups(ENHANCER_FILE)

    enh_out = os.path.join(
        REGION_DIR,
        "all_5groups.enhancers.used_for_epigenomic_AS_overlap.tsv"
    )

    enh.to_csv(enh_out, sep="\t", index=False, na_rep="NA")

    enhancer_bed = write_enhancer_bed(enh)

    print("[INFO] Reading target-gene groups...")
    target_gene_regions, target_file, matched_file, unmatched_file = read_target_gene_groups(
        linked_pair_file=LINKED_PAIR_FILE,
        gene_bed=GENE_BED
    )

    target_gene_bed = write_target_gene_bed(target_gene_regions)

    target_gene_bed_copy = os.path.join(
        REGION_DIR,
        "all_5groups.target_gene_regions.used_for_RNASeq_ASE_overlap.bed"
    )

    target_gene_regions.to_csv(
        os.path.join(
            REGION_DIR,
            "all_5groups.target_gene_regions.used_for_RNASeq_ASE_overlap.tsv"
        ),
        sep="\t",
        index=False,
        na_rep="NA"
    )

    shutil.copyfile(target_gene_bed, target_gene_bed_copy)

    print("[INFO] target gene list:", target_file)
    print("[INFO] matched target gene regions:", matched_file)
    print("[INFO] unmatched target genes:", unmatched_file)
    print("[INFO] target gene BED:", target_gene_bed_copy)

    epi_summaries, epi_pairwise, epi_values, epi_match = process_epigenomic_as(
        enhancer_bed=enhancer_bed
    )

    rnaseq_summaries, rnaseq_pairwise, rnaseq_values, rnaseq_match = process_rnaseq_ase(
        target_gene_bed=target_gene_bed
    )

    all_summaries = epi_summaries + rnaseq_summaries
    all_pairwise = epi_pairwise + rnaseq_pairwise
    all_values = epi_values + rnaseq_values
    all_match_rows = epi_match + rnaseq_match

    if len(all_summaries) > 0:
        combined_summary = pd.concat(all_summaries, axis=0, ignore_index=True)

        combined_summary_file = os.path.join(
            SUMMARY_DIR,
            "six_panels.all_markers.{}.{}.group_summary_with_letters.tsv".format(
                ANALYSIS_LEVEL,
                EFFECT_COL
            )
        )

        combined_summary.to_csv(
            combined_summary_file,
            sep="\t",
            index=False,
            na_rep="NA"
        )
    else:
        combined_summary = pd.DataFrame()
        combined_summary_file = "NA"

    if len(all_pairwise) > 0:
        combined_pairwise = pd.concat(all_pairwise, axis=0, ignore_index=True)

        combined_pairwise_file = os.path.join(
            SUMMARY_DIR,
            "six_panels.all_markers.{}.{}.pairwise_mannwhitney_BH.tsv".format(
                ANALYSIS_LEVEL,
                EFFECT_COL
            )
        )

        combined_pairwise.to_csv(
            combined_pairwise_file,
            sep="\t",
            index=False,
            na_rep="NA"
        )
    else:
        combined_pairwise_file = "NA"

    if len(all_values) > 0:
        combined_values = pd.concat(all_values, axis=0, ignore_index=True)

        combined_values_file = os.path.join(
            VALUE_DIR,
            "six_panels.all_markers.{}.{}.values_for_violin.tsv".format(
                ANALYSIS_LEVEL,
                EFFECT_COL
            )
        )

        combined_values.to_csv(
            combined_values_file,
            sep="\t",
            index=False,
            na_rep="NA"
        )
    else:
        combined_values = pd.DataFrame()
        combined_values_file = "NA"

    match_summary = pd.DataFrame(all_match_rows)

    match_summary_file = os.path.join(
        OUT_DIR,
        "six_panels.AS_ASE_overlap_match_summary.tsv"
    )

    match_summary.to_csv(
        match_summary_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    if combined_values.shape[0] > 0 and combined_summary.shape[0] > 0:
        print("[INFO] Plotting six-panel figure...")
        plot_six_panel_figure(
            all_values=combined_values,
            all_summary=combined_summary
        )
    else:
        print("[WARNING] No combined values/summary available, figure skipped.")

    print("\n .")
    print("Output directory:", OUT_DIR)
    print("enhancer used:", enh_out)
    print("target genes used:", matched_file)
    print("unmatched target genes:", unmatched_file)
    print("match summary:", match_summary_file)
    print("combined values:", combined_values_file)
    print("combined summary:", combined_summary_file)
    print("combined pairwise:", combined_pairwise_file)
    print("six-panel figure dir:", FIG_DIR)


if __name__ == "__main__":
    main()
