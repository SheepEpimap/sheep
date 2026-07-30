#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
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
import matplotlib.patheffects as pe

from matplotlib.colors import TwoSlopeNorm
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

from scipy.stats import fisher_exact


# ============================================================
# ============================================================

CHROMBPNET_ANNOT_DIR = "/data/home/sczd644/run/zsw_chrombpnet/snpscore/02_annotation"

ENHANCER_FILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_enhancers.tsv"
)

OUTDIR = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "zutu/ATAC_ChromBPNet_motif_5groups_sig_only_BC_combined"
)

TMPDIR = os.path.join(OUTDIR, "tmp_bed")
PAIR_DIR = os.path.join(OUTDIR, "overlap_pairs")
FIGDIR = os.path.join(OUTDIR, "figures")
TABLEDIR = os.path.join(OUTDIR, "tables")

for d in [OUTDIR, TMPDIR, PAIR_DIR, FIGDIR, TABLEDIR]:
    os.makedirs(d, exist_ok=True)


# ============================================================
# ============================================================

LOGFC_PVAL_CUTOFF = 0.05

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

GROUP_COLOR = "#0072B2"

LEFT_YLABEL_X = -0.13


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


def parse_bool(x):
    """
      True/False,1/0,yes/no.
    """
    if pd.isna(x):
        return False

    s = str(x).strip()

    if s in ["True", "true", "TRUE", "1", "1.0", "Yes", "yes", "Y", "y"]:
        return True

    return False


def parse_motif_list(x):
    """
      hits_motifs.
    hits_motifs = "-"   motif,  list.
    """
    if pd.isna(x):
        return []

    s = str(x).strip()

    empty_values = {
        "",
        ".",
        "-",
        "NA",
        "NaN",
        "nan",
        "None",
        "none",
        "False",
        "false",
        "[]",
    }

    if s in empty_values:
        return []

    s = s.replace("[", "")
    s = s.replace("]", "")
    s = s.replace("'", "")
    s = s.replace('"', "")

    for sep in [";", "|"]:
        s = s.replace(sep, ",")

    parts = [p.strip() for p in s.split(",")]

    motifs = []

    for p in parts:
        if p in empty_values:
            continue

        motifs.append(p)

    motifs = sorted(set(motifs))

    return motifs


def infer_has_motif(row):
    """
      motif overlap.
     :
    1. hits_motifs  ;
    2.   hits_overlap  ,  hits_overlap   True.
    """
    motifs = row.get("motif_list", [])

    if not isinstance(motifs, list):
        motifs = parse_motif_list(motifs)

    if len(motifs) == 0:
        return False

    if "hits_overlap" in row.index:
        return parse_bool(row["hits_overlap"])

    return True


# ============================================================
# 4. read ChromBPNet annotation
# ============================================================

def parse_annotation_filename(filepath):
    """
    file :
    abomasum.annotations.tsv
    cerebral-cortex.annotations.tsv
    """
    name = os.path.basename(filepath)

    if not name.endswith(".annotations.tsv"):
        raise ValueError("  annotation file : {}".format(name))

    tissue = name.replace(".annotations.tsv", "")

    return tissue


def read_one_chrombpnet_annotation(filepath):
    """
    read  ChromBPNet annotation file.
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

    if "hits_motifs" in df.columns:
        df["motif_list"] = df["hits_motifs"].apply(parse_motif_list)
    else:
        df["hits_motifs"] = "-"
        df["motif_list"] = [[] for _ in range(df.shape[0])]

    if "hits_overlap" not in df.columns:
        df["hits_overlap"] = False

    df["has_motif_overlap"] = df.apply(infer_has_motif, axis=1)

    df["motif_list_str"] = df["motif_list"].apply(
        lambda x: ",".join(x) if isinstance(x, list) and len(x) > 0 else "NA"
    )

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
    files = sorted(glob.glob(os.path.join(annotation_dir, "*.annotations.tsv")))

    if len(files) == 0:
        raise FileNotFoundError("  ChromBPNet annotation file: {}".format(annotation_dir))

    all_list = []
    count_rows = []

    for f in files:
        tissue = parse_annotation_filename(f)
        df = read_one_chrombpnet_annotation(f)

        count_rows.append({
            "file": f,
            "tissue": tissue,
            "rows_after_parsing": df.shape[0],
            "rows_logfc_pval_lt_0p05": int((df["logfc.pval"] < LOGFC_PVAL_CUTOFF).sum()) if df.shape[0] > 0 else 0,
            "rows_with_true_motif": int(df["has_motif_overlap"].sum()) if df.shape[0] > 0 else 0,
            "rows_sig_with_true_motif": int(((df["logfc.pval"] < LOGFC_PVAL_CUTOFF) & (df["has_motif_overlap"] == True)).sum()) if df.shape[0] > 0 else 0,
        })

        if df.shape[0] > 0:
            all_list.append(df)
            print("[INFO] read annotation {}: {} rows".format(os.path.basename(f), df.shape[0]))

    count_file = os.path.join(TABLEDIR, "ChromBPNet_annotation_input_count.tsv")
    pd.DataFrame(count_rows).to_csv(count_file, sep="\t", index=False, na_rep="NA")

    if len(all_list) == 0:
        raise ValueError("  ChromBPNet annotation file .")

    anno = pd.concat(all_list, axis=0, ignore_index=True)
    anno = sort_bed_df(anno)

    return anno


# ============================================================
# 5. read enhancer grouping
# ============================================================

def read_enhancer_groups(enhancer_file):
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
        subset=["tissue_count_group", "enhancer_id", "chr", "start", "end"],
        keep="first"
    ).copy()

    enh["enhancer_length"] = enh["end"] - enh["start"]
    enh["group_label"] = enh["tissue_count_group"].map(GROUP_LABEL)

    enh = sort_bed_df(enh)

    return enh


# ============================================================
# ============================================================

def write_chrombpnet_sig_bed(anno_sig):
    """
      logfc.pval < 0.05   ATAC-AS   BED.
    """
    bed = anno_sig[
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

    outfile = os.path.join(
        TMPDIR,
        "significant_ATAC_ChromBPNet_AS_sites.logfc_pval_lt_0p05.bed"
    )

    bed.to_csv(outfile, sep="\t", header=False, index=False)

    return outfile


def write_enhancer_bed(enh):
    """
      enhancer BED.
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

    outfile = os.path.join(TMPDIR, "all_5groups.enhancers.bed")
    bed.to_csv(outfile, sep="\t", header=False, index=False)

    return outfile


def run_bedtools_intersect(chrombpnet_sig_bed, enhancer_bed):
    """
      ATAC-AS   enhancer   overlap.
    """
    outfile = os.path.join(
        PAIR_DIR,
        "significant_ATAC_ChromBPNet_overlap_5groups_enhancer_pairs.raw.tsv"
    )

    cmd = [
        "bedtools",
        "intersect",
        "-a",
        chrombpnet_sig_bed,
        "-b",
        enhancer_bed,
        "-wa",
        "-wb",
    ]

    with open(outfile, "w") as out:
        subprocess.run(
            cmd,
            stdout=out,
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

    if (not os.path.exists(outfile)) or os.path.getsize(outfile) == 0:
        pair = pd.DataFrame(columns=names)
    else:
        pair = pd.read_csv(
            outfile,
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


def merge_overlap_with_full_annotation(overlap_pairs, anno_sig):
    """
    overlap results  chrombpnet_record_id   ChromBPNet annotation.
    """
    if overlap_pairs.shape[0] == 0:
        return overlap_pairs.copy()

    merged = overlap_pairs.merge(
        anno_sig,
        on="chrombpnet_record_id",
        how="left",
        suffixes=("", "_chrombpnet")
    )

    merged["has_chrombpnet_annotation"] = merged["logfc.pval"].notna()
    merged["logfc.pval"] = pd.to_numeric(merged["logfc.pval"], errors="coerce")
    merged["is_logfc_pval_lt_0p05"] = merged["logfc.pval"] < LOGFC_PVAL_CUTOFF

    if "has_motif_overlap" not in merged.columns:
        merged["has_motif_overlap"] = False

    if "motif_list_str" not in merged.columns:
        merged["motif_list_str"] = "NA"

    return merged


# ============================================================
# ============================================================

def build_sig_event_level(merged):
    """
      tissue + variant_id + enhancer group  .
    """
    event = merged[
        (merged["has_chrombpnet_annotation"] == True) &
        (merged["is_logfc_pval_lt_0p05"] == True)
    ].drop_duplicates(
        subset=["tissue", "variant_id", "tissue_count_group"],
        keep="first"
    ).copy()

    event["has_motif_overlap"] = event["has_motif_overlap"].fillna(False).astype(bool)

    return event


def summarize_motif_overlap_by_group_sig_only(merged, enh):
    """
      ATAC-AS  statistics motif overlap.

     :
          overlap   enhancer   logfc.pval < 0.05   ATAC-AS  .

     :
          ATAC-AS   motif overlap  .
    """
    event = build_sig_event_level(merged)

    rows = []

    for group in GROUP_ORDER:
        sub_sig = event[event["tissue_count_group"] == group].copy()
        group_enh = enh[enh["tissue_count_group"] == group].copy()

        enhancer_total = group_enh["enhancer_id"].nunique()

        enhancer_total_bp = (
            group_enh[["enhancer_id", "enhancer_length"]]
            .drop_duplicates()
            ["enhancer_length"]
            .sum()
        )

        n_sig = sub_sig["chrombpnet_record_id"].nunique()

        n_sig_with_motif = sub_sig[
            sub_sig["has_motif_overlap"] == True
        ]["chrombpnet_record_id"].nunique()

        motif_percent_sig = n_sig_with_motif / n_sig * 100 if n_sig > 0 else np.nan

        sig_density = (
            n_sig / enhancer_total_bp * 1e6
            if enhancer_total_bp > 0 else np.nan
        )

        sig_motif_density = (
            n_sig_with_motif / enhancer_total_bp * 1e6
            if enhancer_total_bp > 0 else np.nan
        )

        rows.append({
            "tissue_count_group": group,
            "tissue_count_range": GROUP_LABEL[group],
            "enhancer_total": enhancer_total,
            "enhancer_total_bp": enhancer_total_bp,
            "significant_events_logfc_pval_lt_0p05": n_sig,
            "significant_events_with_motif": n_sig_with_motif,
            "significant_events_without_motif": n_sig - n_sig_with_motif,
            "significant_events_with_motif_percent": motif_percent_sig,
            "significant_events_density_per_Mb": sig_density,
            "significant_motif_events_density_per_Mb": sig_motif_density,
        })

    return pd.DataFrame(rows)


def fisher_motif_overlap_enrichment_by_group_sig_only(motif_summary):
    """
      ATAC-AS   Fisher  .

      group:
        motif-overlap significant AS vs non-motif significant AS

      groups:
        motif-overlap significant AS vs non-motif significant AS

     :
        alternative="greater"
    """
    rows = []

    for group in GROUP_ORDER:
        this = motif_summary[motif_summary["tissue_count_group"] == group].copy()
        rest = motif_summary[motif_summary["tissue_count_group"] != group].copy()

        if this.shape[0] == 0 or rest.shape[0] == 0:
            rows.append({
                "tissue_count_group": group,
                "fisher_or": np.nan,
                "log2_fisher_or": np.nan,
                "fisher_p_enrichment": np.nan,
            })
            continue

        a = int(this["significant_events_with_motif"].sum())
        b = int(this["significant_events_without_motif"].sum())
        c = int(rest["significant_events_with_motif"].sum())
        d = int(rest["significant_events_without_motif"].sum())

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
            "a_this_group_motif_overlap": a,
            "b_this_group_non_motif": b,
            "c_other_groups_motif_overlap": c,
            "d_other_groups_non_motif": d,
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

def draw_vector_log2or_heatmap(ax_heat, mat, qmat, group_labels, x, norm, cmap, vmax):
    """
      Rectangle   heatmap.

     :
    1.   imshow output  image/raster  ;
    2.   PDF   Illustrator   heatmap  ;
    3.   heatmap   Illustrator  .
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
        ["Motif-overlap"],
        fontsize=9
    )

    ax_heat.tick_params(
        axis="y",
        left=False,
        right=False,
        labelleft=True,
        pad=7
    )

    ax_heat.set_ylabel("")

    for spine in ax_heat.spines.values():
        spine.set_visible(False)


def draw_vector_horizontal_colorbar(ax_cbar, cmap, norm, vmin, vmax, n_steps=256):
    """
      Rectangle   colorbar.

    n_steps  ,colorbar  ;
    n_steps  ,Illustrator  .
      256, .
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


# ============================================================
# ============================================================

def plot_combined_BC_motif(motif_summary):
    """
      b/c  :

     :
        motif-overlap significant ATAC-AS sites per Mb enhancer

     :
        motif-overlap enrichment heatmap

     :
        1.   heatmap   x  ,x = 1,2,3,4,5;
        2. heatmap   0.5-5.5;
        3.   x  ;
        4.   panel a/d;
        5.  ;
        6. grouping  heatmap  ;
        7. colorbar  ;
        8. heatmap  ,  log2OR  ;
        9. heatmap   colorbar   Rectangle  ,  imshow/fig.colorbar.
    """
    group_labels = [GROUP_LABEL[g] for g in GROUP_ORDER]
    x = np.arange(1, len(GROUP_ORDER) + 1)

    sub = motif_summary.set_index("tissue_count_group")

    fig = plt.figure(figsize=(5.8, 5.0))

    left_x0 = 0.16
    left_w = 0.78

    cbar_y0 = 0.10
    cbar_h = 0.035

    gap_heat_cbar = 0.145

    heat_y0 = cbar_y0 + cbar_h + gap_heat_cbar

    fig_w, fig_h = fig.get_size_inches()
    heat_h = (left_w * fig_w / len(GROUP_ORDER)) / fig_h

    gap_line_heat = 0.015

    line_y0 = heat_y0 + heat_h + gap_line_heat
    line_h = 0.90 - line_y0

    ax_line = fig.add_axes([left_x0, line_y0, left_w, line_h])
    ax_heat = fig.add_axes([left_x0, heat_y0, left_w, heat_h], sharex=ax_line)
    ax_cbar = fig.add_axes([left_x0, cbar_y0, left_w, cbar_h])

    # ========================================================
    # ========================================================

    y_density = [
        sub.loc[g, "significant_motif_events_density_per_Mb"]
        if g in sub.index else np.nan
        for g in GROUP_ORDER
    ]

    ax_line.plot(
        x,
        y_density,
        marker="o",
        linewidth=2.2,
        markersize=5,
        color=GROUP_COLOR,
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
        "Motif-overlap significant ATAC-AS sites per Mb",
        fontsize=10
    )
    ax_line.yaxis.set_label_coords(LEFT_YLABEL_X, 0.5)

    nature_style_ax(ax_line)

    # ========================================================
    # ========================================================

    mat = pd.DataFrame(
        index=["Motif-overlap"],
        columns=group_labels,
        dtype=float
    )

    qmat = pd.DataFrame(
        index=["Motif-overlap"],
        columns=group_labels,
        dtype=float
    )

    for _, row in motif_summary.iterrows():
        group_label = row["tissue_count_range"]
        mat.loc["Motif-overlap", group_label] = row["log2_fisher_or"]
        qmat.loc["Motif-overlap", group_label] = row["fisher_q_BH_enrichment"]

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

    pdf = os.path.join(
        FIGDIR,
        "combined_BC_motif_overlap_density_and_enrichment_EnhA.pdf"
    )

    png = os.path.join(
        FIGDIR,
        "combined_BC_motif_overlap_density_and_enrichment_EnhA.png"
    )

    svg = os.path.join(
        FIGDIR,
        "combined_BC_motif_overlap_density_and_enrichment_EnhA.svg"
    )

    plt.savefig(
        pdf,
        bbox_inches="tight",
        dpi=600,
        facecolor="white",
        transparent=False
    )

    plt.savefig(
        png,
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        transparent=False
    )

    plt.savefig(
        svg,
        bbox_inches="tight",
        facecolor="white",
        transparent=False
    )

    plt.close()

    print("[INFO] Combined B/C motif figure:")
    print(pdf)
    print(png)
    print(svg)


# ============================================================
# 10. main workflow
# ============================================================

def main():
    check_bedtools()

    print("[INFO] Reading ChromBPNet ATAC annotation files...")
    anno_all = read_all_chrombpnet_annotations(CHROMBPNET_ANNOT_DIR)

    anno_all_file = os.path.join(
        TABLEDIR,
        "all_tissue.ATAC_ChromBPNet_annotations.all_sites.with_motif.tsv"
    )

    anno_all.to_csv(
        anno_all_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    # ========================================================
    # ========================================================

    anno_sig = anno_all[
        anno_all["logfc.pval"] < LOGFC_PVAL_CUTOFF
    ].copy()

    anno_sig = sort_bed_df(anno_sig)

    anno_sig_file = os.path.join(
        TABLEDIR,
        "all_tissue.ATAC_ChromBPNet_annotations.logfc_pval_lt_0p05.with_motif.tsv"
    )

    anno_sig.to_csv(
        anno_sig_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    print("[INFO] all annotation rows:", anno_all.shape[0])
    print("[INFO] significant annotation rows logfc.pval < 0.05:", anno_sig.shape[0])
    print("[INFO] significant rows with true motif:", int(anno_sig["has_motif_overlap"].sum()))
    print("[INFO] tissues in significant rows:", anno_sig["tissue"].nunique())

    if anno_sig.shape[0] == 0:
        raise ValueError("  logfc.pval < 0.05   ATAC-AS  , .")

    print("[INFO] Reading enhancer groups...")
    enh = read_enhancer_groups(ENHANCER_FILE)

    enh_file = os.path.join(
        TABLEDIR,
        "all_5groups.enhancers.used_for_sig_ATAC_ChromBPNet_motif_overlap.tsv"
    )

    enh.to_csv(
        enh_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    print("[INFO] enhancer rows:", enh.shape[0])

    print("[INFO] Writing BED files...")
    chrombpnet_sig_bed = write_chrombpnet_sig_bed(anno_sig)
    enhancer_bed = write_enhancer_bed(enh)

    print("[INFO] Running bedtools intersect with significant ATAC-AS sites only...")
    overlap_pairs = run_bedtools_intersect(
        chrombpnet_sig_bed,
        enhancer_bed
    )

    overlap_pair_file = os.path.join(
        TABLEDIR,
        "significant_ATAC_ChromBPNet_overlap_5groups.enhancer_pairs.before_full_annotation.tsv"
    )

    overlap_pairs.to_csv(
        overlap_pair_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    print("[INFO] significant overlap pairs:", overlap_pairs.shape[0])

    print("[INFO] Merging overlap pairs with full significant ChromBPNet annotation...")
    merged = merge_overlap_with_full_annotation(
        overlap_pairs,
        anno_sig
    )

    all_overlap_full_file = os.path.join(
        TABLEDIR,
        "significant_ATAC_ChromBPNet_overlap_5groups.enhancer_pairs.with_full_annotation_and_motif.tsv"
    )

    merged.to_csv(
        all_overlap_full_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    # ========================================================
    # ========================================================

    sig = merged[
        (merged["has_chrombpnet_annotation"] == True) &
        (merged["is_logfc_pval_lt_0p05"] == True)
    ].copy()

    front_cols = [
        "tissue",
        "variant_id",
        "chrombpnet_record_id",
        "as_chr",
        "as_start",
        "as_end",
        "as_pos",
        "allele1",
        "allele2",
        "enhancer_id",
        "enh_chr",
        "enh_start",
        "enh_end",
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_length",
        "logfc",
        "abs_logfc",
        "logfc.pval",
        "abs_logfc.pval",
        "jsd",
        "jsd.pval",
        "has_motif_overlap",
        "hits_overlap",
        "hits_motifs",
        "motif_list_str",
    ]

    front_cols = [c for c in front_cols if c in sig.columns]
    rest_cols = [c for c in sig.columns if c not in front_cols]

    sig = sig[front_cols + rest_cols].copy()

    sig_file = os.path.join(
        TABLEDIR,
        "significant_logfc_pval_lt_0p05.ATAC_overlap_enhancer.ChromBPNet_full.with_motif.tsv"
    )

    sig.to_csv(
        sig_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    sig_unique = sig.drop_duplicates(
        subset=["tissue", "variant_id", "tissue_count_group"],
        keep="first"
    ).copy()

    sig_unique_file = os.path.join(
        TABLEDIR,
        "significant_logfc_pval_lt_0p05.ATAC_overlap_enhancer.ChromBPNet_unique_events.with_motif.tsv"
    )

    sig_unique.to_csv(
        sig_unique_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    # ========================================================
    # motif summary
    # ========================================================

    print("[INFO] Summarizing motif overlap among significant ATAC-AS sites only...")

    motif_summary = summarize_motif_overlap_by_group_sig_only(
        merged,
        enh
    )

    motif_fisher = fisher_motif_overlap_enrichment_by_group_sig_only(
        motif_summary
    )

    motif_summary = motif_summary.merge(
        motif_fisher,
        on="tissue_count_group",
        how="left"
    )

    motif_summary_file = os.path.join(
        TABLEDIR,
        "ATAC_ChromBPNet_motif_overlap_5groups.significant_only.summary.tsv"
    )

    motif_summary.to_csv(
        motif_summary_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    # ========================================================
    # match summary
    # ========================================================

    event_level = build_sig_event_level(merged)

    match_summary = pd.DataFrame({
        "item": [
            "ChromBPNet_ATAC_annotation_rows_all",
            "ChromBPNet_ATAC_annotation_rows_logfc_pval_lt_0p05",
            "ChromBPNet_ATAC_sig_rows_with_true_motif",
            "enhancer_rows_input",
            "significant_overlap_pairs_before_full_annotation",
            "significant_overlap_pairs_with_full_annotation",
            "significant_overlap_unique_events",
            "significant_overlap_unique_events_with_motif",
            "significant_unique_events_output",
        ],
        "n": [
            anno_all.shape[0],
            anno_sig.shape[0],
            int(anno_sig["has_motif_overlap"].sum()),
            enh.shape[0],
            overlap_pairs.shape[0],
            int(merged["has_chrombpnet_annotation"].sum()),
            event_level["chrombpnet_record_id"].nunique(),
            event_level[event_level["has_motif_overlap"] == True]["chrombpnet_record_id"].nunique(),
            sig_unique.shape[0],
        ]
    })

    match_file = os.path.join(
        TABLEDIR,
        "ATAC_ChromBPNet_motif_overlap_match_summary.significant_only.tsv"
    )

    match_summary.to_csv(
        match_file,
        sep="\t",
        index=False
    )

    # ========================================================
    # ========================================================

    print("[INFO] Plotting combined B/C motif figure only...")
    plot_combined_BC_motif(motif_summary)

    print("")
    print(" . output b/c  , output a/d.")
    print("Output directory:", OUTDIR)

    print("")
    print(" Summary table:")
    print(motif_summary_file)
    print(match_file)
    print(sig_file)
    print(sig_unique_file)
    print(all_overlap_full_file)
    print(overlap_pair_file)
    print(anno_sig_file)
    print(enh_file)

    print("")
    print(" :")
    print(os.path.join(
        FIGDIR,
        "combined_BC_motif_overlap_density_and_enrichment_EnhA.pdf"
    ))
    print(os.path.join(
        FIGDIR,
        "combined_BC_motif_overlap_density_and_enrichment_EnhA.png"
    ))
    print(os.path.join(
        FIGDIR,
        "combined_BC_motif_overlap_density_and_enrichment_EnhA.svg"
    ))


if __name__ == "__main__":
    main()
