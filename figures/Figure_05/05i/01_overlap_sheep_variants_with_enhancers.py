#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import re
import gzip
import bisect
import math
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
from scipy.stats import kruskal, rankdata, norm


# ============================================================
# 1. inputoutput
# ============================================================

VCF = "/vol2/zhangshiwen/GWAS/locus/chrAuto/chrAuto.vcf.gz"

ENHANCER_FILE = (
    "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/"
    "E5_5groups_simple_target_gene/all_5groups.simple_enhancers.tsv"
)

OUTDIR = "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/enhancer_5group_snp_overlap_python"

FIGDIR = os.path.join(OUTDIR, "figures")
STATDIR = os.path.join(OUTDIR, "statistics")
DETAILDIR = os.path.join(OUTDIR, "overlap_details")

for d in [OUTDIR, FIGDIR, STATDIR, DETAILDIR]:
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
    "G1_1_tissue": "1",
    "G2_2_5_tissues": "2-5",
    "G3_6_10_tissues": "6-10",
    "G4_11_20_tissues": "11-20",
    "G5_21_43_tissues": "21-43",
}

GROUP_COLORS = {
    "G1_1_tissue": "#4C78A8",
    "G2_2_5_tissues": "#F58518",
    "G3_6_10_tissues": "#54A24B",
    "G4_11_20_tissues": "#B279A2",
    "G5_21_43_tissues": "#E45756",
}

ALPHA = 0.05

VARIANT_LEVEL = "site"

WRITE_ALL_SNP_BED = True

PLOT_LOG10_BURDEN = True


# ============================================================
# 3. helper functions
# ============================================================

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


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan


def bh_adjust(pvalues):
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


def format_pvalue(p):
    if pd.isna(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return "{:.1e}".format(p)
    return "{:.3f}".format(p)


def nature_style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10, length=4, width=1.0)
    ax.grid(False)


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


# ============================================================
# ============================================================

def read_enhancers(path):
    if not os.path.exists(path):
        raise FileNotFoundError("Not found enhancer file: {}".format(path))

    enh = pd.read_csv(path, sep="\t", dtype=str)

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

    enh = enh.dropna(subset=["chr", "start", "end", "tissue_count_group"]).copy()
    enh["start"] = enh["start"].astype(int)
    enh["end"] = enh["end"].astype(int)

    enh = enh[enh["end"] > enh["start"]].copy()
    enh = enh[enh["tissue_count_group"].isin(GROUP_ORDER)].copy()

    enh = enh.drop_duplicates(
        subset=["tissue_count_group", "chr", "start", "end"],
        keep="first"
    ).copy()

    enh["enhancer_length"] = enh["end"] - enh["start"]
    enh["group_label"] = enh["tissue_count_group"].map(GROUP_LABEL)

    enh["_chr_order"] = enh["chr"].map(chr_sort_key)
    enh = enh.sort_values(["_chr_order", "start", "end"])
    enh = enh.drop(columns=["_chr_order"])

    enh["enhancer_key"] = (
        enh["tissue_count_group"].astype(str) + "|" +
        enh["chr"].astype(str) + ":" +
        enh["start"].astype(str) + "-" +
        enh["end"].astype(str)
    )

    out = os.path.join(OUTDIR, "01_enhancers.used.no_merge.tsv")
    enh.to_csv(out, sep="\t", index=False, na_rep="NA")

    return enh


def build_interval_index(enh):
    """
      per chromosome interval index.
      enhancer   merge  ,  start  .
     , /  overlap  .
    """
    index = {}

    for chrom, sub in enh.groupby("chr", sort=False):
        rows = []

        for _, r in sub.iterrows():
            rows.append({
                "chr": r["chr"],
                "start": int(r["start"]),
                "end": int(r["end"]),
                "enhancer_id": r["enhancer_id"],
                "enhancer_key": r["enhancer_key"],
                "tissue_count_group": r["tissue_count_group"],
                "tissue_count_range": r["tissue_count_range"],
                "group_label": r["group_label"],
                "n_tissues": int(r["n_tissues"]) if pd.notna(r["n_tissues"]) else np.nan,
                "enhancer_length": int(r["enhancer_length"]),
            })

        rows = sorted(rows, key=lambda x: (x["start"], x["end"]))
        starts = [r["start"] for r in rows]

        index[chrom] = {
            "starts": starts,
            "rows": rows,
        }

    return index


def find_overlapping_enhancers(interval_index, chrom, pos0):
    """
    SNP BED   [pos0, pos0+1)
    enhancer overlap  :
    enh.start < pos0+1 and enh.end > pos0
      SNP,  enh.start <= pos0 < enh.end.
    """
    chrom = normalize_chr(chrom)

    if chrom not in interval_index:
        return []

    starts = interval_index[chrom]["starts"]
    rows = interval_index[chrom]["rows"]

    idx = bisect.bisect_right(starts, pos0)

    hits = []

    j = idx - 1
    while j >= 0:
        r = rows[j]
        if r["end"] <= pos0:
            break
        if r["start"] <= pos0 < r["end"]:
            hits.append(r)
        j -= 1

    j = idx
    while j < len(rows):
        r = rows[j]
        if r["start"] > pos0:
            break
        if r["start"] <= pos0 < r["end"]:
            hits.append(r)
        j += 1

    return hits


# ============================================================
# ============================================================

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def make_site_variant_id(chrom, pos):
    return "{}_{}".format(normalize_chr(chrom), pos)


def make_allele_variant_id(chrom, pos, ref, alt):
    return "{}_{}_{}_{}".format(normalize_chr(chrom), pos, ref, alt)


def parse_vcf_snps(vcf_path):
    """
     read VCF,  SNP site   SNP allele.

      generator, :
    {
      chrom, pos, start, end, ref, alt, vcf_id, variant_id
    }
    """
    with open_vcf(vcf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom = normalize_chr(fields[0])
            pos = int(fields[1])
            vcf_id = fields[2]
            ref = fields[3]
            alt_field = fields[4]

            if len(ref) != 1:
                continue

            alts = [a for a in alt_field.split(",") if a not in [".", "*", "<*>"]]
            snp_alts = [a for a in alts if len(a) == 1]

            if len(snp_alts) == 0:
                continue

            start = pos - 1
            end = pos

            if VARIANT_LEVEL == "site":
                yield {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "pos": pos,
                    "ref": ref,
                    "alt": ",".join(snp_alts),
                    "vcf_id": vcf_id,
                    "variant_id": make_site_variant_id(chrom, pos),
                }

            elif VARIANT_LEVEL == "allele":
                for alt in snp_alts:
                    yield {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "vcf_id": vcf_id,
                        "variant_id": make_allele_variant_id(chrom, pos, ref, alt),
                    }
            else:
                raise ValueError("VARIANT_LEVEL   site   allele.")


def scan_vcf_and_overlap(vcf_path, interval_index, enh):
    all_snp_bed_file = os.path.join(DETAILDIR, "02_vcf_snps.used.bed")
    unique_overlap_file = os.path.join(DETAILDIR, "03_unique_snp_overlap_by_group.tsv")
    detail_overlap_file = os.path.join(DETAILDIR, "04_snp_enhancer_overlap_details.tsv")

    all_snp_out = open(all_snp_bed_file, "w") if WRITE_ALL_SNP_BED else None

    unique_out = open(unique_overlap_file, "w")
    detail_out = open(detail_overlap_file, "w")

    if all_snp_out is not None:
        all_snp_out.write(
            "chr\tstart\tend\tvariant_id\tvcf_id\tref\talt\n"
        )

    unique_out.write(
        "group\tgroup_label\tvariant_chr\tvariant_start\tvariant_end\tvariant_id\tvcf_id\tref\talt\n"
    )

    detail_out.write(
        "group\tgroup_label\tvariant_chr\tvariant_start\tvariant_end\tvariant_id\tvcf_id\tref\talt\t"
        "enh_chr\tenh_start\tenh_end\tenhancer_id\tenhancer_key\tn_tissues\tenhancer_length\n"
    )

    total_snp_n = 0

    group_snp_sets = {g: set() for g in GROUP_ORDER}

    enhancer_snp_sets = {k: set() for k in enh["enhancer_key"].tolist()}

    for snp in parse_vcf_snps(vcf_path):
        total_snp_n += 1

        if all_snp_out is not None:
            all_snp_out.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    snp["chrom"],
                    snp["start"],
                    snp["end"],
                    snp["variant_id"],
                    snp["vcf_id"],
                    snp["ref"],
                    snp["alt"],
                )
            )

        hits = find_overlapping_enhancers(
            interval_index=interval_index,
            chrom=snp["chrom"],
            pos0=snp["start"]
        )

        if len(hits) == 0:
            continue

        groups_hit_this_snp = set()

        for h in hits:
            g = h["tissue_count_group"]
            group_snp_sets[g].add(snp["variant_id"])
            enhancer_snp_sets[h["enhancer_key"]].add(snp["variant_id"])

            if g not in groups_hit_this_snp:
                unique_out.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        g,
                        GROUP_LABEL[g],
                        snp["chrom"],
                        snp["start"],
                        snp["end"],
                        snp["variant_id"],
                        snp["vcf_id"],
                        snp["ref"],
                        snp["alt"],
                    )
                )
                groups_hit_this_snp.add(g)

            detail_out.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    g,
                    GROUP_LABEL[g],
                    snp["chrom"],
                    snp["start"],
                    snp["end"],
                    snp["variant_id"],
                    snp["vcf_id"],
                    snp["ref"],
                    snp["alt"],
                    h["chr"],
                    h["start"],
                    h["end"],
                    h["enhancer_id"],
                    h["enhancer_key"],
                    h["n_tissues"],
                    h["enhancer_length"],
                )
            )

    if all_snp_out is not None:
        all_snp_out.close()

    unique_out.close()
    detail_out.close()

    return {
        "total_snp_n": total_snp_n,
        "group_snp_sets": group_snp_sets,
        "enhancer_snp_sets": enhancer_snp_sets,
        "all_snp_bed_file": all_snp_bed_file if WRITE_ALL_SNP_BED else "NA",
        "unique_overlap_file": unique_overlap_file,
        "detail_overlap_file": detail_overlap_file,
    }


# ============================================================
# ============================================================

def build_summary_and_burden(enh, overlap_result):
    total_snp_n = overlap_result["total_snp_n"]
    group_snp_sets = overlap_result["group_snp_sets"]
    enhancer_snp_sets = overlap_result["enhancer_snp_sets"]

    burden_rows = []

    for _, r in enh.iterrows():
        k = r["enhancer_key"]
        snp_count = len(enhancer_snp_sets.get(k, set()))

        burden_rows.append({
            "group": r["tissue_count_group"],
            "group_label": r["group_label"],
            "chr": r["chr"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "enhancer_id": r["enhancer_id"],
            "enhancer_key": r["enhancer_key"],
            "n_tissues": r["n_tissues"],
            "enhancer_length": int(r["enhancer_length"]),
            "snp_count": snp_count,
            "snp_count_per_kb": snp_count / int(r["enhancer_length"]) * 1000 if int(r["enhancer_length"]) > 0 else np.nan,
            "log10_snp_count_plus1": np.log10(snp_count + 1),
            "has_snp": snp_count > 0,
        })

    burden = pd.DataFrame(burden_rows)

    summary_rows = []

    for g in GROUP_ORDER:
        sub_enh = enh[enh["tissue_count_group"] == g].copy()
        sub_burden = burden[burden["group"] == g].copy()

        enhancer_n = sub_enh.shape[0]
        total_bp = int(sub_enh["enhancer_length"].sum())
        total_mb = total_bp / 1_000_000 if total_bp > 0 else np.nan

        overlap_unique_snp_n = len(group_snp_sets[g])
        snp_density_per_mb = overlap_unique_snp_n / total_bp * 1_000_000 if total_bp > 0 else np.nan

        pct_total_snps = overlap_unique_snp_n / total_snp_n * 100 if total_snp_n > 0 else np.nan

        enhancers_with_snp_n = int((sub_burden["snp_count"] > 0).sum())
        enhancers_with_snp_pct = enhancers_with_snp_n / enhancer_n * 100 if enhancer_n > 0 else np.nan

        summary_rows.append({
            "group": g,
            "group_label": GROUP_LABEL[g],
            "enhancer_n": enhancer_n,
            "total_enhancer_bp_no_merge": total_bp,
            "total_enhancer_mb_no_merge": total_mb,
            "overlap_unique_snp_n": overlap_unique_snp_n,
            "overlap_unique_snp_per_mb": snp_density_per_mb,
            "pct_total_snps": pct_total_snps,
            "enhancers_with_snp_n": enhancers_with_snp_n,
            "enhancers_with_snp_pct": enhancers_with_snp_pct,
            "mean_snp_count_per_enhancer_all": sub_burden["snp_count"].mean() if sub_burden.shape[0] > 0 else np.nan,
            "median_snp_count_per_enhancer_all": sub_burden["snp_count"].median() if sub_burden.shape[0] > 0 else np.nan,
            "mean_snp_count_per_kb_all": sub_burden["snp_count_per_kb"].mean() if sub_burden.shape[0] > 0 else np.nan,
            "median_snp_count_per_kb_all": sub_burden["snp_count_per_kb"].median() if sub_burden.shape[0] > 0 else np.nan,
        })

    summary = pd.DataFrame(summary_rows)

    summary["group"] = pd.Categorical(summary["group"], categories=GROUP_ORDER, ordered=True)
    summary = summary.sort_values("group").copy()
    summary["group"] = summary["group"].astype(str)

    burden["group"] = pd.Categorical(burden["group"], categories=GROUP_ORDER, ordered=True)
    burden = burden.sort_values(["group", "chr", "start", "end"]).copy()
    burden["group"] = burden["group"].astype(str)

    summary_file = os.path.join(OUTDIR, "05_group_summary.no_enhancer_merge.tsv")
    burden_file = os.path.join(OUTDIR, "06_enhancer_snp_burden.all_enhancers.tsv")

    summary.to_csv(summary_file, sep="\t", index=False, na_rep="NA")
    burden.to_csv(burden_file, sep="\t", index=False, na_rep="NA")

    return summary, burden, summary_file, burden_file


# ============================================================
# 7. Kruskal + Dunn post hoc
# ============================================================

def dunn_posthoc_bh(values, groups, group_order):
    """
      Dunn test pairwise comparison.

    input:
    values:
    groups: grouping
    group_order:  grouping

    output:
    pairwise Dunn test   + adjusted P matrix
    """
    df = pd.DataFrame({
        "value": values,
        "group": groups,
    }).dropna().copy()

    df = df[df["group"].isin(group_order)].copy()

    if df.shape[0] == 0 or df["group"].nunique() < 2:
        empty = pd.DataFrame(columns=[
            "group1", "group2", "group1_label", "group2_label",
            "mean_rank1", "mean_rank2", "z", "p_raw", "p_adj_BH"
        ])
        pmat = pd.DataFrame(np.nan, index=group_order, columns=group_order)
        return empty, pmat

    df["rank"] = rankdata(df["value"].values, method="average")

    N = df.shape[0]

    # tie correction
    _, counts = np.unique(df["value"].values, return_counts=True)

    if N > 1:
        tie_sum = np.sum(counts**3 - counts)
        tie_correction = 1.0 - tie_sum / (N**3 - N)
    else:
        tie_correction = 1.0

    if tie_correction <= 0:
        tie_correction = np.nan

    group_stats = {}

    for g in group_order:
        sub = df[df["group"] == g]
        group_stats[g] = {
            "n": sub.shape[0],
            "mean_rank": sub["rank"].mean() if sub.shape[0] > 0 else np.nan,
        }

    rows = []

    for i in range(len(group_order)):
        for j in range(i + 1, len(group_order)):
            g1 = group_order[i]
            g2 = group_order[j]

            n1 = group_stats[g1]["n"]
            n2 = group_stats[g2]["n"]
            r1 = group_stats[g1]["mean_rank"]
            r2 = group_stats[g2]["mean_rank"]

            if n1 == 0 or n2 == 0 or pd.isna(tie_correction):
                z = np.nan
                p = np.nan
            else:
                se = math.sqrt(
                    (N * (N + 1) / 12.0) *
                    tie_correction *
                    (1.0 / n1 + 1.0 / n2)
                )

                if se == 0:
                    z = np.nan
                    p = np.nan
                else:
                    z = (r1 - r2) / se
                    p = 2.0 * norm.sf(abs(z))

            rows.append({
                "group1": g1,
                "group2": g2,
                "group1_label": GROUP_LABEL[g1],
                "group2_label": GROUP_LABEL[g2],
                "n1": n1,
                "n2": n2,
                "mean_rank1": r1,
                "mean_rank2": r2,
                "z": z,
                "p_raw": p,
            })

    pairwise = pd.DataFrame(rows)
    pairwise["p_adj_BH"] = bh_adjust(pairwise["p_raw"].values)

    pmat = pd.DataFrame(1.0, index=group_order, columns=group_order, dtype=float)

    for _, row in pairwise.iterrows():
        g1 = row["group1"]
        g2 = row["group2"]
        p = row["p_adj_BH"]

        if pd.notna(p):
            pmat.loc[g1, g2] = p
            pmat.loc[g2, g1] = p
        else:
            pmat.loc[g1, g2] = np.nan
            pmat.loc[g2, g1] = np.nan

    return pairwise, pmat


def compact_letters_from_pmat(pmat, means, alpha=0.05):
    """
      pairwise adjusted P  .
      mean  .
    """
    groups = list(pmat.index)

    order = sorted(
        groups,
        key=lambda g: means.get(g, -np.inf) if pd.notna(means.get(g, np.nan)) else -np.inf,
        reverse=True
    )

    letters = {g: "" for g in groups}
    letter_pool = list("abcdefghijklmnopqrstuvwxyz")

    for g in order:
        if pd.isna(means.get(g, np.nan)):
            letters[g] = "NA"
            continue

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


def run_kruskal_dunn(burden):
    burden = burden.copy()

    if PLOT_LOG10_BURDEN:
        value_col = "log10_snp_count_plus1"
        value_label = "log10(SNP count + 1)"
    else:
        value_col = "snp_count"
        value_label = "SNP count"

    arrays = []
    valid_groups = []

    means = {}

    for g in GROUP_ORDER:
        vals = burden.loc[burden["group"] == g, value_col].dropna().astype(float).values

        means[g] = np.mean(vals) if len(vals) > 0 else np.nan

        if len(vals) > 0:
            arrays.append(vals)
            valid_groups.append(g)

    if len(valid_groups) >= 2:
        try:
            kw = kruskal(*arrays)
            kruskal_H = kw.statistic
            kruskal_p = kw.pvalue
        except Exception:
            kruskal_H = np.nan
            kruskal_p = np.nan
    else:
        kruskal_H = np.nan
        kruskal_p = np.nan

    pairwise, pmat = dunn_posthoc_bh(
        values=burden[value_col].values,
        groups=burden["group"].values,
        group_order=GROUP_ORDER
    )

    letters = compact_letters_from_pmat(
        pmat=pmat,
        means=means,
        alpha=ALPHA
    )

    stat_rows = []

    for g in GROUP_ORDER:
        vals_raw = burden.loc[burden["group"] == g, "snp_count"].dropna().astype(float)
        vals_plot = burden.loc[burden["group"] == g, value_col].dropna().astype(float)

        stat_rows.append({
            "group": g,
            "group_label": GROUP_LABEL[g],
            "n_enhancers": int(vals_raw.shape[0]),
            "test_value_col": value_col,
            "test_value_label": value_label,
            "mean_snp_count": vals_raw.mean() if vals_raw.shape[0] > 0 else np.nan,
            "median_snp_count": vals_raw.median() if vals_raw.shape[0] > 0 else np.nan,
            "q25_snp_count": vals_raw.quantile(0.25) if vals_raw.shape[0] > 0 else np.nan,
            "q75_snp_count": vals_raw.quantile(0.75) if vals_raw.shape[0] > 0 else np.nan,
            "mean_test_value": vals_plot.mean() if vals_plot.shape[0] > 0 else np.nan,
            "median_test_value": vals_plot.median() if vals_plot.shape[0] > 0 else np.nan,
            "kruskal_H": kruskal_H,
            "kruskal_p": kruskal_p,
            "letter": letters.get(g, "NA"),
            "letter_source": "Dunn_test_BH_adjusted_P",
        })

    stat_summary = pd.DataFrame(stat_rows)

    stat_summary_file = os.path.join(STATDIR, "07_snp_burden.kruskal_dunn.summary.tsv")
    pairwise_file = os.path.join(STATDIR, "08_snp_burden.dunn_pairwise_BH.tsv")
    pmat_file = os.path.join(STATDIR, "09_snp_burden.dunn_pairwise_BH_pmatrix.tsv")

    stat_summary.to_csv(stat_summary_file, sep="\t", index=False, na_rep="NA")
    pairwise.to_csv(pairwise_file, sep="\t", index=False, na_rep="NA")
    pmat.to_csv(pmat_file, sep="\t", na_rep="NA")

    return stat_summary, pairwise, pmat, stat_summary_file, pairwise_file, pmat_file


# ============================================================
# ============================================================

def plot_big_figure(summary, burden, stat_summary):
    summary = summary.copy()
    burden = burden.copy()

    summary["group"] = pd.Categorical(summary["group"], categories=GROUP_ORDER, ordered=True)
    summary = summary.sort_values("group").copy()
    summary["group"] = summary["group"].astype(str)

    burden["group"] = pd.Categorical(burden["group"], categories=GROUP_ORDER, ordered=True)
    burden = burden.sort_values("group").copy()
    burden["group"] = burden["group"].astype(str)

    stat_summary["group"] = pd.Categorical(stat_summary["group"], categories=GROUP_ORDER, ordered=True)
    stat_summary = stat_summary.sort_values("group").copy()
    stat_summary["group"] = stat_summary["group"].astype(str)

    labels = [GROUP_LABEL[g] for g in GROUP_ORDER]
    colors = [GROUP_COLORS[g] for g in GROUP_ORDER]
    x = np.arange(len(GROUP_ORDER))

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5))

    ax1, ax2, ax3, ax4 = axes.flatten()

    # --------------------------------------------------------
    # A. Unique SNP count
    # --------------------------------------------------------
    y1 = summary["overlap_unique_snp_n"].values

    ax1.bar(x, y1, color=colors, edgecolor="black", linewidth=0.7)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_xlabel("Enhancer tissue-count group")
    ax1.set_ylabel("Unique overlapping SNP count")
    ax1.set_title("Unique SNPs overlapping enhancer groups", fontsize=11)
    nature_style_ax(ax1)

    for i, v in enumerate(y1):
        ax1.text(i, v, "{:,}".format(int(v)), ha="center", va="bottom", fontsize=8, rotation=0)

    ax1.text(-0.16, 1.08, "a", transform=ax1.transAxes, fontsize=14, fontweight="bold")

    # --------------------------------------------------------
    # B. SNP density per Mb
    # --------------------------------------------------------
    y2 = summary["overlap_unique_snp_per_mb"].values

    ax2.plot(x, y2, marker="o", color="black", linewidth=1.8, markersize=5)
    ax2.scatter(x, y2, color=colors, edgecolor="black", linewidth=0.7, s=55, zorder=3)

    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_xlabel("Enhancer tissue-count group")
    ax2.set_ylabel("Unique SNPs per Mb enhancer")
    ax2.set_title("SNP density normalized by enhancer length", fontsize=11)
    nature_style_ax(ax2)

    for i, v in enumerate(y2):
        ax2.text(i, v, "{:.2f}".format(v), ha="center", va="bottom", fontsize=8)

    ax2.text(-0.16, 1.08, "b", transform=ax2.transAxes, fontsize=14, fontweight="bold")

    # --------------------------------------------------------
    # C. Fraction of enhancers with SNP
    # --------------------------------------------------------
    y3 = summary["enhancers_with_snp_pct"].values

    ax3.bar(x, y3, color=colors, edgecolor="black", linewidth=0.7)
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels)
    ax3.set_xlabel("Enhancer tissue-count group")
    ax3.set_ylabel("Enhancers with ≥1 SNP (%)")
    ax3.set_title("Fraction of enhancers overlapped by SNPs", fontsize=11)
    nature_style_ax(ax3)

    for i, v in enumerate(y3):
        ax3.text(i, v, "{:.1f}%".format(v), ha="center", va="bottom", fontsize=8)

    ax3.text(-0.16, 1.08, "c", transform=ax3.transAxes, fontsize=14, fontweight="bold")

    # --------------------------------------------------------
    # D. SNP burden per enhancer with Kruskal + Dunn letters
    # --------------------------------------------------------
    if PLOT_LOG10_BURDEN:
        y_col = "log10_snp_count_plus1"
        y_label = "log10(SNP count per enhancer + 1)"
    else:
        y_col = "snp_count"
        y_label = "SNP count per enhancer"

    violin_data = []
    violin_pos = []
    box_data = []
    box_pos = []
    mean_x = []
    mean_y = []

    for i, g in enumerate(GROUP_ORDER, start=1):
        vals = burden.loc[burden["group"] == g, y_col].dropna().astype(float).values

        if len(vals) >= 2:
            violin_data.append(vals)
            violin_pos.append(i)

        if len(vals) > 0:
            box_data.append(vals)
            box_pos.append(i)
            mean_x.append(i)
            mean_y.append(np.mean(vals))

    if len(violin_data) > 0:
        parts = ax4.violinplot(
            violin_data,
            positions=violin_pos,
            widths=0.78,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

        pos_to_group = {i + 1: g for i, g in enumerate(GROUP_ORDER)}

        for body, pos in zip(parts["bodies"], violin_pos):
            g = pos_to_group[pos]
            body.set_facecolor(GROUP_COLORS[g])
            body.set_edgecolor("none")
            body.set_alpha(0.75)

    if len(box_data) > 0:
        ax4.boxplot(
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

    if len(mean_x) > 0:
        ax4.scatter(mean_x, mean_y, marker="D", color="black", s=18, zorder=4)

    ax4.set_xticks(np.arange(1, len(GROUP_ORDER) + 1))
    ax4.set_xticklabels(labels)
    ax4.set_xlabel("Enhancer tissue-count group")
    ax4.set_ylabel(y_label)
    ax4.set_title("SNP burden per enhancer", fontsize=11)
    nature_style_ax(ax4)

    all_vals = burden[y_col].dropna().astype(float)

    if all_vals.shape[0] > 0:
        ymin = all_vals.min()
        ymax = all_vals.max()
        yr = ymax - ymin if ymax > ymin else 1.0

        letter_y = ymax + 0.06 * yr
        ax4.set_ylim(ymin - 0.04 * yr, ymax + 0.23 * yr)

        letter_dict = dict(zip(stat_summary["group"], stat_summary["letter"]))

        for i, g in enumerate(GROUP_ORDER, start=1):
            letter = letter_dict.get(g, "NA")
            if letter != "NA":
                ax4.text(
                    i,
                    letter_y,
                    str(letter),
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    color="black"
                )

    if stat_summary.shape[0] > 0:
        kp = stat_summary["kruskal_p"].iloc[0]
        ax4.text(
            0.03,
            0.97,
            "Kruskal-Wallis P = {}".format(format_pvalue(kp)),
            transform=ax4.transAxes,
            ha="left",
            va="top",
            fontsize=8
        )

    ax4.text(-0.16, 1.08, "d", transform=ax4.transAxes, fontsize=14, fontweight="bold")

    plt.tight_layout(w_pad=2.4, h_pad=2.2)

    prefix = "enhancer_5group_snp_overlap.kruskal_dunn"

    out_pdf = os.path.join(FIGDIR, prefix + ".pdf")
    out_png = os.path.join(FIGDIR, prefix + ".png")
    out_svg = os.path.join(FIGDIR, prefix + ".svg")

    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()

    return out_pdf, out_png, out_svg


# ============================================================
# ============================================================

def main():
    if not os.path.exists(VCF):
        raise FileNotFoundError("VCF not found: {}".format(VCF))

    print("[INFO] Reading enhancer groups...")
    enh = read_enhancers(ENHANCER_FILE)

    print("[INFO] Enhancer rows:", enh.shape[0])

    interval_index = build_interval_index(enh)

    print("[INFO] Scanning VCF and overlapping SNPs with enhancers...")
    overlap_result = scan_vcf_and_overlap(
        vcf_path=VCF,
        interval_index=interval_index,
        enh=enh
    )

    print("[INFO] Total SNPs parsed from VCF:", overlap_result["total_snp_n"])

    print("[INFO] Building summary and enhancer burden table...")
    summary, burden, summary_file, burden_file = build_summary_and_burden(
        enh=enh,
        overlap_result=overlap_result
    )

    print("[INFO] Running Kruskal-Wallis + Dunn post hoc...")
    stat_summary, pairwise, pmat, stat_summary_file, pairwise_file, pmat_file = run_kruskal_dunn(
        burden=burden
    )

    print("[INFO] Plotting...")
    out_pdf, out_png, out_svg = plot_big_figure(
        summary=summary,
        burden=burden,
        stat_summary=stat_summary
    )

    print("\n .")
    print("Output directory:", OUTDIR)
    print("  enhancer file:", os.path.join(OUTDIR, "01_enhancers.used.no_merge.tsv"))
    print("  SNP BED:", overlap_result["all_snp_bed_file"])
    print("unique SNP overlap:", overlap_result["unique_overlap_file"])
    print("SNP-enhancer overlap details:", overlap_result["detail_overlap_file"])
    print("group summary:", summary_file)
    print("enhancer SNP burden:", burden_file)
    print("Kruskal + Dunn summary:", stat_summary_file)
    print("Dunn pairwise:", pairwise_file)
    print("Dunn adjusted-P matrix:", pmat_file)
    print("Figure PDF:", out_pdf)
    print("Figure PNG:", out_png)
    print("Figure SVG:", out_svg)


if __name__ == "__main__":
    main()
