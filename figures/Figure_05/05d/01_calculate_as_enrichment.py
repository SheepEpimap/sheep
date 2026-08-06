#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import re
import gzip
import tarfile
import tempfile
import subprocess
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_tissue_list(path: Path):
    with open(path) as f:
        tissues = [x.strip() for x in f if x.strip()]
    return tissues


def to_bool_series(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.lower()
    return s.isin(["true", "1", "yes", "y", "t"])


def motif_overlap_mask(df: pd.DataFrame, mode="either") -> pd.Series:
    """
    mode:
      - either: hits_overlap==True   hits_motifs  -
      - overlap_only:   hits_overlap==True
      - motifs_only:   hits_motifs  -
    """
    has_overlap_col = "hits_overlap" in df.columns
    has_motifs_col = "hits_motifs" in df.columns

    overlap = pd.Series(False, index=df.index)
    motifs = pd.Series(False, index=df.index)

    if has_overlap_col:
        overlap = to_bool_series(df["hits_overlap"])

    if has_motifs_col:
        m = df["hits_motifs"].astype(str).str.strip()
        motifs = ~(m.isin(["", "-", "nan", "NaN", "None", "none"]))

    if mode == "overlap_only":
        return overlap
    elif mode == "motifs_only":
        return motifs
    else:
        return overlap | motifs


def safe_div(a, b):
    if b == 0 or pd.isna(b):
        return np.nan
    return a / b


def find_vcf_member_in_tar(vcf_tar: Path):
    with tarfile.open(vcf_tar, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and re.search(r"\.vcf(\.gz)?$", m.name)]
        if not members:
            raise FileNotFoundError(f"{vcf_tar}   .vcf   .vcf.gz file")
        if len(members) > 1:
            print(f"[WARN] {vcf_tar}   VCF  , : {members[0].name}")
        return members[0].name


def count_total_vcf_sites_from_tar(vcf_tar: Path):
    member_name = find_vcf_member_in_tar(vcf_tar)

    n_sites = 0
    with tarfile.open(vcf_tar, "r:gz") as tar:
        f = tar.extractfile(member_name)
        if f is None:
            raise RuntimeError(f"  {vcf_tar}   {member_name}")

        if member_name.endswith(".gz"):
            with gzip.GzipFile(fileobj=f) as gz:
                for raw in gz:
                    if not raw.startswith(b"#"):
                        n_sites += 1
        else:
            for raw in f:
                if not raw.startswith(b"#"):
                    n_sites += 1

    return n_sites, member_name


def read_variant_5col(path: Path):
    """
    read 5  file:
    chr pos allele1 allele2 variant_id
    """
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] < 5:
        raise ValueError(f"{path}   5")
    df = df.iloc[:, :5].copy()
    df.columns = ["chr", "pos", "allele1", "allele2", "variant_id"]
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"]).copy()
    df["pos"] = df["pos"].astype(int)
    df = df.drop_duplicates(subset=["variant_id"]).copy()
    return df


def write_variant_bed(df: pd.DataFrame, out_bed: Path):
    """
      5   BED:
    chr start end variant_id
    SNP   1bp  :[pos-1, pos)
    """
    tmp = df[["chr", "pos", "variant_id"]].copy()
    tmp["start"] = tmp["pos"] - 1
    tmp["start"] = tmp["start"].clip(lower=0)
    tmp["end"] = tmp["pos"]
    tmp = tmp[["chr", "start", "end", "variant_id"]]
    tmp.to_csv(out_bed, sep="\t", header=False, index=False)


def count_bed_overlap_variants(variant_df: pd.DataFrame, motif_bed: Path):
    """
      bedtools intersect -u   variant_df   overlap motif_bed
      overlap   variant_id
    """
    if variant_df.empty:
        return 0

    if not motif_bed.exists() or motif_bed.stat().st_size == 0:
        return 0

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        var_bed = tmpdir / "variants.bed"
        ov_bed = tmpdir / "overlap.bed"

        write_variant_bed(variant_df, var_bed)

        cmd = [
            "bedtools", "intersect",
            "-u",
            "-a", str(var_bed),
            "-b", str(motif_bed)
        ]

        with open(ov_bed, "w") as fout:
            subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, text=True, check=True)

        if not ov_bed.exists() or ov_bed.stat().st_size == 0:
            return 0

        ov_df = pd.read_csv(ov_bed, sep="\t", header=None, dtype=str)
        if ov_df.shape[1] < 4:
            return 0

        return ov_df.iloc[:, 3].drop_duplicates().shape[0]


def main():
    parser = argparse.ArgumentParser(
        description="Calculate p1_t, q1_t, p2, p3 and their ratios using strict peak-filtered AS."
    )
    parser.add_argument(
        "--tissue-list",
        default="/data/home/sczd644/run/zsw_chrombpnet/tissue.txt",
        help=" file"
    )
    parser.add_argument(
        "--vcf-tar",
        default="/vol2/00_panlab_rawdata/01_gwasdata/WGS_235raw/35_high_depth/35_sub_indv.vcf.tar.gz",
        help="  VCF tar.gz file, statistics all VCF  "
    )
    parser.add_argument(
        "--as-in-peak-dir",
        default="/data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched",
        help="AS_in_peak filedirectory,file : ATAC_${tissue}.AS_in_peak.tsv"
    )
    parser.add_argument(
        "--nonas-in-peak-dir",
        default="/data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_peak",
        help="non-AS_in_peak filedirectory,file : ATAC_${tissue}_non_AS_in_peak.tsv"
    )
    parser.add_argument(
        "--ann-dir",
        default="/data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched/AS/02",
        help="AS  filedirectory,file : ${tissue}.annotations.tsv"
    )
    parser.add_argument(
        "--finemo-base",
        default="/data/home/sczd644/run/zsw_chrombpnet/finemo",
        help="finemo  directory,motif file : ${finemo-base}/${tissue}_finemo/hits_tf.bed"
    )
    parser.add_argument(
        "--logfc-pval-threshold",
        type=float,
        default=0.05,
        help="logfc.pval  ,  0.05"
    )
    parser.add_argument(
        "--motif-mode",
        default="either",
        choices=["either", "overlap_only", "motifs_only"],
        help="AS  motif overlap  "
    )
    parser.add_argument(
        "--outdir",
        default="/data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched/p1t_q1t_p2_p3_ratios_strict_peak",
        help="Output directory"
    )
    args = parser.parse_args()

    tissue_list = Path(args.tissue_list)
    vcf_tar = Path(args.vcf_tar)
    as_in_peak_dir = Path(args.as_in_peak_dir)
    nonas_in_peak_dir = Path(args.nonas_in_peak_dir)
    ann_dir = Path(args.ann_dir)
    finemo_base = Path(args.finemo_base)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # check bedtools
    cmd_check = subprocess.run(
        ["bash", "-lc", "command -v /public/home/zhangshiwen2/anaconda3/envs/abcmodel/bin/bedtools >/dev/null 2>&1"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    if cmd_check.returncode != 0:
        raise RuntimeError("  bedtools,  PATH")

    tissues = read_tissue_list(tissue_list)

    if not vcf_tar.exists():
        raise FileNotFoundError(f"VCF not found tar.gz file: {vcf_tar}")

    n_total_vcf_sites, used_member = count_total_vcf_sites_from_tar(vcf_tar)
    print(f"[INFO] VCF  = {n_total_vcf_sites}")
    print(f"[INFO]  VCF  = {used_member}")

    rows = []

    for tissue in tissues:
        as_in_peak_file = as_in_peak_dir / f"ATAC_{tissue}.AS_in_peak.tsv"
        nonas_in_peak_file = nonas_in_peak_dir / f"ATAC_{tissue}_non_AS_in_peak.tsv"
        ann_file = ann_dir / f"{tissue}.annotations.tsv"
        motif_bed = finemo_base / f"{tissue}_finemo" / "hits_tf.bed"

        if not as_in_peak_file.exists():
            print(f"[WARN]   {tissue}:   AS_in_peak file {as_in_peak_file}")
            continue
        if not nonas_in_peak_file.exists():
            print(f"[WARN]   {tissue}:   non-AS_in_peak file {nonas_in_peak_file}")
            continue
        if not ann_file.exists():
            print(f"[WARN]   {tissue}:  file {ann_file}")
            continue
        if not motif_bed.exists():
            print(f"[WARN]   {tissue}:   motif file {motif_bed}")
            continue

        # -------------------------
        # -------------------------
        as_peak_df = read_variant_5col(as_in_peak_file)
        n_as_total = len(as_peak_df)
        as_peak_ids = set(as_peak_df["variant_id"])

        # -------------------------
        # 2. test het = AS_in_peak ∪ nonAS_in_peak
        # -------------------------
        nonas_peak_df = read_variant_5col(nonas_in_peak_file)

        test_het_df = pd.concat([as_peak_df, nonas_peak_df], axis=0)
        test_het_df = test_het_df.drop_duplicates(subset=["variant_id"]).copy()

        n_as_in_peak = len(as_peak_df)
        n_nonas_in_peak = len(nonas_peak_df)
        n_test_het = len(test_het_df)

        # p1_t = test het / all VCF
        p1_t = safe_div(n_test_het, n_total_vcf_sites)

        # -------------------------
        # 3. q1_t = motif among test het
        # -------------------------
        n_test_het_motif = count_bed_overlap_variants(test_het_df, motif_bed)
        q1_t = safe_div(n_test_het_motif, n_test_het)

        # -------------------------
        # -------------------------
        ann_df = pd.read_csv(ann_file, sep="\t", dtype=str)
        required_cols = ["variant_id", "logfc.pval"]
        missing = [c for c in required_cols if c not in ann_df.columns]
        if missing:
            print(f"[WARN]   {tissue}: {ann_file}   {','.join(missing)}")
            continue

        ann_df = ann_df.drop_duplicates(subset=["variant_id"]).copy()
        ann_df = ann_df[ann_df["variant_id"].isin(as_peak_ids)].copy()

        n_as_annotated = len(ann_df)
        annotation_coverage_of_AS = safe_div(n_as_annotated, n_as_total)

        motif_mask = motif_overlap_mask(ann_df, mode=args.motif_mode)
        as_motif_df = ann_df.loc[motif_mask].copy()
        n_as_motif = len(as_motif_df)
        p2 = safe_div(n_as_motif, n_as_total)

        ann_df["logfc.pval_num"] = pd.to_numeric(ann_df["logfc.pval"], errors="coerce")
        sig_mask = ann_df["logfc.pval_num"].notna() & (ann_df["logfc.pval_num"] < args.logfc_pval_threshold)
        sig_df = ann_df.loc[sig_mask].copy()
        n_sig_as = len(sig_df)

        sig_motif_df = sig_df.loc[motif_overlap_mask(sig_df, mode=args.motif_mode)].copy()
        n_sig_as_motif = len(sig_motif_df)
        p3 = safe_div(n_sig_as_motif, n_sig_as)

        # -------------------------
        # -------------------------
        p2_over_p1_t = safe_div(p2, p1_t)
        p3_over_p1_t = safe_div(p3, p1_t)
        p2_over_q1_t = safe_div(p2, q1_t)
        p3_over_q1_t = safe_div(p3, q1_t)

        rows.append({
            "tissue": tissue,
            "vcf_member_used": used_member,

            "n_total_vcf_sites": n_total_vcf_sites,

            "n_AS_total_strictPeak": n_as_total,
            "n_AS_annotated_strictPeak": n_as_annotated,
            "annotation_coverage_of_AS_strictPeak": annotation_coverage_of_AS,

            "n_AS_in_peak": n_as_in_peak,
            "n_nonAS_in_peak": n_nonas_in_peak,
            "n_test_het": n_test_het,
            "p1_t_testHet_over_allVCF": p1_t,

            "n_test_het_motif": n_test_het_motif,
            "q1_t_motifAmongTestHet": q1_t,

            "n_AS_motif_overlap_strictPeak": n_as_motif,
            "p2_ASmotif_over_AS_strictPeak": p2,

            "n_sig_AS_logfc_pval_strictPeak": n_sig_as,
            "n_sig_AS_logfc_pval_and_motif_overlap_strictPeak": n_sig_as_motif,
            "p3_sigASmotif_over_sigAS_strictPeak": p3,

            "p2_over_p1_t": p2_over_p1_t,
            "p3_over_p1_t": p3_over_p1_t,
            "p2_over_q1_t": p2_over_q1_t,
            "p3_over_q1_t": p3_over_q1_t,

            "logfc_pval_threshold": args.logfc_pval_threshold,
            "motif_mode": args.motif_mode,
            "motif_bed": str(motif_bed),
        })

        print(
            f"[INFO] {tissue}: "
            f"AS_strictPeak={n_as_total}, "
            f"annotated_AS_strictPeak={n_as_annotated}, "
            f"test_het={n_test_het}, "
            f"test_het_motif={n_test_het_motif}, "
            f"AS_motif={n_as_motif}, "
            f"sigAS={n_sig_as}, "
            f"sigAS_motif={n_sig_as_motif}"
        )

    res_df = pd.DataFrame(rows)

    if res_df.empty:
        raise RuntimeError(" results, checkinput file .")

    res_df = res_df.sort_values(
        by=["p3_over_q1_t", "p2_over_q1_t", "tissue"],
        ascending=[False, False, True],
        na_position="last"
    ).reset_index(drop=True)

    out_file = outdir / "p1t_q1t_p2_p3_ratios_strictPeak.tsv"
    out_simple = outdir / "p1t_q1t_p2_p3_ratios_strictPeak.simple.tsv"

    res_df.to_csv(out_file, sep="\t", index=False)

    simple_cols = [
        "tissue",
        "n_total_vcf_sites",
        "n_test_het",
        "p1_t_testHet_over_allVCF",
        "n_test_het_motif",
        "q1_t_motifAmongTestHet",
        "n_AS_total_strictPeak",
        "n_AS_motif_overlap_strictPeak",
        "p2_ASmotif_over_AS_strictPeak",
        "n_sig_AS_logfc_pval_strictPeak",
        "n_sig_AS_logfc_pval_and_motif_overlap_strictPeak",
        "p3_sigASmotif_over_sigAS_strictPeak",
        "p2_over_p1_t",
        "p3_over_p1_t",
        "p2_over_q1_t",
        "p3_over_q1_t",
    ]
    res_df[simple_cols].to_csv(out_simple, sep="\t", index=False)

    print("[DONE]  ")
    print(f"[OUT] {out_file}")
    print(f"[OUT] {out_simple}")


if __name__ == "__main__":
    main()
