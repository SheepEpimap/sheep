#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_GWAS_ROOT = Path("/vol2/00_panlab_rawdata/01_gwasdata/00_GWAS_result")
DEFAULT_REG_DIR = Path(
    "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/"
    "Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Permutation enrichment for GWAS STDERR tables against TSR E5 BED files."
    )
    parser.add_argument("--gwas-root", type=Path, default=DEFAULT_GWAS_ROOT)
    parser.add_argument("--reg-dir", type=Path, default=DEFAULT_REG_DIR)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nperm", type=int, default=10000)
    parser.add_argument(
        "--gwas-pattern",
        default="*/STDERR_*.tbl",
        help="Pattern under --gwas-root used to find GWAS tables.",
    )
    parser.add_argument(
        "--reg-pattern",
        default="TSR_*_E5_id.bed",
        help="Pattern under --reg-dir used to find BED files.",
    )
    return parser.parse_args()


def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)

    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = ranked[i] * n / rank
        prev = min(prev, value)
        adjusted[i] = min(prev, 1.0)

    result = np.empty(n, dtype=float)
    result[order] = adjusted
    return result


def load_gwas_table(path):
    df = pd.read_csv(path, sep="\t", dtype={"MarkerName": str})
    required = {"MarkerName", "Effect", "P-value"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    marker_parts = df["MarkerName"].str.extract(r"^(?P<chrom>[^_]+)_(?P<pos>\d+)$")
    keep = marker_parts["chrom"].notna() & marker_parts["pos"].notna()
    if not keep.any():
        raise ValueError(f"{path} does not contain parseable MarkerName values")

    df = df.loc[keep, ["MarkerName", "Effect", "P-value"]].copy()
    marker_parts = marker_parts.loc[keep].copy()
    df["chrom"] = "chr" + marker_parts["chrom"].astype(str)
    df["start"] = marker_parts["pos"].astype(np.int64)
    df["end"] = df["start"] + 1
    df["SNP"] = df["MarkerName"]
    df["b"] = pd.to_numeric(df["Effect"], errors="coerce")
    df["p_value"] = pd.to_numeric(df["P-value"], errors="coerce")
    df = df.dropna(subset=["b"]).copy()
    df["b2"] = df["b"].astype(float) ** 2
    df = df.sort_values(["chrom", "start", "SNP"]).reset_index(drop=True)
    return df[["chrom", "start", "end", "SNP", "b", "b2", "p_value"]]


def load_region_bed(path):
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2])
    df.columns = ["chrom", "start", "end"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"]).copy()
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    df = df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    return df


def build_overlap_mask(gwas_df, reg_df):
    mask = np.zeros(len(gwas_df), dtype=bool)

    for chrom, gwas_idx in gwas_df.groupby("chrom").groups.items():
        reg_chr = reg_df[reg_df["chrom"] == chrom]
        if reg_chr.empty:
            continue

        starts = reg_chr["start"].to_numpy()
        ends = reg_chr["end"].to_numpy()
        points = gwas_df.loc[gwas_idx, "start"].to_numpy()

        hit_idx = np.searchsorted(starts, points, side="right") - 1
        valid = hit_idx >= 0
        valid[valid] &= points[valid] < ends[hit_idx[valid]]
        mask[np.asarray(gwas_idx)] = valid

    return mask


def empirical_pvalue(stat, mask, nperm, rng):
    n = stat.size
    observed = stat[mask].sum()
    nset = int(mask.sum())
    if n == 0 or nset == 0:
        return observed, nset, 1.0

    shifts = rng.integers(0, n, size=nperm)
    perm_sums = np.empty(nperm, dtype=float)
    hit_indices = np.flatnonzero(mask)

    for i, shift in enumerate(shifts):
        shifted_idx = (hit_indices + shift) % n
        perm_sums[i] = stat[shifted_idx].sum()

    greater_equal = np.count_nonzero(perm_sums >= observed)
    p_value = (greater_equal + 1.0) / (nperm + 1.0)
    return observed, nset, p_value


def tissue_name_from_bed(path):
    name = path.name
    if not name.startswith("TSR_") or not name.endswith("_E5_id.bed"):
        return path.stem
    return name[len("TSR_") : -len("_E5_id.bed")]


def trait_name_from_gwas(path):
    return path.parent.name


def run_enrichment(gwas_paths, reg_paths, nperm):
    rng = np.random.default_rng(20260310)
    reg_cache = []
    for reg_path in reg_paths:
        reg_cache.append((tissue_name_from_bed(reg_path), load_region_bed(reg_path)))

    records = []
    for gwas_path in gwas_paths:
        trait = trait_name_from_gwas(gwas_path)
        gwas_df = load_gwas_table(gwas_path)
        stat = gwas_df["b2"].to_numpy()

        for tissue, reg_df in reg_cache:
            mask = build_overlap_mask(gwas_df, reg_df)
            observed, nset, p_value = empirical_pvalue(stat, mask, nperm=nperm, rng=rng)
            records.append(
                {
                    "trait": trait,
                    "tissue": tissue,
                    "nset": nset,
                    "sum_b2": observed,
                    "p": p_value,
                    "log10p": -math.log10(p_value) if p_value > 0 else np.inf,
                }
            )

    result = pd.DataFrame(records)
    result["p_fdr"] = bh_fdr(result["p"].to_numpy())
    result["significant_fdr_0.05"] = result["p_fdr"] < 0.05
    result["significant_p_0.05"] = result["p"] < 0.05
    return result.sort_values(["p_fdr", "p", "trait", "tissue"]).reset_index(drop=True)


def main():
    args = parse_args()
    gwas_paths = sorted(args.gwas_root.glob(args.gwas_pattern))
    reg_paths = sorted(args.reg_dir.glob(args.reg_pattern))

    if not gwas_paths:
        raise SystemExit(f"No GWAS files found under {args.gwas_root} with {args.gwas_pattern}")
    if not reg_paths:
        raise SystemExit(f"No BED files found under {args.reg_dir} with {args.reg_pattern}")

    result = run_enrichment(gwas_paths, reg_paths, nperm=args.nperm)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
