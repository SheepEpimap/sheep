#!/usr/bin/env python3
"""Filter Juicer samples and plot QC metrics for Figure 2j-2k inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STRICT = {
    "min_valid_pairs": 50_000_000,
    "min_valid_rate": 0.60,
    "max_dup_rate": 0.25,
    "min_cis_fraction": 0.55,
    "min_cis_trans_ratio": 1.20,
    "min_long_cis_fraction": 0.50,
}

RELAXED = {
    "min_valid_pairs": 30_000_000,
    "min_valid_rate": 0.60,
    "max_dup_rate": 0.30,
    "min_cis_fraction": 0.50,
    "min_cis_trans_ratio": 1.00,
    "min_long_cis_fraction": 0.35,
}

DISTANCE_ORDER = [
    "0_10kb",
    "10_50kb",
    "50_100kb",
    "100_500kb",
    "500kb_1Mb",
    "1_5Mb",
    "gt5Mb",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--basic-qc",
        type=Path,
        default=Path("juicer_basic_qc.tsv"),
    )
    parser.add_argument(
        "--distance-decay",
        type=Path,
        default=Path("juicer_distance_decay.tsv"),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("hic_qc_filter_out"),
    )
    return parser.parse_args()


def failure_reasons(row: pd.Series, limits: dict[str, float]) -> str:
    checks = [
        (
            "valid_pairs_merged_nodups",
            limits["min_valid_pairs"],
            "min",
            "valid_pairs",
        ),
        (
            "valid_rate_vs_raw",
            limits["min_valid_rate"],
            "min",
            "valid_rate",
        ),
        (
            "dup_rate_vs_valid_plus_dups",
            limits["max_dup_rate"],
            "max",
            "dup_rate",
        ),
        (
            "cis_fraction",
            limits["min_cis_fraction"],
            "min",
            "cis_fraction",
        ),
        (
            "cis_trans_ratio",
            limits["min_cis_trans_ratio"],
            "min",
            "cis_trans_ratio",
        ),
        (
            "long_cis_fraction_recalc",
            limits["min_long_cis_fraction"],
            "min",
            "long_cis_fraction",
        ),
    ]

    reasons: list[str] = []
    for column, cutoff, direction, label in checks:
        value = row[column]
        if pd.isna(value):
            reasons.append(f"{label}=NA")
        elif direction == "min" and value < cutoff:
            reasons.append(f"{label}<{cutoff}")
        elif direction == "max" and value > cutoff:
            reasons.append(f"{label}>{cutoff}")
    return ";".join(reasons)


def group_positions(groups: pd.Series) -> tuple[list[str], np.ndarray]:
    labels = list(pd.unique(groups))
    mapping = {label: index for index, label in enumerate(labels)}
    positions = np.array([mapping[label] for label in groups], dtype=float)
    jitter = np.random.default_rng(123).normal(0, 0.07, len(positions))
    return labels, positions + jitter


def plot_metric(
    frame: pd.DataFrame,
    metric: str,
    ylabel: str,
    output: Path,
    cutoff: float,
    log_scale: bool = False,
) -> None:
    plot_frame = frame.reset_index(drop=True)
    groups, x = group_positions(plot_frame["group"])
    colors = np.where(plot_frame["pass_strict"], "#3B82F6", "#DC2626")

    fig, axis = plt.subplots(figsize=(10, 5))
    axis.scatter(
        x,
        plot_frame[metric],
        c=colors,
        s=38,
        alpha=0.85,
        edgecolors="black",
        linewidths=0.3,
    )
    axis.axhline(cutoff, linestyle="--", linewidth=1, color="black")
    if log_scale:
        axis.set_yscale("log")
    axis.set_xticks(range(len(groups)))
    axis.set_xticklabels(groups, rotation=35, ha="right")
    axis.set_xlabel("Tissue group")
    axis.set_ylabel(ylabel)
    axis.grid(axis="y", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    basic = pd.read_csv(args.basic_qc, sep="\t")
    decay = pd.read_csv(args.distance_decay, sep="\t")

    required = {
        "sample",
        "group",
        "valid_pairs_merged_nodups",
        "valid_rate_vs_raw",
        "dup_rate_vs_valid_plus_dups",
        "cis_pairs",
        "cis_fraction",
        "cis_trans_ratio",
        "cis_short_lt20kb",
    }
    missing = sorted(required.difference(basic.columns))
    if missing:
        raise ValueError(f"Missing Juicer QC columns: {', '.join(missing)}")

    numeric = list(required.difference({"sample", "group"}))
    for column in numeric:
        basic[column] = pd.to_numeric(basic[column], errors="coerce")

    basic["valid_pairs_M"] = basic["valid_pairs_merged_nodups"] / 1e6
    basic["dup_rate_pct"] = basic["dup_rate_vs_valid_plus_dups"] * 100
    basic["cis_fraction_pct"] = basic["cis_fraction"] * 100
    basic["cis_long_ge20kb_recalc"] = (
        basic["cis_pairs"] - basic["cis_short_lt20kb"]
    ).where(lambda values: values >= 0)
    basic["long_cis_fraction_recalc"] = (
        basic["cis_long_ge20kb_recalc"] / basic["cis_pairs"]
    )
    basic["long_cis_fraction_pct"] = (
        basic["long_cis_fraction_recalc"] * 100
    )

    basic["strict_fail_reasons"] = basic.apply(
        failure_reasons,
        axis=1,
        limits=STRICT,
    )
    basic["relaxed_fail_reasons"] = basic.apply(
        failure_reasons,
        axis=1,
        limits=RELAXED,
    )
    basic["pass_strict"] = basic["strict_fail_reasons"].eq("")
    basic["pass_relaxed"] = basic["relaxed_fail_reasons"].eq("")
    basic["qc_recommendation"] = np.select(
        [basic["pass_strict"], basic["pass_relaxed"]],
        ["KEEP_STRICT", "KEEP_RELAXED"],
        default="DROP",
    )

    for column in ["n_pairs", "mean_distance"]:
        decay[column] = pd.to_numeric(decay[column], errors="coerce")
    decay["bin"] = pd.Categorical(
        decay["bin"],
        categories=DISTANCE_ORDER,
        ordered=True,
    )
    decay["total_cis_pairs"] = decay.groupby("sample")["n_pairs"].transform(
        "sum"
    )
    decay["pair_fraction"] = decay["n_pairs"] / decay["total_cis_pairs"]

    basic.to_csv(
        args.outdir / "hic_qc_with_filter.tsv",
        sep="\t",
        index=False,
    )
    decay.to_csv(
        args.outdir / "distance_decay_fraction.tsv",
        sep="\t",
        index=False,
    )

    for filename, mask in [
        ("hic_qc_keep_strict_samples.txt", basic["pass_strict"]),
        ("hic_qc_keep_relaxed_samples.txt", basic["pass_relaxed"]),
        ("hic_qc_drop_samples.txt", ~basic["pass_relaxed"]),
    ]:
        basic.loc[mask, "sample"].to_csv(
            args.outdir / filename,
            index=False,
            header=False,
        )

    group_summary = basic.groupby("group", as_index=False).agg(
        n_samples=("sample", "count"),
        n_keep_strict=("pass_strict", "sum"),
        n_keep_relaxed=("pass_relaxed", "sum"),
        median_valid_pairs_M=("valid_pairs_M", "median"),
        median_valid_rate=("valid_rate_vs_raw", "median"),
        median_dup_rate=("dup_rate_vs_valid_plus_dups", "median"),
        median_cis_fraction=("cis_fraction", "median"),
        median_long_cis_fraction=("long_cis_fraction_recalc", "median"),
    )
    group_summary.to_csv(
        args.outdir / "hic_qc_group_summary.tsv",
        sep="\t",
        index=False,
    )

    plot_metric(
        basic,
        "valid_pairs_M",
        "Valid read pairs (million)",
        args.outdir / "01_valid_pairs_by_group.pdf",
        STRICT["min_valid_pairs"] / 1e6,
        log_scale=True,
    )
    plot_metric(
        basic,
        "dup_rate_pct",
        "Duplicate rate (%)",
        args.outdir / "02_duplicate_rate_by_group.pdf",
        STRICT["max_dup_rate"] * 100,
    )
    plot_metric(
        basic,
        "cis_fraction_pct",
        "Cis fraction (%)",
        args.outdir / "03_cis_fraction_by_group.pdf",
        STRICT["min_cis_fraction"] * 100,
    )
    plot_metric(
        basic,
        "long_cis_fraction_pct",
        "Long-range cis fraction (>=20 kb, %)",
        args.outdir / "04_long_cis_fraction_by_group.pdf",
        STRICT["min_long_cis_fraction"] * 100,
    )

    fig, axis = plt.subplots(figsize=(10, 6))
    dropped = set(
        basic.loc[basic["qc_recommendation"].eq("DROP"), "sample"]
    )
    for sample, subset in decay.groupby("sample", observed=True):
        subset = subset.sort_values("bin")
        axis.plot(
            subset["bin"].astype(str),
            subset["pair_fraction"],
            marker="o",
            linewidth=0.8 if sample in dropped else 1.2,
            linestyle="--" if sample in dropped else "-",
            alpha=0.35 if sample in dropped else 0.85,
        )
    axis.set_yscale("log")
    axis.set_xlabel("Genomic distance bin")
    axis.set_ylabel("Fraction of cis pairs")
    axis.tick_params(axis="x", rotation=35)
    axis.grid(True, linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(args.outdir / "05_distance_decay.pdf")
    plt.close(fig)

    print(basic["qc_recommendation"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
