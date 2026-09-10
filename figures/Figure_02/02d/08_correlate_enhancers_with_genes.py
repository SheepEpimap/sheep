#!/usr/bin/env python3
"""Link candidate enhancers to genes by cross-tissue signal correlation.

For each expressed gene, find enhancers whose midpoint is near its transcription
start site, correlate H3K27ac signal with RNA expression across shared samples,
and apply Benjamini-Hochberg correction to all tested pairs.

The fifth positional argument is retained for compatibility with the original
Figure 2d command line; pass ``-`` when no state-label file is required.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re

import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests


ENHANCER_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_expression_1", type=Path)
    parser.add_argument("gene_expression_2", type=Path)
    parser.add_argument("enhancer_h3k27ac_csv", type=Path)
    parser.add_argument("tss_bed", type=Path)
    parser.add_argument("state_labels", help="Legacy placeholder; use '-' if unused")
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("--window-bp", type=int, default=500_000)
    parser.add_argument("--min-dynamic-range", type=float, default=6.0)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args()


def dynamic_range_filter(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    dynamic_range = numeric.max(axis=1) / (numeric.min(axis=1) + 0.0001)
    return numeric.loc[dynamic_range > threshold]


def read_inputs(args: argparse.Namespace):
    genes_1 = pd.read_csv(args.gene_expression_1, sep="\t", index_col=0)
    genes_2 = pd.read_csv(args.gene_expression_2, sep="\t", index_col=0)
    genes = pd.concat([genes_1, genes_2], axis=1, join="inner")
    enhancers = pd.read_csv(args.enhancer_h3k27ac_csv, index_col=0)

    samples = [sample for sample in genes.columns if sample in enhancers.columns]
    if len(samples) < 3:
        raise ValueError("At least three shared samples are required")
    genes = dynamic_range_filter(genes.loc[:, samples], args.min_dynamic_range)
    enhancers = dynamic_range_filter(enhancers.loc[:, samples], args.min_dynamic_range)

    tss = pd.read_csv(
        args.tss_bed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "gene_id", "gene_name", "strand"],
    ).drop_duplicates("gene_id", keep="first").set_index("gene_id")
    return genes, enhancers, tss, samples


def enhancer_locations(enhancers: pd.DataFrame) -> dict[str, tuple[str, int]]:
    locations = {}
    for enhancer_id in enhancers.index:
        match = ENHANCER_RE.match(str(enhancer_id))
        if match is None:
            raise ValueError(f"Invalid enhancer coordinate: {enhancer_id!r}")
        midpoint = (int(match.group("start")) + int(match.group("end"))) // 2
        locations[str(enhancer_id)] = (match.group("chrom"), midpoint)
    return locations


def correlate_gene(gene_id, genes, enhancers, tss, locations, samples, window_bp):
    if gene_id not in tss.index:
        return []
    chrom = str(tss.at[gene_id, "chrom"])
    position = int(tss.at[gene_id, "start"])
    records = []
    for enhancer_id, (enhancer_chrom, midpoint) in locations.items():
        distance = abs(position - midpoint)
        if enhancer_chrom != chrom or distance > window_bp:
            continue
        paired = pd.concat(
            [genes.loc[gene_id, samples], enhancers.loc[enhancer_id, samples]],
            axis=1,
            keys=["expression", "h3k27ac"],
        ).dropna()
        if len(paired) < 3 or paired["expression"].nunique() < 2 or paired["h3k27ac"].nunique() < 2:
            continue
        correlation, p_value = pearsonr(paired["expression"], paired["h3k27ac"])
        records.append(
            {
                "gene": gene_id,
                "enhancer": enhancer_id,
                "pearson_r": correlation,
                "p_value": p_value,
                "distance_bp": distance,
                "n_samples": len(paired),
            }
        )
    return records


def run(args: argparse.Namespace) -> pd.DataFrame:
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    genes, enhancers, tss, samples = read_inputs(args)
    locations = enhancer_locations(enhancers)

    def analyse(gene_id):
        return correlate_gene(gene_id, genes, enhancers, tss, locations, samples, args.window_bp)

    if args.workers == 1:
        nested = [analyse(gene_id) for gene_id in genes.index]
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            nested = list(executor.map(analyse, genes.index))

    records = [record for gene_records in nested for record in gene_records]
    columns = ["gene", "enhancer", "pearson_r", "p_value", "distance_bp", "n_samples", "q_value"]
    if not records:
        return pd.DataFrame(columns=columns)
    result = pd.DataFrame.from_records(records)
    result["q_value"] = multipletests(result["p_value"], method="fdr_bh")[1]
    return result.loc[:, columns].sort_values(["gene", "enhancer"]).reset_index(drop=True)


def main() -> None:
    args = parse_args()
    result = run(args)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Tested {len(result)} enhancer-gene pairs; results written to {args.output_tsv}")


if __name__ == "__main__":
    main()
