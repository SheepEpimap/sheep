#!/usr/bin/env python3
"""Run and verify the synthetic enhancer-target gene analysis."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[2]
DEMO = ROOT / "tests" / "demo"
SCRIPT = ROOT / "figures" / "Figure_02" / "02d" / "08_correlate_enhancers_with_genes.py"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify(observed_path: Path) -> None:
    expected = read_rows(DEMO / "expected_output" / "enhancer_gene_correlations.tsv")
    observed = read_rows(observed_path)
    if len(observed) != len(expected):
        raise AssertionError(f"Expected {len(expected)} pairs, observed {len(observed)}")
    for observed_row, expected_row in zip(observed, expected):
        for column in ("gene", "enhancer", "distance_bp", "n_samples"):
            if observed_row[column] != expected_row[column]:
                raise AssertionError(f"Mismatch in {column}: {observed_row[column]} != {expected_row[column]}")
        for column in ("pearson_r", "p_value", "q_value"):
            if not math.isclose(float(observed_row[column]), float(expected_row[column]), abs_tol=1e-12):
                raise AssertionError(f"Mismatch in {column}: {observed_row[column]} != {expected_row[column]}")


def main() -> None:
    start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="sheep_epimap_demo_") as temporary_directory:
        output = Path(temporary_directory) / "enhancer_gene_correlations.tsv"
        command = [
            sys.executable,
            str(SCRIPT),
            str(DEMO / "input" / "gene_expression_part1.tsv"),
            str(DEMO / "input" / "gene_expression_part2.tsv"),
            str(DEMO / "input" / "enhancer_h3k27ac.csv"),
            str(DEMO / "input" / "gene_tss.bed"),
            "-",
            str(output),
            "--workers",
            "1",
        ]
        subprocess.run(command, check=True)
        verify(output)
    elapsed = time.perf_counter() - start
    print(f"PASS: enhancer-target gene demo completed in {elapsed:.3f} seconds")


if __name__ == "__main__":
    main()
