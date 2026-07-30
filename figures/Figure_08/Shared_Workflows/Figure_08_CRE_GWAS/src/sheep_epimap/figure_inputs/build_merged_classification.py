#!/usr/bin/env python3
"""Build deduplicated E1-E9 hg38 CRE BEDs and the Figure 8d summary table."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

STATES = [f"E{i}" for i in range(1, 10)]
CRE_TYPES = ["sfCRE", "sdCRE", "soCRE", "ssCRE"]
ACTIVE_CRE_TYPES = ["sfCRE", "sdCRE", "soCRE"]
CRE_PRIORITY = {"sfCRE": 0, "sdCRE": 1, "soCRE": 2}
TISSUE_PAIRS = [
    ("Adipose", "adipose"),
    ("Colon", "colon"),
    ("Cortex", "cerebral-cortex"),
    ("Heart", "heart"),
    ("Liver", "liver"),
    ("Lung", "lung"),
    ("Muscle", "muscle"),
    ("Ovary", "ovary"),
    ("Sintest", "jejunum"),
    ("Spleen", "spleen"),
    ("Stomach", "abomasum"),
    ("Testis", "testis"),
]


def parse_coord_key(value: str) -> tuple[str, int, int] | None:
    parts = value.replace(":", "-").split("-")
    if len(parts) != 3:
        return None
    try:
        return parts[0], int(parts[1]), int(parts[2])
    except ValueError:
        return None


def iter_column_four(path: Path, skip_comments: bool = False):
    if not path.is_file() or path.stat().st_size == 0:
        return
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if skip_comments and line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) >= 4:
                yield columns[3]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_workdir = Path(os.environ["E1E9_WORKDIR"]) if os.environ.get("E1E9_WORKDIR") else None
    parser.add_argument("--workdir", type=Path, default=default_workdir)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.workdir is None:
        parser.error("--workdir or E1E9_WORKDIR is required")
    if args.output_dir is None:
        args.output_dir = args.workdir / "visualization_tables_E1E9_reuse_existing" / "results_hg38"
    return args


def main() -> int:
    args = parse_args()
    class_dir = args.workdir / "classified"
    lifted_dir = args.workdir / "lifted"
    if not class_dir.is_dir() or not lifted_dir.is_dir():
        raise FileNotFoundError(f"Expected classified/ and lifted/ under {args.workdir}")
    existing = list(args.output_dir.glob("E*.bed")) if args.output_dir.exists() else []
    summary = args.output_dir / "merged_classification_summary.tsv"
    if summary.exists():
        existing.append(summary)
    if existing and not args.force:
        raise FileExistsError(
            f"Refusing to overwrite {len(existing)} existing figure-input file(s) under {args.output_dir}"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, int | str]] = []
    for state in STATES:
        best: dict[str, str] = {}
        for human_tissue, sheep_tissue in TISSUE_PAIRS:
            for cre_type in ACTIVE_CRE_TYPES:
                path = class_dir / f"{human_tissue}_{sheep_tissue}_{state}.{cre_type}.bed"
                for key in iter_column_four(path) or ():
                    old = best.get(key)
                    if old is None or CRE_PRIORITY[cre_type] < CRE_PRIORITY[old]:
                        best[key] = cre_type

        ss_keys: set[str] = set()
        for human_tissue, _ in TISSUE_PAIRS:
            path = lifted_dir / f"{human_tissue}_{state}.unmapped"
            for key in iter_column_four(path, skip_comments=True) or ():
                if key not in best:
                    ss_keys.add(key)

        counts = {"sfCRE": 0, "sdCRE": 0, "soCRE": 0, "ssCRE": len(ss_keys)}
        for cre_type in best.values():
            counts[cre_type] += 1
        summary_rows.append({"state": state, **counts, "total": sum(counts.values())})

        for cre_type in CRE_TYPES:
            keys = ss_keys if cre_type == "ssCRE" else {key for key, value in best.items() if value == cre_type}
            output = args.output_dir / f"{state}.{cre_type}.bed"
            with output.open("w", encoding="utf-8", newline="") as handle:
                for key in sorted(keys):
                    coord = parse_coord_key(key)
                    if coord is not None:
                        handle.write(f"{coord[0]}\t{coord[1]}\t{coord[2]}\t{state}\n")

    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["state", *CRE_TYPES, "total"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
