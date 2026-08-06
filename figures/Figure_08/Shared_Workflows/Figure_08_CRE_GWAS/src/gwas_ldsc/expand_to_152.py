#!/usr/bin/env python
"""
Derive a complete 152-trait table with the same columns as selected_41_traits.tsv.
Assign anchor_tissue automatically:
  liver_kidney             → Liver
  immune_blood             → Spleen
  metabolic_lipid          -> Liver
  cardiovascular           → Heart
  reproduction_endocrine   → Ovary
  brain_psych              → Cortex
  multi-tag containing respiratory -> Lung
  multi-tag uses the first recognized group
  other_relevant           -> other_relevant
"""
import argparse
import os
from pathlib import Path

import pandas as pd

GROUP_TO_ANCHOR = {
    "liver_kidney": "Liver",
    "immune_blood": "Spleen",
    "metabolic_lipid": "Liver",
    "cardiovascular": "Heart",
    "reproduction_endocrine": "Ovary",
    "brain_psych": "Cortex",
}


def assign_anchor(groups_str):
    parts = [p.strip() for p in groups_str.split(",")]
    # Use the first recognized group.
    for p in parts:
        if p in GROUP_TO_ANCHOR:
            return GROUP_TO_ANCHOR[p]
    # Assign multi-tag respiratory traits to Lung.
    if any(p == "respiratory" for p in parts):
        return "Lung"
    return "other_relevant"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inventory",
        type=Path,
        default=Path(os.environ["TRAIT_INVENTORY_152"]) if os.environ.get("TRAIT_INVENTORY_152") else None,
        help="Input 152-trait inventory TSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(os.environ["GWAS_BASE"]) / "all_152_traits.tsv" if os.environ.get("GWAS_BASE") else None,
        help="Output all_152_traits.tsv path.",
    )
    args = parser.parse_args()
    if args.inventory is None:
        parser.error("--inventory or TRAIT_INVENTORY_152 is required")
    if args.output is None:
        parser.error("--output or GWAS_BASE is required")
    return args


def main():
    args = parse_args()
    inv = pd.read_csv(args.inventory, sep="\t")
    inv["anchor_tissue"] = inv["groups"].apply(assign_anchor)

    out_cols = ["phenocode", "Description", "groups", "trait_type",
                "anchor_tissue", "filename", "sumstats_path"]
    out = inv[out_cols]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out)} traits to {args.output}")
    print("\nAnchor-tissue distribution:")
    print(out["anchor_tissue"].value_counts().to_string())


if __name__ == "__main__":
    main()
