#!/usr/bin/env python
"""
Extract total h2, mean_chi2 and intercept for each trait from LDSC log files.
Build the trait-QC table used to select traits for final display.

Input: outputs/A1_152/sfCRE_common/<trait>_<annot>_baseline_custom.log
Output: summary/trait_qc_152.tsv
"""
import re
import os
import pandas as pd
from pathlib import Path

BASE = Path(os.environ["GWAS_BASE"])
LOG_DIR = BASE / "outputs" / "A1_152" / "sfCRE_common"
TRAIT_LIST = BASE / "all_152_traits.tsv"
OUT = BASE / "summary" / "trait_qc_152.tsv"

PATTERNS = {
    "h2": r"Total Observed scale h2:\s+([\-\d\.eE]+)\s+\(([\-\d\.eE]+)\)",
    "lambda_gc": r"Lambda GC:\s+([\-\d\.eE]+)",
    "mean_chi2": r"Mean Chi\^2:\s+([\-\d\.eE]+)",
    "intercept": r"Intercept:\s+([\-\d\.eE]+)\s+\(([\-\d\.eE]+)\)",
}


def parse_log(log_path):
    text = log_path.read_text()
    out = {}
    m = re.search(PATTERNS["h2"], text)
    if not m: return None
    out["h2"] = float(m.group(1))
    out["h2_se"] = float(m.group(2))
    out["h2_z"] = out["h2"] / out["h2_se"] if out["h2_se"] > 0 else float("nan")
    m = re.search(PATTERNS["mean_chi2"], text)
    out["mean_chi2"] = float(m.group(1)) if m else float("nan")
    m = re.search(PATTERNS["intercept"], text)
    if m:
        out["intercept"] = float(m.group(1))
        out["intercept_se"] = float(m.group(2))
    else:
        out["intercept"] = float("nan")
        out["intercept_se"] = float("nan")
    return out


def main():
    meta = pd.read_csv(TRAIT_LIST, sep="\t")
    meta["trait_id"] = meta["filename"].str.replace(r"\.tsv\.bgz$", "", regex=True)

    rows = []
    for _, r in meta.iterrows():
        log = LOG_DIR / f"{r['trait_id']}_sfCRE_common_baseline_custom.log"
        if not log.exists():
            rows.append({"trait_id": r["trait_id"], **{k: float("nan") for k in
                ["h2", "h2_se", "h2_z", "mean_chi2", "intercept", "intercept_se"]}})
            continue
        stats = parse_log(log)
        rows.append({"trait_id": r["trait_id"], **(stats or {})})

    qc = pd.DataFrame(rows)
    qc = qc.merge(meta[["trait_id", "Description", "anchor_tissue", "groups"]],
                  on="trait_id", how="left")

    # Three QC flags; True indicates that the criterion passed.
    qc["pass_h2_z"] = qc["h2_z"] > 4
    qc["pass_mean_chi2"] = qc["mean_chi2"] > 1.02
    qc["pass_intercept"] = (qc["intercept"] > 0.9) & (qc["intercept"] < 1.2)
    qc["pass_qc"] = qc["pass_h2_z"] & qc["pass_mean_chi2"] & qc["pass_intercept"]

    # Determine whether an expected anchor tissue is defined.
    anchored_set = {"Liver", "Spleen", "Heart", "Adipose", "Muscle",
                    "Cortex", "Lung", "Ovary"}
    qc["pass_anchor"] = qc["anchor_tissue"].isin(anchored_set)

    qc["pass_final"] = qc["pass_qc"] & qc["pass_anchor"]

    qc = qc.sort_values(["pass_final", "anchor_tissue", "Description"],
                        ascending=[False, True, True])
    qc.to_csv(OUT, sep="\t", index=False)
    print(f"Output: {OUT}")

    print("\nQC summary for 152 traits:")
    print(f"  Passed h2_z > 4: {qc['pass_h2_z'].sum()}")
    print(f"  Passed mean_chi2 > 1.02: {qc['pass_mean_chi2'].sum()}")
    print(f"  Passed intercept in [0.9, 1.2]: {qc['pass_intercept'].sum()}")
    print(f"  Passed all three QC criteria: {qc['pass_qc'].sum()}")
    print(f"  Passed QC with a defined anchor: {qc['pass_final'].sum()}")
    print("\nAnchor distribution among final passing traits:")
    print(qc[qc["pass_final"]]["anchor_tissue"].value_counts().to_string())


if __name__ == "__main__":
    main()
