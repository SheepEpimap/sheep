#!/usr/bin/env python
"""
Identify tissue-specific candidate hits among QC-passing traits without a defined anchor.

Pre-specified rules:
  1. The trait must pass QC: h2_z > 4, mean_chi2 > 1.02 and intercept in [0.9, 1.2].
  2. Apply BH-FDR correction to Coefficient_P_value across 11 tissues within sfCRE.
  3. BH q < 0.05 must occur in only one or two tissues to enforce specificity.
  4. For each candidate tissue, sfCRE signed_log10p must exceed the maximum
     among sdCRE, soCRE and ssCRE.

Output:
  summary/trait_data_driven_hits.tsv  - candidate trait + assigned anchor
"""
import pandas as pd
import numpy as np
import os
from pathlib import Path
from statsmodels.stats.multitest import multipletests

BASE = Path(os.environ["GWAS_BASE"])
QC_FILE = BASE / "summary" / "trait_qc_152.tsv"
LONG_FILE = BASE / "summary" / "A2_B2_cts_152_long.csv"
OUT = BASE / "summary" / "trait_data_driven_hits.tsv"

ANCHORED = {"Liver", "Spleen", "Heart", "Adipose", "Muscle",
            "Cortex", "Lung", "Ovary"}
MAX_TISSUE_HITS = 2
FDR_Q = 0.05


def main():
    qc = pd.read_csv(QC_FILE, sep="\t")
    df = pd.read_csv(LONG_FILE)

    # Candidate pool: QC-passing traits without an anchored tissue.
    pool = qc[qc["pass_qc"] & ~qc["anchor_tissue"].isin(ANCHORED)]
    print(f"Candidate pool: {len(pool)} traits passed QC without a defined anchor")

    hits = []
    for _, row in pool.iterrows():
        tid = row["trait_id"]
        sub = df[(df["trait_id"] == tid) & (df["class"] == "sfCRE")]
        if len(sub) == 0:
            continue

        # Apply BH-FDR across 11 tissues for one trait.
        p = sub["Coefficient_P_value"].values
        rej, q, _, _ = multipletests(p, alpha=FDR_Q, method="fdr_bh")
        sub = sub.copy()
        sub["q_value"] = q
        sub["sig"] = rej

        n_sig = sub["sig"].sum()
        if n_sig == 0 or n_sig > MAX_TISSUE_HITS:
            continue

        # Validate the conservation gradient for each candidate tissue.
        for _, sig_row in sub[sub["sig"]].iterrows():
            t = sig_row["tissue"]
            sf_signal = sig_row["signed_log10p"]
            other = df[(df["trait_id"] == tid) &
                       (df["tissue"] == t) &
                       (df["class"].isin(["sdCRE", "soCRE", "ssCRE"]))]
            if len(other) == 0:
                continue
            max_other = other["signed_log10p"].max()
            if sf_signal <= max_other:
                continue

            hits.append({
                "trait_id": tid,
                "Description": row["Description"],
                "original_anchor": row["anchor_tissue"],
                "candidate_tissue": t,
                "sfCRE_signed_log10p": sf_signal,
                "BH_q": sig_row["q_value"],
                "max_other_class_signal": max_other,
                "h2_z": row["h2_z"],
                "mean_chi2": row["mean_chi2"],
            })

    out_df = pd.DataFrame(hits).sort_values(
        ["candidate_tissue", "BH_q"]
    )
    out_df.to_csv(OUT, sep="\t", index=False)
    print(f"\nFound {len(out_df)} candidate hits from {out_df['trait_id'].nunique()} unique traits")
    print("\nCandidate-tissue distribution:")
    print(out_df["candidate_tissue"].value_counts().to_string())
    print(f"\nOutput: {OUT}")
    print("\nManual review: inspect each candidate Description and confirm that its")
    print("  candidate tissue is biologically plausible. Treat unexpected mappings")
    print("  such as hair colour to cortex or dental traits to stomach with caution.")


if __name__ == "__main__":
    main()
