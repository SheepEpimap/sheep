#!/usr/bin/env python
"""
Diagnose sfCRE-by-tissue tests for 24 candidate traits.
Determine whether traits lack signal or fail because the rules are too stringent.
"""
import pandas as pd
import numpy as np
import os
from pathlib import Path
from statsmodels.stats.multitest import multipletests

BASE = Path(os.environ["GWAS_BASE"])
qc = pd.read_csv(BASE / "summary" / "trait_qc_152.tsv", sep="\t")
df = pd.read_csv(BASE / "summary" / "A2_B2_cts_152_long.csv")

ANCHORED = {"Liver", "Spleen", "Heart", "Adipose", "Muscle",
            "Cortex", "Lung", "Ovary"}

pool = qc[qc["pass_qc"] & ~qc["anchor_tissue"].isin(ANCHORED)]

diag = []
for _, row in pool.iterrows():
    tid = row["trait_id"]
    sub = df[(df["trait_id"] == tid) & (df["class"] == "sfCRE")]
    if len(sub) == 0:
        continue
    p = sub["Coefficient_P_value"].values
    rej_05, q, _, _ = multipletests(p, alpha=0.05, method="fdr_bh")
    rej_10, _, _, _ = multipletests(p, alpha=0.10, method="fdr_bh")
    best_idx = p.argmin()
    diag.append({
        "trait_id": tid,
        "Description": row["Description"],
        "h2_z": row["h2_z"],
        "min_raw_p": p.min(),
        "min_BH_q": q.min(),
        "best_tissue_sfCRE": sub.iloc[best_idx]["tissue"],
        "best_signed_log10p": sub.iloc[best_idx]["signed_log10p"],
        "n_sig_q05": int(rej_05.sum()),
        "n_sig_q10": int(rej_10.sum()),
        "n_nominal_p05": int((p < 0.05).sum()),
    })

diag_df = pd.DataFrame(diag).sort_values("min_BH_q")
out = BASE / "summary" / "trait_data_driven_diagnostic.tsv"
diag_df.to_csv(out, sep="\t", index=False)
print(f"Output: {out}\n")

print("=== Summary of 24 candidate traits, ordered by min_BH_q ===")
print(diag_df.to_string(index=False))

print("\n=== Distribution statistics ===")
print(f"min_BH_q < 0.05: {(diag_df['min_BH_q'] < 0.05).sum()} trait")
print(f"min_BH_q < 0.10: {(diag_df['min_BH_q'] < 0.10).sum()} trait")
print(f"min_BH_q < 0.20: {(diag_df['min_BH_q'] < 0.20).sum()} trait")
print(f"Any raw P < 0.05: {(diag_df['min_raw_p'] < 0.05).sum()} traits")
