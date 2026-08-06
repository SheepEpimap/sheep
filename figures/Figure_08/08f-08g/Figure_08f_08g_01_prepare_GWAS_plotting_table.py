#!/usr/bin/env python
# Figure 8f-8g, step 01: prepare the LDSC plotting table.

import os
import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


BASE = Path(os.environ["GWAS_BASE"])
SUMMARY_DIR = BASE / "summary"
FIG_DIR = BASE / "figures"
SUMMARY_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

CLASSES = ["sfCRE", "sdCRE", "soCRE", "ssCRE", "sf_proj"]

TISSUE_ORDER = [
    "Liver", "Spleen", "Heart", "Adipose", "Muscle",
    "Cortex", "Lung", "Ovary",
    "Colon", "Sintest", "Stomach",
]


def boxplot_labeled(data, labels, **kwargs):
    """Compatible with both old and new matplotlib."""
    params = inspect.signature(plt.boxplot).parameters
    if "tick_labels" in params:
        return plt.boxplot(data, tick_labels=labels, **kwargs)
    else:
        return plt.boxplot(data, labels=labels, **kwargs)


df = pd.read_csv(SUMMARY_DIR / "A2_B2_cts_152_long.csv")

df = df[df["class"].isin(CLASSES)].copy()
df = df[df["tissue"].isin(TISSUE_ORDER)].copy()

df["P"] = df["Coefficient_P_value"].clip(lower=1e-300)
df["signed_log10p"] = -np.log10(df["P"]) * np.sign(df["Coefficient"])
df["signed_z"] = df["Coefficient"] / df["Coefficient_std_error"].replace(0, np.nan)
df["signed_z"] = df["signed_z"].replace([np.inf, -np.inf], np.nan)

# FDR 1: within each class across all trait × tissue tests
df["FDR_by_class"] = np.nan

for c in CLASSES:
    idx = df["class"] == c
    if idx.sum() > 0:
        df.loc[idx, "FDR_by_class"] = multipletests(
            df.loc[idx, "P"],
            method="fdr_bh"
        )[1]

# FDR 2: global across all trait × tissue × class tests
df["FDR_global"] = multipletests(df["P"], method="fdr_bh")[1]

df["sig_pos_by_class"] = (df["Coefficient"] > 0) & (df["FDR_by_class"] < 0.05)
df["sig_pos_global"] = (df["Coefficient"] > 0) & (df["FDR_global"] < 0.05)

df.to_csv(
    SUMMARY_DIR / "A2_B2_all_tissues_all_classes_enrichment_stats_long.csv",
    index=False
)

summary_class = (
    df.groupby("class")
      .agg(
          n_tests=("P", "size"),
          n_traits=("trait_id", "nunique"),
          n_tissues=("tissue", "nunique"),
          n_sig_pos_by_class=("sig_pos_by_class", "sum"),
          n_sig_pos_global=("sig_pos_global", "sum"),
          median_signed_z=("signed_z", "median"),
          mean_signed_z=("signed_z", "mean"),
          median_signed_log10p=("signed_log10p", "median"),
          mean_signed_log10p=("signed_log10p", "mean"),
      )
      .reindex(CLASSES)
      .reset_index()
)

summary_class["sig_pos_fraction_by_class"] = (
    summary_class["n_sig_pos_by_class"] / summary_class["n_tests"]
)

summary_class["sig_pos_fraction_global"] = (
    summary_class["n_sig_pos_global"] / summary_class["n_tests"]
)

summary_class.to_csv(
    SUMMARY_DIR / "A2_B2_all_tissues_all_classes_enrichment_summary_by_class.csv",
    index=False
)

summary_tissue_class = (
    df.groupby(["tissue", "class"])
      .agg(
          n_tests=("P", "size"),
          n_sig_pos_by_class=("sig_pos_by_class", "sum"),
          n_sig_pos_global=("sig_pos_global", "sum"),
          median_signed_z=("signed_z", "median"),
          mean_signed_z=("signed_z", "mean"),
          median_signed_log10p=("signed_log10p", "median"),
          mean_signed_log10p=("signed_log10p", "mean"),
      )
      .reset_index()
)

summary_tissue_class.to_csv(
    SUMMARY_DIR / "A2_B2_all_tissues_all_classes_enrichment_summary_by_tissue_class.csv",
    index=False
)

# Paired comparison across the same trait × tissue pairs
wide_z = df.pivot_table(
    index=["trait_id", "tissue"],
    columns="class",
    values="signed_z",
    aggfunc="first"
).reindex(columns=CLASSES)

tests = []
ref = "sfCRE"

for c in CLASSES:
    if c == ref:
        continue

    pair = wide_z[[ref, c]].dropna()
    diff = pair[c] - pair[ref]

    if len(diff) > 0 and diff.abs().sum() > 0:
        stat, p = wilcoxon(diff)
    else:
        stat, p = np.nan, np.nan

    tests.append({
        "comparison": f"{c} - {ref}",
        "class": c,
        "n_pair": len(diff),
        "median_diff_signed_z": diff.median(),
        "mean_diff_signed_z": diff.mean(),
        "wilcoxon_p": p,
    })

tests = pd.DataFrame(tests)
tests["wilcoxon_FDR"] = multipletests(
    tests["wilcoxon_p"].fillna(1),
    method="fdr_bh"
)[1]

tests.to_csv(
    SUMMARY_DIR / "A2_B2_all_tissues_all_classes_paired_tests_vs_sfCRE.csv",
    index=False
)

# Plot 1: significant positive CTS evidence count by class
plt.figure(figsize=(5.2, 4.0))
plt.bar(summary_class["class"], summary_class["n_sig_pos_by_class"])
plt.ylabel("Number of significant positive tests\n(coef > 0, class-wise FDR < 0.05)")
plt.xlabel("")
plt.xticks(rotation=35, ha="right")
plt.title("All tissues CTS enrichment evidence by CRE class")
plt.tight_layout()
plt.savefig(FIG_DIR / "A2_B2_all_tissues_sig_pos_count_by_class.pdf")
plt.close()

# Plot 2: signed Z-score distribution by class
plot_df = df.dropna(subset=["signed_z"]).copy()
plot_df["class"] = pd.Categorical(plot_df["class"], categories=CLASSES, ordered=True)
plot_df = plot_df.sort_values("class")

data = [
    plot_df.loc[plot_df["class"] == c, "signed_z"].values
    for c in CLASSES
]

plt.figure(figsize=(5.6, 4.2))
boxplot_labeled(data, CLASSES, showfliers=False)

rng = np.random.default_rng(1)

for i, c in enumerate(CLASSES, start=1):
    y = plot_df.loc[plot_df["class"] == c, "signed_z"].values
    x = i + rng.normal(0, 0.035, size=len(y))
    plt.scatter(x, y, s=6, alpha=0.25)

plt.axhline(0, lw=1)
plt.ylabel("Signed Z-score\n(Coefficient / SE)")
plt.xlabel("")
plt.xticks(rotation=35, ha="right")
plt.title("All trait × tissue CTS evidence by CRE class")
plt.tight_layout()
plt.savefig(FIG_DIR / "A2_B2_all_tissues_signed_z_by_class.pdf")
plt.close()

# Plot 3: paired difference relative to sfCRE
diff_rows = []

for c in CLASSES:
    if c == "sfCRE":
        continue

    pair = wide_z[["sfCRE", c]].dropna()

    tmp = pd.DataFrame({
        "trait_id": pair.index.get_level_values("trait_id"),
        "tissue": pair.index.get_level_values("tissue"),
        "class": c,
        "diff_vs_sfCRE": pair[c] - pair["sfCRE"],
    })

    diff_rows.append(tmp)

diff_df = pd.concat(diff_rows, ignore_index=True)

DIFF_CLASSES = ["sdCRE", "soCRE", "ssCRE", "sf_proj"]

diff_df["class"] = pd.Categorical(
    diff_df["class"],
    categories=DIFF_CLASSES,
    ordered=True
)

data = [
    diff_df.loc[diff_df["class"] == c, "diff_vs_sfCRE"].values
    for c in DIFF_CLASSES
]

plt.figure(figsize=(5.4, 4.2))
boxplot_labeled(data, DIFF_CLASSES, showfliers=False)

for i, c in enumerate(DIFF_CLASSES, start=1):
    y = diff_df.loc[diff_df["class"] == c, "diff_vs_sfCRE"].values
    x = i + rng.normal(0, 0.035, size=len(y))
    plt.scatter(x, y, s=6, alpha=0.25)

plt.axhline(0, lw=1)
plt.ylabel("Difference in signed Z-score\n(class - sfCRE)")
plt.xlabel("")
plt.xticks(rotation=35, ha="right")
plt.title("Paired CTS evidence difference across all tissues")
plt.tight_layout()
plt.savefig(FIG_DIR / "A2_B2_all_tissues_signed_z_diff_vs_sfCRE.pdf")
plt.close()

# Plot 4: tissue × class significant positive count heatmap
count_mat = summary_tissue_class.pivot(
    index="tissue",
    columns="class",
    values="n_sig_pos_by_class"
).reindex(index=TISSUE_ORDER, columns=CLASSES)

plt.figure(figsize=(6.2, 4.8))
plt.imshow(count_mat.values, aspect="auto")
plt.colorbar(label="Number of significant positive tests")
plt.xticks(range(len(CLASSES)), CLASSES, rotation=35, ha="right")
plt.yticks(range(len(TISSUE_ORDER)), TISSUE_ORDER)

for i in range(count_mat.shape[0]):
    for j in range(count_mat.shape[1]):
        v = count_mat.iloc[i, j]
        if pd.notna(v):
            plt.text(j, i, int(v), ha="center", va="center", fontsize=8)

plt.title("Significant positive CTS tests by tissue and CRE class")
plt.xlabel("")
plt.ylabel("")
plt.tight_layout()
plt.savefig(FIG_DIR / "A2_B2_all_tissues_sig_pos_count_tissue_class_heatmap.pdf")
plt.close()

# Plot 5: tissue × class median signed Z heatmap
z_mat = summary_tissue_class.pivot(
    index="tissue",
    columns="class",
    values="median_signed_z"
).reindex(index=TISSUE_ORDER, columns=CLASSES)

vmax = np.nanmax(np.abs(z_mat.values))
vmax = max(vmax, 1)

plt.figure(figsize=(6.2, 4.8))
plt.imshow(z_mat.values, aspect="auto", vmin=-vmax, vmax=vmax)
plt.colorbar(label="Median signed Z-score")
plt.xticks(range(len(CLASSES)), CLASSES, rotation=35, ha="right")
plt.yticks(range(len(TISSUE_ORDER)), TISSUE_ORDER)

for i in range(z_mat.shape[0]):
    for j in range(z_mat.shape[1]):
        v = z_mat.iloc[i, j]
        if pd.notna(v):
            plt.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8)

plt.title("Median CTS signed Z-score by tissue and CRE class")
plt.xlabel("")
plt.ylabel("")
plt.tight_layout()
plt.savefig(FIG_DIR / "A2_B2_all_tissues_median_signed_z_tissue_class_heatmap.pdf")
plt.close()

print("Done.")
print("\nClass-level summary:")
print(summary_class.to_string(index=False))

print("\nPaired tests vs sfCRE:")
print(tests.to_string(index=False))
