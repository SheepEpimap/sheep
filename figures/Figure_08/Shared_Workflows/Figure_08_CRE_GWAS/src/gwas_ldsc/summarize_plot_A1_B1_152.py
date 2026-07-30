#!/usr/bin/env python
"""
A1 and B1 LDSC result summarization and heatmap plotting without coefficient columns.
"""
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(os.environ["GWAS_BASE"])
TRAIT_LIST = BASE / "all_152_traits.tsv"
SUMMARY_DIR = BASE / "summary"
FIG_DIR = BASE / "figures"
SUMMARY_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

ANNOT_ORDER = [
    "sfCRE_common", "sfCRE_intermediate", "sfCRE_tissue_specific",
    "sdCRE_common", "sdCRE_intermediate", "sdCRE_tissue_specific",
    "soCRE_common", "soCRE_intermediate", "soCRE_tissue_specific",
    "ssCRE_common", "ssCRE_intermediate", "ssCRE_tissue_specific",
    "sf_proj_common", "sf_proj_intermediate", "sf_proj_tissue_specific",
]

TISSUE_ORDER = [
    "Liver", "Spleen", "Heart", "Adipose", "Muscle",
    "Cortex", "Lung", "Ovary", "negative_control",
    "other_relevant",
]


def parse_one_results(results_file):
    df = pd.read_csv(results_file, sep="\t")
    return df.iloc[-1].to_dict()


def collect_all():
    records = []
    for source_dir in ["A1_152", "B1_152"]:
        for f in (BASE / "outputs" / source_dir).glob("*/*_baseline_custom.results"):
            annot_id = f.parent.name
            trait_id = f.name.replace(f"_{annot_id}_baseline_custom.results", "")
            row = parse_one_results(f)
            row["source"] = source_dir
            row["annot_id"] = annot_id
            row["trait_id"] = trait_id
            records.append(row)

    df = pd.DataFrame(records)
    if len(df) == 0:
        sys.exit("[ERROR] No .results files were found")

    meta = pd.read_csv(TRAIT_LIST, sep="\t")
    meta["trait_id"] = meta["filename"].str.replace(r"\.tsv\.bgz$", "", regex=True)
    df = df.merge(
        meta[["trait_id", "Description", "anchor_tissue"]],
        on="trait_id", how="left"
    )

    p = df["Enrichment_p"].clip(lower=1e-300)
    df["signed_log10p"] = -np.log10(p) * np.sign(df["Enrichment"] - 1)

    print(f"      Observed .results columns: {list(df.columns)}")
    return df


def make_pivot(df, value_col):
    trait_meta = (
        df[["trait_id", "Description", "anchor_tissue"]]
        .drop_duplicates()
        .copy()
    )
    trait_meta["tissue_rank"] = trait_meta["anchor_tissue"].apply(
        lambda x: TISSUE_ORDER.index(x) if x in TISSUE_ORDER else 99
    )
    trait_meta = trait_meta.sort_values(["tissue_rank", "Description"])
    trait_order = trait_meta["trait_id"].tolist()

    annot_present = [a for a in ANNOT_ORDER if a in df["annot_id"].unique()]

    pivot = df.pivot(index="trait_id", columns="annot_id", values=value_col)
    pivot = pivot.reindex(index=trait_order, columns=annot_present)

    desc_map = trait_meta.set_index("trait_id")["Description"].to_dict()
    tissue_map = trait_meta.set_index("trait_id")["anchor_tissue"].to_dict()
    pivot.index = [
        f"[{tissue_map[t][:3]}] {desc_map[t][:45]}"
        for t in pivot.index
    ]
    return pivot, trait_meta


def plot_heatmap(pivot, title, outfile, cmap, center,
                 vmin=None, vmax=None, cbar_label="", trait_meta=None):
    n_row, n_col = pivot.shape
    fig, ax = plt.subplots(figsize=(max(8, n_col * 0.8), max(8, n_row * 0.28)))
    sns.heatmap(
        pivot, ax=ax, cmap=cmap, center=center, vmin=vmin, vmax=vmax,
        cbar_kws={"label": cbar_label, "shrink": 0.5},
        linewidths=0.3, linecolor="white",
        xticklabels=True, yticklabels=True,
    )
    if trait_meta is not None:
        tm = trait_meta.copy()
        tm["tissue_rank"] = tm["anchor_tissue"].apply(
            lambda x: TISSUE_ORDER.index(x) if x in TISSUE_ORDER else 99
        )
        tm = tm.sort_values(["tissue_rank", "Description"]).reset_index(drop=True)
        change_points = tm.index[tm["anchor_tissue"] != tm["anchor_tissue"].shift()].tolist()
        for cp in change_points[1:]:
            ax.axhline(cp, color="black", lw=0.8)
    cre_breaks = [3, 6, 9, 12]
    for cb in cre_breaks:
        if cb < n_col:
            ax.axvline(cb, color="black", lw=0.8)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_title(title, fontsize=12, pad=10)
    ax.set_xlabel(""); ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  -> {outfile}")


def main():
    print("[1/3] Scanning .results files...")
    df = collect_all()
    print(f"      Total rows: {len(df)}; expected 615 (492 A1 + 123 B1)")

    # Retain only columns that are present.
    keep_cols = ["source", "trait_id", "Description", "anchor_tissue", "annot_id"]
    for c in ["Prop._SNPs", "Prop._h2", "Prop._h2_std_error",
              "Enrichment", "Enrichment_std_error", "Enrichment_p",
              "signed_log10p"]:
        if c in df.columns:
            keep_cols.append(c)
    # Retain coefficient columns when reruns used --print-coefficients.
    for c in ["Coefficient", "Coefficient_std_error", "Coefficient_z-score"]:
        if c in df.columns:
            keep_cols.append(c)

    df_long = df[keep_cols]
    long_csv = SUMMARY_DIR / "A1_B1_152_long.csv"
    df_long.to_csv(long_csv, index=False)
    print(f"  -> {long_csv}")

    print("[2/3] Building wide-format matrices...")
    pivot_e, trait_meta = make_pivot(df, "Enrichment")
    pivot_p, _ = make_pivot(df, "signed_log10p")
    pivot_e.to_csv(SUMMARY_DIR / "A1_B1_152_wide_enrich.csv")
    pivot_p.to_csv(SUMMARY_DIR / "A1_B1_152_wide_signedlog10p.csv")
    print(f"  -> {SUMMARY_DIR / 'A1_B1_152_wide_enrich.csv'}")
    print(f"  -> {SUMMARY_DIR / 'A1_B1_152_wide_signedlog10p.csv'}")

    print("[3/3] Drawing heatmaps...")
    enrich_max = np.nanpercentile(pivot_e.values, 95)
    plot_heatmap(
        pivot_e,
        title="Marginal enrichment (Prop._h2 / Prop._SNPs)",
        outfile=FIG_DIR / "heatmap_enrichment_152.png",
        cmap="RdBu_r", center=1,
        vmin=0, vmax=max(2, enrich_max),
        cbar_label="Enrichment",
        trait_meta=trait_meta,
    )
    plot_heatmap(
        pivot_p,
        title="Signed -log10(Enrichment_p)  (positive: Enrichment > 1)",
        outfile=FIG_DIR / "heatmap_signedlog10p_152.png",
        cmap="RdBu_r", center=0,
        vmin=-5, vmax=5,
        cbar_label="signed -log10(p)",
        trait_meta=trait_meta,
    )

    print("\nCompleted. Summary:")
    print(f"  CSV:   {SUMMARY_DIR}")
    print(f"  Figs:  {FIG_DIR}")


if __name__ == "__main__":
    main()
