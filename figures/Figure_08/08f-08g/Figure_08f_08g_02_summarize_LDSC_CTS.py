#!/usr/bin/env python
# Figure 8f-8g, step 02: summarize LDSC CTS results.
"""
A2 and B2 LDSC CTS result summarization and heatmap plotting.

Inputs:
  outputs/A2_152/<class>/<trait>_<class>_cts.cell_type_results.txt
  outputs/B2_152/sf_proj/<trait>_sf_proj_cts.cell_type_results.txt
  all_152_traits.tsv

LDSC CTS output columns:
  Name, Coefficient, Coefficient_std_error, Coefficient_P_value

Interpretation:
  - P_value is one-tailed (H1: coefficient > 0).
  - A positive coefficient indicates that the tissue carries more SNP heritability
    than the union of all tissues from the same CRE class.
  - Raw coefficients should not be compared directly across tissues or classes.
  - The heatmaps primarily show LDSC enrichment evidence patterns.

Outputs:
  summary/A2_B2_cts_152_long.csv
  summary/A2_B2_cts_152_<class>_signedlog10p.csv
  summary/A2_B2_cts_152_anchor_matched.csv
  summary/A2_B2_cts_152_alltissue_allclass_signedlog10p.csv
  summary/A2_B2_cts_152_alltissue_allclass_*_for_cluster.csv
  summary/A2_B2_cts_152_clustered_*_trait_order.csv
  summary/A2_B2_cts_152_clustered_*_tissue_class_order.csv

  figures/cts_<class>_heatmap_152.png
  figures/cts_anchor_matched_152.png
  figures/cts_all_tissues_all_classes_152.png
  figures/cts_all_tissues_all_classes_152_clustered_signed_log10p.png
  figures/cts_all_tissues_all_classes_152_clustered_signed_z.png
  figures/cts_all_tissues_all_classes_152_clustered_signed_z_rowscaled.png
"""

import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


BASE = Path(os.environ["GWAS_BASE"])
TRAIT_LIST = BASE / "all_152_traits.tsv"
SUMMARY_DIR = BASE / "summary"
FIG_DIR = BASE / "figures"
SUMMARY_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

CLASSES = ["sfCRE", "sdCRE", "soCRE", "ssCRE", "sf_proj"]

# Tissue order on the CTS x axis.
TISSUE_ORDER = [
    "Liver", "Spleen", "Heart", "Adipose", "Muscle",
    "Cortex", "Lung", "Ovary",
    "Colon", "Sintest", "Stomach",
]

# Anchor-tissue order used to group traits on the y axis.
TRAIT_ANCHOR_ORDER = [
    "Liver", "Spleen", "Heart", "Adipose", "Muscle",
    "Cortex", "Lung", "Ovary",
    "negative_control",
    "other_relevant",
]


def safe_short_label(trait_id, desc_map, tissue_map, max_len=42):
    """Build robust trait labels when Description or anchor_tissue is missing."""
    desc = desc_map.get(trait_id, trait_id)
    anchor = tissue_map.get(trait_id, "NA")

    if pd.isna(desc):
        desc = trait_id
    if pd.isna(anchor):
        anchor = "NA"

    return f"[{str(anchor)[:3]}] {str(desc)[:max_len]}"


def collect_all():
    records = []

    for c in CLASSES:
        if c == "sf_proj":
            search_dir = BASE / "outputs" / "B2_152" / "sf_proj"
        else:
            search_dir = BASE / "outputs" / "A2_152" / c

        if not search_dir.exists():
            print(f"  [WARN] Directory not found; skipping: {search_dir}")
            continue

        files = sorted(search_dir.glob("*_cts.cell_type_results.txt"))
        if len(files) == 0:
            print(f"  [WARN] No result files found: {search_dir}")
            continue

        for f in files:
            fname = f.name.replace(".cell_type_results.txt", "")
            trait_id = fname.replace(f"_{c}_cts", "")

            try:
                tdf = pd.read_csv(f, sep="\t")
            except Exception as e:
                print(f"  [WARN] Read failed; skipping: {f}")
                print(f"         {e}")
                continue

            required_cols = {
                "Name",
                "Coefficient",
                "Coefficient_std_error",
                "Coefficient_P_value",
            }
            missing = required_cols - set(tdf.columns)
            if missing:
                print(f"  [WARN] Missing columns {missing}; skipping: {f}")
                continue

            for _, r in tdf.iterrows():
                records.append({
                    "class": c,
                    "trait_id": trait_id,
                    "tissue": r["Name"],
                    "Coefficient": r["Coefficient"],
                    "Coefficient_std_error": r["Coefficient_std_error"],
                    "Coefficient_P_value": r["Coefficient_P_value"],
                })

    df = pd.DataFrame(records)

    if len(df) == 0:
        sys.exit("[ERROR] No .cell_type_results.txt files were found")

    # Add trait metadata.
    meta = pd.read_csv(TRAIT_LIST, sep="\t")
    meta["trait_id"] = meta["filename"].str.replace(r"\.tsv\.bgz$", "", regex=True)

    df = df.merge(
        meta[["trait_id", "Description", "anchor_tissue"]],
        on="trait_id",
        how="left",
    )

    # signed -log10(P)
    p = df["Coefficient_P_value"].clip(lower=1e-300)
    df["neg_log10p"] = -np.log10(p)
    df["signed_log10p"] = df["neg_log10p"] * np.sign(df["Coefficient"])

    # signed Z-score
    # This is more suitable for clustering than the raw coefficient because it is standardized by SE.
    se = df["Coefficient_std_error"].replace(0, np.nan)
    df["signed_z"] = df["Coefficient"] / se
    df["signed_z"] = df["signed_z"].replace([np.inf, -np.inf], np.nan)

    return df


def order_traits(df):
    """Order traits by anchor tissue and then by Description within each group."""
    tm = (
        df[["trait_id", "Description", "anchor_tissue"]]
        .drop_duplicates()
        .copy()
    )

    tm["rank"] = tm["anchor_tissue"].apply(
        lambda x: TRAIT_ANCHOR_ORDER.index(x) if x in TRAIT_ANCHOR_ORDER else 99
    )

    tm = tm.sort_values(["rank", "Description"]).reset_index(drop=True)
    return tm


def plot_one_class_heatmap(df, c, trait_meta):
    sub = df[df["class"] == c].copy()

    if len(sub) == 0:
        print(f"  [WARN] No data for {c}")
        return

    pivot = sub.pivot_table(
        index="trait_id",
        columns="tissue",
        values="signed_log10p",
        aggfunc="first",
    )

    tissues_present = [t for t in TISSUE_ORDER if t in pivot.columns]

    pivot = pivot.reindex(
        index=trait_meta["trait_id"].tolist(),
        columns=tissues_present,
    )

    desc_map = trait_meta.set_index("trait_id")["Description"].to_dict()
    tissue_map = trait_meta.set_index("trait_id")["anchor_tissue"].to_dict()

    pivot.index = [
        safe_short_label(t, desc_map, tissue_map, max_len=42)
        for t in pivot.index
    ]

    pivot.to_csv(SUMMARY_DIR / f"A2_B2_cts_152_{c}_signedlog10p.csv")

    n_row, n_col = pivot.shape
    fig, ax = plt.subplots(
        figsize=(max(7, n_col * 0.85), max(8, n_row * 0.28))
    )

    sns.heatmap(
        pivot,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        cbar_kws={"label": "signed -log10(P)", "shrink": 0.5},
        linewidths=0.3,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
    )

    # Highlight tissue-matched cells.
    trait_ids = trait_meta["trait_id"].tolist()
    raw_tissue_map = trait_meta.set_index("trait_id")["anchor_tissue"].to_dict()

    for i, trait_id in enumerate(trait_ids):
        anchor = raw_tissue_map.get(trait_id, None)
        if anchor in tissues_present:
            j = tissues_present.index(anchor)
            rect = mpatches.Rectangle(
                (j, i),
                1,
                1,
                linewidth=2.0,
                edgecolor="black",
                facecolor="none",
            )
            ax.add_patch(rect)

    tm = trait_meta.copy()
    change_points = tm.index[
        tm["anchor_tissue"] != tm["anchor_tissue"].shift()
    ].tolist()

    for cp in change_points[1:]:
        ax.axhline(cp, color="gray", lw=0.8)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)

    ax.set_title(
        f"CTS conditional enrichment: {c}\n"
        f"black box = tissue-matched annotation",
        fontsize=12,
        pad=10,
    )
    ax.set_xlabel("focal tissue (vs same-class union)", fontsize=10)
    ax.set_ylabel("")

    plt.tight_layout()
    out = FIG_DIR / f"cts_{c}_heatmap_152.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"  -> {out}")


def plot_anchor_matched_summary(df, trait_meta):
    """
    For each trait, extract signed_log10p for the anchor tissue across five classes.
    Return a trait-by-class matrix.
    """
    rows = []

    tissue_map = trait_meta.set_index("trait_id")["anchor_tissue"].to_dict()
    desc_map = trait_meta.set_index("trait_id")["Description"].to_dict()

    for trait_id in trait_meta["trait_id"]:
        anchor = tissue_map.get(trait_id, None)

        if anchor in {"negative_control", "other_relevant"}:
            continue

        for c in CLASSES:
            sub = df[
                (df["class"] == c) &
                (df["trait_id"] == trait_id) &
                (df["tissue"] == anchor)
            ]

            v = sub["signed_log10p"].iloc[0] if len(sub) else np.nan

            rows.append({
                "trait_id": trait_id,
                "Description": desc_map.get(trait_id, trait_id),
                "anchor_tissue": anchor,
                "class": c,
                "signed_log10p": v,
            })

    sm = pd.DataFrame(rows)

    pivot = sm.pivot(
        index="trait_id",
        columns="class",
        values="signed_log10p",
    )

    pivot = pivot.reindex(columns=CLASSES)

    tm_sub = trait_meta[
        ~trait_meta["anchor_tissue"].isin(["negative_control", "other_relevant"])
    ].copy()

    pivot = pivot.reindex(index=tm_sub["trait_id"].tolist())

    desc_map = tm_sub.set_index("trait_id")["Description"].to_dict()
    tissue_map = tm_sub.set_index("trait_id")["anchor_tissue"].to_dict()

    pivot.index = [
        safe_short_label(t, desc_map, tissue_map, max_len=42)
        for t in pivot.index
    ]

    pivot.to_csv(SUMMARY_DIR / "A2_B2_cts_152_anchor_matched.csv")

    n_row, n_col = pivot.shape
    fig, ax = plt.subplots(figsize=(6, max(8, n_row * 0.30)))

    sns.heatmap(
        pivot,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        cbar_kws={"label": "signed -log10(P) at anchor tissue", "shrink": 0.5},
        linewidths=0.3,
        linecolor="white",
    )

    tm_idx = tm_sub.reset_index(drop=True)
    cp = tm_idx.index[
        tm_idx["anchor_tissue"] != tm_idx["anchor_tissue"].shift()
    ].tolist()

    for c_ in cp[1:]:
        ax.axhline(c_, color="gray", lw=0.8)

    # Separate the four A2 classes from B2 sf_proj.
    if "sf_proj" in CLASSES:
        sf_proj_pos = CLASSES.index("sf_proj")
        ax.axvline(sf_proj_pos, color="black", lw=1.0)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)

    ax.set_title(
        "Anchor-matched CTS signal across conservation classes\n"
        "expected: sfCRE > sdCRE > soCRE > ssCRE; sf_proj ≈ sfCRE",
        fontsize=11,
        pad=10,
    )
    ax.set_xlabel("conservation class", fontsize=10)
    ax.set_ylabel("")

    plt.tight_layout()
    out = FIG_DIR / "cts_anchor_matched_152.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"  -> {out}")


def build_all_classes_all_tissues_matrix(df, trait_meta, value_col="signed_log10p"):
    """
    Build the combined matrix:
      rows = traits
      columns = tissue by class
      values = value_col

    Valid value_col choices:
      - signed_log10p
      - signed_z
    """
    sub = df.copy()
    sub["tissue_class"] = sub["tissue"] + "|" + sub["class"]

    pivot = sub.pivot_table(
        index="trait_id",
        columns="tissue_class",
        values=value_col,
        aggfunc="first",
    )

    # Order rows by anchor_tissue and Description.
    trait_order = trait_meta["trait_id"].tolist()
    pivot = pivot.reindex(index=trait_order)

    # Order columns by tissue and then class.
    ordered_cols = []

    for t in TISSUE_ORDER:
        for c in CLASSES:
            tc = f"{t}|{c}"
            if tc in pivot.columns:
                ordered_cols.append(tc)

    pivot = pivot.reindex(columns=ordered_cols)

    # Build readable row labels.
    desc_map = trait_meta.set_index("trait_id")["Description"].to_dict()
    tissue_map = trait_meta.set_index("trait_id")["anchor_tissue"].to_dict()

    pivot.index = [
        safe_short_label(t, desc_map, tissue_map, max_len=42)
        for t in pivot.index
    ]

    return pivot


def plot_all_classes_all_tissues(df, trait_meta):
    """
    Draw the ordered combined heatmap:
    rows = traits
    columns = tissue by class
    values = signed_log10p

    Example column order:
    Liver|sfCRE, Liver|sdCRE, Liver|soCRE, Liver|ssCRE, Liver|sf_proj,
    Spleen|sfCRE, ...
    """
    pivot = build_all_classes_all_tissues_matrix(
        df,
        trait_meta,
        value_col="signed_log10p",
    )

    pivot.to_csv(
        SUMMARY_DIR / "A2_B2_cts_152_alltissue_allclass_signedlog10p.csv"
    )

    n_row, n_col = pivot.shape

    fig, ax = plt.subplots(
        figsize=(max(14, n_col * 0.42), max(8, n_row * 0.28))
    )

    sns.heatmap(
        pivot,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        cbar_kws={"label": "signed -log10(P)", "shrink": 0.5},
        linewidths=0.2,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
    )

    # Add horizontal separators between trait anchor groups.
    tm = trait_meta.copy().reset_index(drop=True)

    cps = tm.index[
        tm["anchor_tissue"] != tm["anchor_tissue"].shift()
    ].tolist()

    for cp in cps[1:]:
        ax.axhline(cp, color="gray", lw=0.8)

    # Add vertical separators between tissues.
    col_count = 0

    for t in TISSUE_ORDER[:-1]:
        k = sum(1 for c in CLASSES if f"{t}|{c}" in pivot.columns)
        col_count += k
        ax.axvline(col_count, color="black", lw=1.0)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)

    ax.set_title(
        "CTS global overview: all traits × all tissue-class combinations",
        fontsize=12,
    )
    ax.set_xlabel("tissue | conservation class")
    ax.set_ylabel("")

    plt.tight_layout()
    out = FIG_DIR / "cts_all_tissues_all_classes_152.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"  -> {out}")


def plot_all_classes_all_tissues_clustered(
    df,
    trait_meta,
    value_col="signed_log10p",
    row_scale=False,
    clip_value=3,
):
    """
    Draw the clustered combined heatmap:
    rows = traits
    columns = tissue by class

    Three output variants are generated:
      1. signed_log10p:
         Shows the significance-evidence pattern of LDSC enrichment.
      2. signed_z:
         Uses coefficient divided by SE and is more suitable for clustering.
      3. signed_z + row_scale:
         Reduces differences in GWAS power and emphasizes within-trait preferences.

    Notes:
      - clustermap cannot process NaN values, so missing values are filled with zero.
      - Zero approximates the absence of a directional signal for signed_log10p and signed_z.
      - Clustered heatmaps are exploratory and do not imply that raw coefficients are directly comparable.
    """
    if value_col not in {"signed_log10p", "signed_z"}:
        raise ValueError("value_col must be 'signed_log10p' or 'signed_z'")

    pivot = build_all_classes_all_tissues_matrix(
        df,
        trait_meta,
        value_col=value_col,
    )

    if value_col == "signed_log10p":
        base_suffix = "signed_log10p"
        label = "signed -log10(P)"
    elif value_col == "signed_z":
        base_suffix = "signed_z"
        label = "signed Z-score"
    else:
        base_suffix = value_col
        label = value_col

    if row_scale:
        suffix = f"{base_suffix}_rowscaled"
    else:
        suffix = base_suffix

    # Save the unclustered matrix before missing-value filling.
    pivot.to_csv(
        SUMMARY_DIR / f"A2_B2_cts_152_alltissue_allclass_{suffix}_for_cluster.csv"
    )

    # clustermap cannot contain NaN values.
    mat = pivot.dropna(axis=0, how="all").dropna(axis=1, how="all").fillna(0)

    if mat.shape[0] < 2 or mat.shape[1] < 2:
        print(f"  [WARN] The {suffix} matrix is too small; skipping clustering")
        return

    if row_scale:
        # Apply within-trait z-scores to emphasize relative preferences.
        row_mean = mat.mean(axis=1)
        row_std = mat.std(axis=1).replace(0, np.nan)

        mat_for_cluster = mat.sub(row_mean, axis=0).div(row_std, axis=0)
        mat_for_cluster = mat_for_cluster.replace([np.inf, -np.inf], np.nan).fillna(0)

        heatmap_label = f"row-scaled {label}"
        title_extra = "row-scaled within each trait"
        vmax = clip_value
        vmin = -clip_value
    else:
        mat_for_cluster = mat.copy()
        heatmap_label = f"{label}, clipped to [-{clip_value}, {clip_value}]"
        title_extra = "unscaled"
        vmax = clip_value
        vmin = -clip_value

    # Limit the dynamic range so extreme P or Z values do not dominate distances.
    mat_for_cluster = mat_for_cluster.clip(lower=-clip_value, upper=clip_value)

    # Save the matrix actually used for clustering.
    mat_for_cluster.to_csv(
        SUMMARY_DIR / f"A2_B2_cts_152_alltissue_allclass_{suffix}_cluster_input.csv"
    )

    g = sns.clustermap(
        mat_for_cluster,
        cmap="RdBu_r",
        center=0,
        vmin=vmin,
        vmax=vmax,
        method="average",
        metric="euclidean",
        figsize=(
            max(14, mat_for_cluster.shape[1] * 0.42),
            max(8, mat_for_cluster.shape[0] * 0.26),
        ),
        linewidths=0.1,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": heatmap_label},
        dendrogram_ratio=(0.12, 0.12),
        colors_ratio=0.02,
    )

    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(),
        rotation=90,
        fontsize=7,
    )
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(),
        rotation=0,
        fontsize=5,
    )

    g.ax_heatmap.set_xlabel("clustered tissue | conservation class")
    g.ax_heatmap.set_ylabel("clustered traits")

    g.fig.suptitle(
        "Clustered CTS global overview: all traits × all tissue-class combinations\n"
        f"{title_extra}; clustering reflects LDSC evidence pattern, not raw effect-size comparability",
        fontsize=12,
        y=1.03,
    )

    out = FIG_DIR / f"cts_all_tissues_all_classes_152_clustered_{suffix}.png"
    plt.savefig(out, dpi=220, bbox_inches="tight")
    plt.close()

    print(f"  -> {out}")

    # Save clustered row and column orders for interpretation.
    row_order = g.dendrogram_row.reordered_ind
    col_order = g.dendrogram_col.reordered_ind

    clustered_rows = mat_for_cluster.index[row_order]
    clustered_cols = mat_for_cluster.columns[col_order]

    pd.Series(clustered_rows, name="clustered_trait").to_csv(
        SUMMARY_DIR / f"A2_B2_cts_152_clustered_{suffix}_trait_order.csv",
        index=False,
    )

    pd.Series(clustered_cols, name="clustered_tissue_class").to_csv(
        SUMMARY_DIR / f"A2_B2_cts_152_clustered_{suffix}_tissue_class_order.csv",
        index=False,
    )


def main():
    print("[1/5] Scanning .cell_type_results.txt files...")
    df = collect_all()

    print(f"      Total rows: {len(df)}")
    print("      Class distribution:")
    print(df["class"].value_counts().to_string())

    long_csv = SUMMARY_DIR / "A2_B2_cts_152_long.csv"
    df.to_csv(long_csv, index=False)
    print(f"  -> {long_csv}")

    print("[2/5] Drawing five per-class trait-by-tissue heatmaps...")
    trait_meta = order_traits(df)

    for c in CLASSES:
        plot_one_class_heatmap(df, c, trait_meta)

    print("[3/5] Drawing the anchor-matched summary...")
    plot_anchor_matched_summary(df, trait_meta)

    print("[4/5] Drawing the ordered all-traits by class-tissue heatmap...")
    plot_all_classes_all_tissues(df, trait_meta)

    print("[5/5] Drawing clustered all-traits by class-tissue heatmaps...")

    # Version 1: signed -log10(P).
    # Interpretation: significance evidence across trait-class-tissue combinations.
    plot_all_classes_all_tissues_clustered(
        df,
        trait_meta,
        value_col="signed_log10p",
        row_scale=False,
        clip_value=3,
    )

    # Version 2: signed Z-score.
    # Interpretation: coefficient divided by SE, which is preferable to raw coefficients for clustering.
    plot_all_classes_all_tissues_clustered(
        df,
        trait_meta,
        value_col="signed_z",
        row_scale=False,
        clip_value=5,
    )

    # Version 3: row-scaled signed Z-score.
    # Interpretation: reduces cross-trait GWAS power differences and highlights
    # relative tissue-class preferences within each trait.
    plot_all_classes_all_tissues_clustered(
        df,
        trait_meta,
        value_col="signed_z",
        row_scale=True,
        clip_value=3,
    )

    print("\nCompleted. Key outputs:")
    print(f"  CSV:   {SUMMARY_DIR}")
    print(f"  Figs:  {FIG_DIR}")
    print("  -> ordered heatmap:")
    print("     figures/cts_all_tissues_all_classes_152.png")
    print("  -> clustered heatmaps:")
    print("     figures/cts_all_tissues_all_classes_152_clustered_signed_log10p.png")
    print("     figures/cts_all_tissues_all_classes_152_clustered_signed_z.png")
    print("     figures/cts_all_tissues_all_classes_152_clustered_signed_z_rowscaled.png")
    print("  -> anchor-matched summary:")
    print("     figures/cts_anchor_matched_152.png")


if __name__ == "__main__":
    main()
