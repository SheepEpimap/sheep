#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import os
import pandas as pd


# ============================================================
# ============================================================

ENH_DIR = "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/enhancers_by_tissue_count"

TARGET_FILE = "/vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_output_E5_confident.tsv"

OUTDIR = "/vol2/zhangshiwen/sheep_cor/enhancer_tissue_count_distribution/E5_5groups_simple_target_gene"

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# ============================================================

GROUPS = {
    "G1_1_tissue": list(range(1, 2)),
    "G2_2_5_tissues": list(range(2, 6)),
    "G3_6_10_tissues": list(range(6, 11)),
    "G4_11_20_tissues": list(range(11, 21)),
    "G5_21_43_tissues": list(range(21, 44)),
}

GROUP_LABEL = {
    "G1_1_tissue": "1",
    "G2_2_5_tissues": "2-5",
    "G3_6_10_tissues": "6-10",
    "G4_11_20_tissues": "11-20",
    "G5_21_43_tissues": "21-43",
}


# ============================================================
# 3. read enhancer-target gene file
# ============================================================

target_cols = [
    "enhancer_id",
    "spearman_rho",
    "pvalue",
    "gene",
    "distance",
    "qvalue"
]

target = pd.read_csv(
    TARGET_FILE,
    sep="\t",
    header=None,
    names=target_cols,
    dtype=str
)

for col in target.columns:
    target[col] = target[col].astype(str).str.strip()

for col in ["spearman_rho", "pvalue", "distance", "qvalue"]:
    target[col] = pd.to_numeric(target[col], errors="coerce")

target = target.drop_duplicates(
    subset=["enhancer_id", "gene", "spearman_rho", "pvalue", "distance", "qvalue"],
    keep="first"
).copy()


# ============================================================
# ============================================================

def read_one_enhancer_file(n_tissues, group_name):
    """
    read enhancers_present_in_xx_tissue(s).tsv file.
      chr/start/end/n_tissues/enhancer_id/group  .
     43 0/1 .
    """

    if n_tissues == 1:
        filename = f"enhancers_present_in_{n_tissues:02d}_tissue.tsv"
    else:
        filename = f"enhancers_present_in_{n_tissues:02d}_tissues.tsv"

    filepath = os.path.join(ENH_DIR, filename)

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Not foundfile: {filepath}")

    df = pd.read_csv(filepath, sep="\t", header=0, dtype=str)
    df.columns = [x.strip() for x in df.columns]

    for col in ["chr", "start", "end"]:
        if col not in df.columns:
            raise ValueError(f"{filepath} Missing required columns: {col}")

    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")

    if "n_tissues" in df.columns:
        df["n_tissues"] = pd.to_numeric(df["n_tissues"], errors="coerce").astype("Int64")
    else:
        df["n_tissues"] = int(n_tissues)

    df["tissue_count_group"] = group_name
    df["tissue_count_range"] = GROUP_LABEL[group_name]

    df["enhancer_id"] = (
        df["chr"].astype(str)
        + ":"
        + df["start"].astype(str)
        + "-"
        + df["end"].astype(str)
    )

    df = df[
        [
            "tissue_count_group",
            "tissue_count_range",
            "n_tissues",
            "enhancer_id",
            "chr",
            "start",
            "end"
        ]
    ].copy()

    return df


# ============================================================
# ============================================================

all_enhancer_list = []
summary_input_rows = []

for group_name, n_list in GROUPS.items():

    group_df_list = []

    for n in n_list:
        df_n = read_one_enhancer_file(n, group_name)
        group_df_list.append(df_n)

        summary_input_rows.append({
            "tissue_count_group": group_name,
            "tissue_count_range": GROUP_LABEL[group_name],
            "n_tissues": n,
            "enhancer_rows": df_n.shape[0],
            "unique_enhancer_id": df_n["enhancer_id"].nunique()
        })

    group_enh = pd.concat(group_df_list, axis=0, ignore_index=True)

    group_file = os.path.join(
        OUTDIR,
        f"{group_name}.simple_enhancers.tsv"
    )

    group_enh.to_csv(group_file, sep="\t", index=False, na_rep="NA")

    all_enhancer_list.append(group_enh)


all_enhancers = pd.concat(all_enhancer_list, axis=0, ignore_index=True)

all_enhancers_file = os.path.join(
    OUTDIR,
    "all_5groups.simple_enhancers.tsv"
)

all_enhancers.to_csv(all_enhancers_file, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

left_join = all_enhancers.merge(
    target,
    on="enhancer_id",
    how="left"
)

left_join_simple = left_join[
    [
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_id",
        "chr",
        "start",
        "end",
        "gene",
        "spearman_rho",
        "pvalue",
        "distance",
        "qvalue"
    ]
].copy()

left_join_file = os.path.join(
    OUTDIR,
    "all_5groups.simple_left_join_all_enhancers_with_target_genes.tsv"
)

left_join_simple.to_csv(left_join_file, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

linked_pairs = all_enhancers.merge(
    target,
    on="enhancer_id",
    how="inner"
)

linked_pairs_simple = linked_pairs[
    [
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_id",
        "chr",
        "start",
        "end",
        "gene",
        "spearman_rho",
        "pvalue",
        "distance",
        "qvalue"
    ]
].copy()

linked_pairs_simple = linked_pairs_simple.drop_duplicates(
    subset=[
        "tissue_count_group",
        "n_tissues",
        "enhancer_id",
        "gene",
        "spearman_rho",
        "pvalue",
        "distance",
        "qvalue"
    ],
    keep="first"
).copy()

linked_pairs_file = os.path.join(
    OUTDIR,
    "all_5groups.simple_linked_enhancer_gene_pairs.tsv"
)

linked_pairs_simple.to_csv(linked_pairs_file, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

PER_GROUP_DIR = os.path.join(OUTDIR, "per_group_simple_linked_pairs")
os.makedirs(PER_GROUP_DIR, exist_ok=True)

for group_name in GROUPS.keys():
    sub = linked_pairs_simple[
        linked_pairs_simple["tissue_count_group"] == group_name
    ].copy()

    outfile = os.path.join(
        PER_GROUP_DIR,
        f"{group_name}.simple_linked_enhancer_gene_pairs.tsv"
    )

    sub.to_csv(outfile, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

if linked_pairs_simple.shape[0] > 0:
    linked_collapsed = (
        linked_pairs_simple
        .dropna(subset=["gene"])
        .groupby(
            [
                "tissue_count_group",
                "tissue_count_range",
                "n_tissues",
                "enhancer_id",
                "chr",
                "start",
                "end"
            ],
            dropna=False
        )
        .agg(
            target_gene_count=("gene", "nunique"),
            target_gene_list=("gene", lambda x: ",".join(sorted(set(x.dropna().astype(str))))),
            enhancer_gene_pair_count=("gene", "count")
        )
        .reset_index()
    )
else:
    linked_collapsed = pd.DataFrame(
        columns=[
            "tissue_count_group",
            "tissue_count_range",
            "n_tissues",
            "enhancer_id",
            "chr",
            "start",
            "end",
            "target_gene_count",
            "target_gene_list",
            "enhancer_gene_pair_count"
        ]
    )

all_enhancer_unique = all_enhancers.drop_duplicates(
    subset=[
        "tissue_count_group",
        "n_tissues",
        "enhancer_id"
    ],
    keep="first"
).copy()

enhancer_level = all_enhancer_unique.merge(
    linked_collapsed,
    on=[
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_id",
        "chr",
        "start",
        "end"
    ],
    how="left"
)

enhancer_level["target_gene_count"] = enhancer_level["target_gene_count"].fillna(0).astype(int)
enhancer_level["enhancer_gene_pair_count"] = enhancer_level["enhancer_gene_pair_count"].fillna(0).astype(int)
enhancer_level["target_gene_list"] = enhancer_level["target_gene_list"].fillna("NA")

enhancer_level = enhancer_level[
    [
        "tissue_count_group",
        "tissue_count_range",
        "n_tissues",
        "enhancer_id",
        "chr",
        "start",
        "end",
        "target_gene_count",
        "target_gene_list",
        "enhancer_gene_pair_count"
    ]
].copy()

enhancer_level_file = os.path.join(
    OUTDIR,
    "all_5groups.simple_enhancer_level_target_gene_count.tsv"
)

enhancer_level.to_csv(enhancer_level_file, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

PER_GROUP_ENH_LEVEL_DIR = os.path.join(OUTDIR, "per_group_simple_enhancer_level")
os.makedirs(PER_GROUP_ENH_LEVEL_DIR, exist_ok=True)

for group_name in GROUPS.keys():
    sub = enhancer_level[
        enhancer_level["tissue_count_group"] == group_name
    ].copy()

    outfile = os.path.join(
        PER_GROUP_ENH_LEVEL_DIR,
        f"{group_name}.simple_enhancer_level_target_gene_count.tsv"
    )

    sub.to_csv(outfile, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

summary_rows = []

for group_name in GROUPS.keys():

    enh_g = all_enhancers[
        all_enhancers["tissue_count_group"] == group_name
    ].copy()

    enh_level_g = enhancer_level[
        enhancer_level["tissue_count_group"] == group_name
    ].copy()

    pair_g = linked_pairs_simple[
        linked_pairs_simple["tissue_count_group"] == group_name
    ].copy()

    total_enhancer_rows = enh_g.shape[0]
    unique_enhancers = enh_level_g["enhancer_id"].nunique()
    enhancers_with_target = enh_level_g[enh_level_g["target_gene_count"] > 0]["enhancer_id"].nunique()
    enhancer_gene_pairs = pair_g.shape[0]
    unique_target_genes = pair_g["gene"].nunique() if pair_g.shape[0] > 0 else 0

    summary_rows.append({
        "tissue_count_group": group_name,
        "tissue_count_range": GROUP_LABEL[group_name],
        "total_enhancer_rows": total_enhancer_rows,
        "unique_enhancers": unique_enhancers,
        "enhancers_with_target": enhancers_with_target,
        "percent_enhancers_with_target": (
            enhancers_with_target / unique_enhancers * 100
            if unique_enhancers > 0 else 0
        ),
        "enhancer_gene_pairs": enhancer_gene_pairs,
        "unique_target_genes": unique_target_genes,
        "mean_target_genes_per_enhancer": (
            enh_level_g["target_gene_count"].mean()
            if enh_level_g.shape[0] > 0 else 0
        ),
        "median_target_genes_per_enhancer": (
            enh_level_g["target_gene_count"].median()
            if enh_level_g.shape[0] > 0 else 0
        ),
        "mean_target_genes_per_linked_enhancer": (
            enh_level_g.loc[enh_level_g["target_gene_count"] > 0, "target_gene_count"].mean()
            if enhancers_with_target > 0 else 0
        )
    })

summary = pd.DataFrame(summary_rows)

summary_file = os.path.join(
    OUTDIR,
    "summary.simple_5groups_target_gene.tsv"
)

summary.to_csv(summary_file, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

input_count = pd.DataFrame(summary_input_rows)

input_count_file = os.path.join(
    OUTDIR,
    "check.input_43_files.simple_count.tsv"
)

input_count.to_csv(input_count_file, sep="\t", index=False, na_rep="NA")


# ============================================================
# ============================================================

print(" .")
print(f"Output directory: {OUTDIR}")
print("")
print(" resultsfile:")
print(f"1.   enhancer  : {all_enhancers_file}")
print(f"2.   enhancer   left join  : {left_join_file}")
print(f"3.   enhancer-gene pair  : {linked_pairs_file}")
print(f"4.   enhancer  : {enhancer_level_file}")
print(f"5.  statistics summary: {summary_file}")
print("")
print(" results:")
print(f"enhancer-gene pair: {PER_GROUP_DIR}")
print(f"enhancer-level target count: {PER_GROUP_ENH_LEVEL_DIR}")
