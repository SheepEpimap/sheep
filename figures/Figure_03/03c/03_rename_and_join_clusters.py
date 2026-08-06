#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
"""
 ID JASPAR file
 cluster_1, cluster_2, ..., cluster_N
 all_jaspar_with_cluster.tsv
"""

import pandas as pd
import os

def rename_clusters_and_merge():
    """
     ID file
    """
    cluster_file_path = "/data/home/sczd644/run/zsw_chrombpnet/cluster/cross_tissue_cluster_analysis/cluster_results/cluster_key.txt"
    annotation_file_path = "/data/home/sczd644/run/zsw_chrombpnet/all_jaspar_with_cluster.tsv"
    output_file_path = "/data/home/sczd644/run/zsw_chrombpnet/cluster/all_jaspar_with_renamed_clusters.tsv"

    if not os.path.exists(cluster_file_path):
        print(f" : Not found file {cluster_file_path}")
        return

    if not os.path.exists(annotation_file_path):
        print(f" : Not found file {annotation_file_path}")
        return

    print(" ...")

    print("1. read file ID...")
    cluster_data = []
    with open(cluster_file_path, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines, 1):
        line = line.strip()
        if not line or '\t' not in line:
            continue

        parts = line.split('\t')
        if len(parts) < 2:
            continue

        old_cluster_id = parts[0].strip()
        patterns_str = parts[1].strip()

        new_cluster_id = f"cluster_{i}"

        patterns = [p.strip() for p in patterns_str.split(',') if p.strip()]

        for pattern in patterns:
            cluster_data.append({
                'pattern_identifier': pattern,
                'old_cluster_id': old_cluster_id,
                'new_cluster_id': new_cluster_id
            })

    df_cluster_map = pd.DataFrame(cluster_data)
    print(f"   : {len(df_cluster_map)}  , {df_cluster_map['new_cluster_id'].nunique()}  ")

    print("2. readJASPAR file...")
    df_annotation = pd.read_csv(annotation_file_path, sep='\t', usecols=range(7))
    print(f"  read : {len(df_annotation)}  , {df_annotation.columns.tolist()}  ")

    df_annotation['pattern_identifier'] = df_annotation['tissue'] + '__' + df_annotation['pattern']

    print("3.  ...")
    df_merged = pd.merge(df_annotation, df_cluster_map[['pattern_identifier', 'new_cluster_id']],
                        on='pattern_identifier', how='left')

    df_merged = df_merged.drop('pattern_identifier', axis=1)

    matched_count = df_merged['new_cluster_id'].notna().sum()
    unmatched_count = df_merged['new_cluster_id'].isna().sum()

    print(f"4.  results: {matched_count}  , {unmatched_count}  ")

    df_merged.to_csv(output_file_path, sep='\t', index=False)
    print(f"5.  results save: {output_file_path}")

    print("\n" + "=" * 50)
    print(" !")
    print("=" * 50)

    cluster_sizes = df_merged['new_cluster_id'].value_counts()
    print(f"\n statistics ( 10 ):")
    for cluster_id, size in cluster_sizes.head(10).items():
        print(f"  {cluster_id}: {size}  ")

    tissue_cluster_counts = df_merged.groupby('tissue')['new_cluster_id'].nunique().sort_values(ascending=False)
    print(f"\n  ( 10 ):")
    for tissue, count in tissue_cluster_counts.head(10).items():
        print(f"  {tissue}: {count}  ")

    rename_map_file = "/data/home/sczd644/run/zsw_chrombpnet/cluster/cluster_rename_mapping.tsv"
    df_cluster_map[['old_cluster_id', 'new_cluster_id']].drop_duplicates().to_csv(rename_map_file, sep='\t', index=False)
    print(f"\n save: {rename_map_file}")

    return df_merged

def analyze_clustering_quality():
    """
     , 700+ motif 27
    """
    print("\n" + "=" * 50)
    print(" ")
    print("=" * 50)

    merged_file = "/data/home/sczd644/run/zsw_chrombpnet/cluster/all_jaspar_with_renamed_clusters.tsv"

    if not os.path.exists(merged_file):
        print(" : Not found file")
        return

    df = pd.read_csv(merged_file, sep='\t')
    df_with_cluster = df[df['new_cluster_id'].notna()]

    total_motifs = len(df)
    clustered_motifs = len(df_with_cluster)
    num_clusters = df_with_cluster['new_cluster_id'].nunique()

    print(f" motif : {total_motifs}")
    print(f" motif : {clustered_motifs}")
    print(f" : {num_clusters}")
    print(f" : {clustered_motifs/total_motifs*100:.1f}%")

    cluster_sizes = df_with_cluster['new_cluster_id'].value_counts()
    print(f"\n :")
    print(f"   : {cluster_sizes.max()}  motif")
    print(f"   : {cluster_sizes.min()}  motif")
    print(f"   : {cluster_sizes.mean():.1f}  motif")
    print(f"   : {cluster_sizes.median()}  motif")

    large_clusters = cluster_sizes[cluster_sizes > 50]
    if len(large_clusters) > 0:
        print(f"\n :   {len(large_clusters)}   (>50 motif):")
        for cluster_id, size in large_clusters.items():
            print(f"  {cluster_id}: {size}  motif")

    print(f"\n :")
    print(f"  1.  motif , motif ")
    print(f"  2.  ")
    print(f"  3.  ")
    print(f"  4.  ")

if __name__ == '__main__':
    df_result = rename_clusters_and_merge()

    analyze_clustering_quality()
