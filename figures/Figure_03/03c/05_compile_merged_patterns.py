#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# Purpose: Compile merged motifs
# This script combines several h5 modisco objects into one.
# Specifically, we use the h5 modisco files created from merging modisco
# motifs within gimme cluster patterns of motifs.

import sys
import h5py as h5
import numpy as np
import pandas as pd
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Compile merged modisco patterns from cluster directories.")
    parser.add_argument("--merged-dir", type=str, required=True,
                       help="Directory containing merged results (with cluster_1, cluster_2, etc. subdirectories)")
    parser.add_argument("--out-dir", type=str, required=True,
                       help="Output directory for compiled results")

    args = parser.parse_args()

    # Validate inputs
    assert os.path.exists(args.merged_dir), f"Merged directory {args.merged_dir} not found."
    os.makedirs(args.out_dir, exist_ok=True)

    return args

def main():
    args = parse_args()

    # Set up output paths
    compiled_modisco_h5_path = os.path.join(args.out_dir, "modisco_compiled.h5")
    compiled_modisco_tsv_path = os.path.join(args.out_dir, "modisco_compiled.tsv")

    print("Starting compilation of merged motifs...")
    print(f"Input directory: {args.merged_dir}")
    print(f"Output directory: {args.out_dir}")

    # COMPILE PATTERNS ----------------------------------------------
    # Create a dict from cluster to modisco merged file
    input_files = {}

    # Find all cluster directories
    for item in os.listdir(args.merged_dir):
        cluster_dir = os.path.join(args.merged_dir, item)
        if os.path.isdir(cluster_dir) and item.startswith("cluster_"):
            h5_file = os.path.join(cluster_dir, "merged_patterns.h5")
            if os.path.exists(h5_file):
                input_files[item] = h5_file

    print(f"Found {len(input_files)} cluster directories with merged_patterns.h5")

    # Open the compiled modisco h5
    with h5.File(compiled_modisco_h5_path, "w") as compiled_modisco:

        # Set up groups in h5 file
        h5_pattern_groups = {
            "pos_patterns": compiled_modisco.create_group('pos_patterns'),
            "neg_patterns": compiled_modisco.create_group('neg_patterns')
        }

        # Loop over the input files containing merged patterns derived from collapsing/merging similar patterns
        # within each gimme cluster
        for cluster_name, h5_file_path in input_files.items():

            print(f"Processing: {cluster_name} | input file: {h5_file_path}")

            try:
                with h5.File(h5_file_path) as modisco_obj:

                    # Figure out what class of patterns the object contains
                    if "pos_patterns" in modisco_obj:
                        pattern_class = "pos_patterns"
                    elif "neg_patterns" in modisco_obj:
                        pattern_class = "neg_patterns"
                    else:
                        print(f"  WARNING: No pattern class found in {h5_file_path}, skipping...")
                        continue

                    current_group = h5_pattern_groups[pattern_class]

                    # Iterate over patterns in the obj and add them to the compiled object
                    for pattern in modisco_obj[pattern_class].keys():

                        # Create the pattern group
                        pattern_group = current_group.create_group(pattern)

                        # Copy datasets
                        pattern_group.create_dataset(
                            "contrib_scores",
                            data=modisco_obj[pattern_class][pattern]["contrib_scores"]
                        )
                        pattern_group.create_dataset(
                            "sequence",
                            data=modisco_obj[pattern_class][pattern]["sequence"]
                        )
                        pattern_group.create_dataset(
                            "hypothetical_contribs",
                            data=modisco_obj[pattern_class][pattern]["hypothetical_contribs"]
                        )

                        # Add in placeholder value for number of seqlets
                        seqlets_group = pattern_group.create_group("seqlets")
                        seqlets_group.create_dataset("n_seqlets", data=np.array([1]))

            except Exception as e:
                print(f"  ERROR processing {cluster_name}: {str(e)}")
                continue

    print("Combining TSV reports...")

    tsv_files = []
    for item in os.listdir(args.merged_dir):
        cluster_dir = os.path.join(args.merged_dir, item)
        if os.path.isdir(cluster_dir) and item.startswith("cluster_"):
            tsv_file = os.path.join(cluster_dir, "merge_report_tf.tsv")
            if os.path.exists(tsv_file):
                tsv_files.append(tsv_file)
            else:
                fallback_tsv = os.path.join(cluster_dir, "merge_report.tsv")
                if os.path.exists(fallback_tsv):
                    print(f"  WARNING: Using fallback file {fallback_tsv} (no TF information)")
                    tsv_files.append(fallback_tsv)

    if tsv_files:
        dfs = []
        for file in tsv_files:
            try:
                df = pd.read_csv(file, sep="\t")
                if 'TF' not in df.columns:
                    df['TF'] = 'NA'
                dfs.append(df)
                print(f"  Loaded: {os.path.basename(os.path.dirname(file))}/{os.path.basename(file)} ({len(df)} rows)")
            except Exception as e:
                print(f"  ERROR reading {file}: {str(e)}")
                continue

        if dfs:
            merged_report_long = pd.concat(dfs, ignore_index=True)

            expected_columns = ['merged_pattern', 'input_tissue', 'input_pattern', 'input_pattern_full',
                              'n_seqlets', 'total_seqlets_in_merged', 'pattern_class', 'TF']

            existing_columns = [col for col in expected_columns if col in merged_report_long.columns]
            merged_report_long = merged_report_long[existing_columns]

            merged_report_long.to_csv(compiled_modisco_tsv_path, sep="\t", index=False)
            print(f"Combined {len(dfs)} TSV reports into {compiled_modisco_tsv_path}")
            print(f"Final compiled TSV has {len(merged_report_long)} rows and {len(merged_report_long.columns)} columns")
            print(f"Columns: {', '.join(merged_report_long.columns)}")
        else:
            print("No valid TSV data found to combine")
    else:
        print("No TSV reports found to combine")

    print("Compilation completed successfully!")
    print(f"Compiled H5 file: {compiled_modisco_h5_path}")
    print(f"Compiled TSV file: {compiled_modisco_tsv_path}")

if __name__ == "__main__":
    main()
