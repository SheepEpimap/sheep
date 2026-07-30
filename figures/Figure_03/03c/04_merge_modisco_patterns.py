#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
import h5py as h5
import numpy as np
import pandas as pd
import os
import argparse
from datetime import datetime
from modiscolite.aggregator import SimilarPatternsCollapser
from modiscolite.core import TrackSet, Seqlet, SeqletSet
from modiscolite.report import _plot_weights, path_to_image_html

def parse_args():
    parser = argparse.ArgumentParser(description="Merge modisco patterns across tissues by cluster.")
    parser.add_argument("--cluster-tsv", type=str, required=True,
                       help="all_jaspar_with_cluster.tsv file containing cluster information")
    parser.add_argument("--tissues-file", type=str, required=True,
                       help="tissue.txt file containing list of tissue names")
    parser.add_argument("--modisco-dir", type=str, required=True,
                       help="Directory containing modisco h5 files for each tissue")
    parser.add_argument("--contribs-dir", type=str, required=True,
                       help="Directory containing contribution scores h5 files for each tissue")
    parser.add_argument("--out-dir", type=str, required=True,
                       help="Output directory for merged patterns and reports")
    parser.add_argument("--batch", type=int, required=False, default=None,
                       help="Batch number for parallel processing (optional)")
    parser.add_argument("--debug", action="store_true", default=False,
                       help="Debug mode - process only first cluster")

    args = parser.parse_args()

    # Validate inputs
    assert os.path.exists(args.cluster_tsv), f"Cluster TSV file {args.cluster_tsv} not found."
    assert os.path.exists(args.tissues_file), f"Tissues file {args.tissues_file} not found."
    assert os.path.exists(args.modisco_dir), f"Modisco directory {args.modisco_dir} not found."
    assert os.path.exists(args.contribs_dir), f"Contribs directory {args.contribs_dir} not found."
    os.makedirs(args.out_dir, exist_ok=True)

    return args

def load_tissues_and_clusters(args):
    """Load tissue list and cluster information"""

    # Load tissue list
    with open(args.tissues_file, 'r') as f:
        tissues = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(tissues)} tissues: {tissues[:5]}...")

    # Load cluster information
    cluster_df = pd.read_csv(args.cluster_tsv, sep='\t')
    print(f"Loaded cluster data with {len(cluster_df)} patterns")

    # Group by cluster
    clusters_to_merge = []
    for cluster_id, group in cluster_df.groupby('new_cluster_id'):
        components = []
        for _, row in group.iterrows():
            tissue = row['tissue']
            pattern_full = row['pattern_full']
            # Extract pattern name from pattern_full (e.g., "abomasum_pattern_0" -> "pattern_0")
            pattern_name = pattern_full.replace(f"{tissue}_", "")

            components.append({
                "tissue": tissue,
                "pattern_name": pattern_name,
                "pattern_full": pattern_full,
                "seqlets_count": row['seqlets'],
                "motif": row['motif'],
                "TF": row['TF']
            })

        clusters_to_merge.append({
            "cluster_id": cluster_id,
            "components": components
        })

    # Filter for batch if specified
    if args.batch is not None:
        # Simple batching: take every nth cluster based on batch number
        clusters_to_merge = [clusters_to_merge[i] for i in range(len(clusters_to_merge))
                           if i % args.num_batches == args.batch]
        print(f"Processing batch {args.batch}: {len(clusters_to_merge)} clusters")

    # Debug mode: process only first cluster
    if args.debug:
        clusters_to_merge = clusters_to_merge[:1]
        print("Debug mode: processing first cluster only")

    print(f"Total clusters to process: {len(clusters_to_merge)}")
    return clusters_to_merge

def get_file_paths(tissue, modisco_dir, contribs_dir):
    """Get file paths for a given tissue"""
    # Modisco file - based on your example: hippocampus.modisco_results.h5
    modisco_file = f"{modisco_dir}/{tissue}.modisco_results.h5"

    # Contribs file - based on your example: cerebral-cortex_chrombpnet_contribs.counts_scores.h5
    contribs_file = f"{contribs_dir}/{tissue}_chrombpnet_contribs.counts_scores.h5"

    return modisco_file, contribs_file

def calculate_window_offsets(center: int, window_size: int) -> tuple:
    return (center - window_size // 2, center + window_size // 2)

def load_pattern_data(cluster, modisco_dir, contribs_dir):
    """Load pattern data for all components in a cluster"""

    # Merge parameters (same as original)
    input_window = 400
    min_overlap = 0.7
    prob_and_pertrack_sim_merge_thresholds = [(0.8,0.8), (0.5, 0.85), (0.2, 70.9)]
    prob_and_pertrack_sim_dealbreaker_thresholds = [(0.4, 0.75), (0.2,0.8), (0.1, 0.85), (0.0,0.9)]
    min_frac = 0.2
    min_num = 30
    flank_to_add = 5
    window_size = 20
    max_seqlets_subsample = 300

    # Storage for combined data
    union_onehot = []
    union_hypscores = []
    union_projscores = []
    union_patterncoords = []
    all_pattern_idxs = []

    exampleidx_offset = 0
    pattern_class = "pos_patterns"  # Based on your data

    print(f"\tLoading data for {len(cluster['components'])} patterns in cluster {cluster['cluster_id']}")

    # Process each tissue-pattern combination
    for component in cluster['components']:
        tissue = component['tissue']
        pattern_name = component['pattern_name']
        pattern_full = component['pattern_full']

        print(f"\t\tProcessing {tissue} - {pattern_name}")

        modisco_file, contribs_file = get_file_paths(tissue, modisco_dir, contribs_dir)

        # Check if files exist
        if not os.path.exists(modisco_file):
            print(f"\t\tWARNING: Modisco file not found: {modisco_file}")
            continue
        if not os.path.exists(contribs_file):
            print(f"\t\tWARNING: Contribs file not found: {contribs_file}")
            continue

        try:
            with h5.File(contribs_file, 'r') as contribs_fh, h5.File(modisco_file, 'r') as modisco_fh:

                # Load contribution scores data
                onehot = contribs_fh["raw"]["seq"][:]
                hypscores = contribs_fh["shap"]["seq"][:]
                projscores = contribs_fh["projected_shap"]["seq"][:]

                print(f"\t\t\tData shapes - onehot: {onehot.shape}, hypscores: {hypscores.shape}, projscores: {projscores.shape}")

                # Load pattern information from modisco file
                pattern_path = f"{pattern_class}/{pattern_name}"
                if pattern_path not in modisco_fh:
                    print(f"\t\tWARNING: Pattern {pattern_path} not found in {modisco_file}")
                    continue

                seqlets_grp = modisco_fh[f"{pattern_path}/seqlets"]
                pattern_exampleidxs = np.array(seqlets_grp['example_idx'])

                if len(pattern_exampleidxs) == 0:
                    print(f"\t\tWARNING: No seqlets found for {pattern_path}")
                    continue

                # Get unique example indices that have seqlets
                surviving_indices = sorted(list(set(pattern_exampleidxs)))
                print(f"\t\t\t{len(surviving_indices)} indices had seqlets out of {len(onehot)}")

                # Extract window around center for surviving sequences
                center = onehot.shape[2] // 2
                start, end = calculate_window_offsets(center, input_window)

                for idx in surviving_indices:
                    union_onehot.append(onehot[idx, :, start:end])
                    union_hypscores.append(hypscores[idx, :, start:end])
                    union_projscores.append(projscores[idx, :, start:end])

                # Create index remapping
                idx_remapping = dict(zip(surviving_indices, np.arange(len(surviving_indices))))

                # Prepare seqlet coordinates with remapped indices
                pattern_remapped_exampleidxs = np.array([
                    (exampleidx_offset + idx_remapping[idx]) for idx in pattern_exampleidxs
                ])
                pattern_start = np.array(seqlets_grp['start'])
                pattern_end = np.array(seqlets_grp['end'])
                pattern_isrevcomp = np.array(seqlets_grp['is_revcomp'])

                union_patterncoords.append((
                    pattern_remapped_exampleidxs,
                    pattern_start,
                    pattern_end,
                    pattern_isrevcomp
                ))

                # Store pattern index information
                all_pattern_idxs.append({
                    "tissue": tissue,
                    "pattern_name": pattern_name,
                    "pattern_full": pattern_full
                })

                exampleidx_offset += len(surviving_indices)

        except Exception as e:
            print(f"\t\tERROR processing {tissue}/{pattern_name}: {str(e)}")
            continue

    return (union_onehot, union_hypscores, union_projscores,
            union_patterncoords, all_pattern_idxs, pattern_class)

def calculate_total_seqlets_per_merged_pattern(pattern_merge_hierarchy, merged_pattern_names):
    """Calculate total seqlets count for each merged pattern"""
    merged_pattern_seqlets = {}

    for i, root_node in enumerate(pattern_merge_hierarchy.root_nodes):
        merged_pattern_name = merged_pattern_names[i]
        total_seqlets = len(root_node.pattern.seqlets)
        merged_pattern_seqlets[merged_pattern_name] = total_seqlets

    return merged_pattern_seqlets

def merge_patterns_for_cluster(cluster, data_tuple, out_dir):
    """Merge patterns for a single cluster and save results"""

    union_onehot, union_hypscores, union_projscores, union_patterncoords, all_pattern_idxs, pattern_class = data_tuple

    if len(union_onehot) == 0:
        print(f"\tNo data to merge for cluster {cluster['cluster_id']}")
        return

    # Create output directory for this cluster
    cluster_out_dir = os.path.join(out_dir, cluster['cluster_id'])
    os.makedirs(cluster_out_dir, exist_ok=True)

    # Reshape data for TrackSet
    union_onehot = np.transpose(np.array(union_onehot), axes=(0, 2, 1))
    union_hypscores = np.transpose(np.array(union_hypscores), axes=(0, 2, 1))
    union_projscores = np.transpose(np.array(union_projscores), axes=(0, 2, 1))

    print(f"\tCombined data shapes - onehot: {union_onehot.shape}, hypscores: {union_hypscores.shape}")

    # Create TrackSet
    track_set = TrackSet(
        one_hot=union_onehot,
        contrib_scores=union_projscores,
        hypothetical_contribs=union_hypscores
    )

    # Calculate background frequency
    bg_freq = np.mean(union_onehot, axis=(0, 1))

    # Create pattern objects
    all_patterns = []
    for (example_idxs, starts, ends, isrevcomps) in union_patterncoords:
        seqlet_coords = [Seqlet(example_idx, start, end, isrevcomp)
                        for (example_idx, start, end, isrevcomp) in zip(
                            example_idxs, starts, ends, isrevcomps)]
        seqlets = track_set.create_seqlets(seqlet_coords)
        pattern = SeqletSet(seqlets)
        all_patterns.append(pattern)
        print(f"\t\tAdded pattern with {len(seqlets)} seqlets")

    # Merge patterns
    print(f"\tMerging {len(all_patterns)} patterns...")
    merged_patterns, pattern_merge_hierarchy = SimilarPatternsCollapser(
        patterns=all_patterns,
        track_set=track_set,
        min_overlap=0.7,
        prob_and_pertrack_sim_merge_thresholds=[(0.8,0.8), (0.5, 0.85), (0.2, 0.9)],
        prob_and_pertrack_sim_dealbreaker_thresholds=[(0.4, 0.75), (0.2,0.8), (0.1, 0.85), (0.0,0.9)],
        min_frac=0.2,
        min_num=30,
        flank_to_add=5,
        window_size=20,
        bg_freq=bg_freq,
        max_seqlets_subsample=300
    )

    print(f"\tFound {len(merged_patterns)} merged patterns after collapsing")

    # Create merged pattern names and calculate total seqlets
    merged_pattern_names = []
    for i in range(len(merged_patterns)):
        merged_pattern_name = f"{cluster['cluster_id']}_merged_{i}"
        merged_pattern_names.append(merged_pattern_name)
        print(f"\t\t{merged_pattern_name}, supported by {len(merged_patterns[i].seqlets)} seqlets")

    # Calculate total seqlets for each merged pattern
    merged_pattern_total_seqlets = calculate_total_seqlets_per_merged_pattern(
        pattern_merge_hierarchy, merged_pattern_names
    )

    # Save merged patterns to H5
    save_merged_patterns_h5(cluster, merged_patterns, pattern_class, cluster_out_dir)

    # Generate reports
    generate_reports(cluster, merged_patterns, pattern_merge_hierarchy,
                    all_patterns, all_pattern_idxs, pattern_class, cluster_out_dir,
                    merged_pattern_total_seqlets)

def save_merged_patterns_h5(cluster, merged_patterns, pattern_class, cluster_out_dir):
    """Save merged patterns to H5 file"""
    h5_path = os.path.join(cluster_out_dir, "merged_patterns.h5")

    with h5.File(h5_path, "w") as f:
        group = f.create_group(pattern_class)
        for i, pattern in enumerate(merged_patterns):
            pattern_name = f"{cluster['cluster_id']}_merged_{i}"
            pattern_group = group.create_group(pattern_name)
            pattern_group.create_dataset("contrib_scores", data=pattern.contrib_scores)
            pattern_group.create_dataset("sequence", data=pattern.sequence)
            pattern_group.create_dataset("hypothetical_contribs", data=pattern.hypothetical_contribs)

    print(f"\tSaved {len(merged_patterns)} merged patterns to {h5_path}")

def find_pattern_matches(node, merged_pattern_name, all_patterns_idxs, merged_logo_path, make_logo=False):
    """Find which original pattern matches the given node"""
    matched_indices = np.nonzero([np.allclose(node.pattern.contrib_scores, x.contrib_scores) for x in all_patterns])[0]
    assert len(matched_indices)==1

    source_celltype = next(iter(all_pattern_idxs[matched_indices[0]].keys()))
    source_pattern = next(iter(all_pattern_idxs[matched_indices[0]].values()))

    if make_logo:
        # produce forward and reverse logos for the input CWM
        logo_path_fwd = f"logos/{source_celltype}_{source_pattern}_fwd.png"
        full_fwd_path = os.path.join(cluster_out_dir, logo_path_fwd)
        _plot_weights(node.pattern.contrib_scores, full_fwd_path)
        input_logo_fwd = f"./{logo_path_fwd}"

        # Reverse logo
        cwm_rev = node.pattern.contrib_scores[::-1, ::-1]
        logo_path_rev = f"logos/{source_celltype}_{source_pattern}_rev.png"
        full_rev_path = os.path.join(cluster_out_dir, logo_path_rev)
        _plot_weights(cwm_rev, full_rev_path)
        input_logo_rev = f"./{logo_path_rev}"
    else:
        input_logo_fwd = None
        input_logo_rev = None

    return [merged_pattern_name, source_celltype, source_pattern, merged_logo_path, input_logo_fwd, input_logo_rev]

def pattern_match_dfs(node, merged_pattern_name, merged_logo_path, make_logo=False,
                     merged_pattern_total_seqlets=None):
    """
    Performs a depth-first search on the PatternMergeHierarchy structure.
    """
    results = []

    # check if the node has child nodes
    if len(node.child_nodes) > 0:
        # recursively visit each child node
        for child in node.child_nodes:
            child_results = pattern_match_dfs(child, merged_pattern_name, merged_logo_path,
                                            make_logo=make_logo,
                                            merged_pattern_total_seqlets=merged_pattern_total_seqlets)
            results.extend(child_results)
    else:
        # we're at a leaf node (empty child_nodes list)
        outs = find_pattern_matches(node, merged_pattern_name, all_pattern_idxs, merged_logo_path, make_logo=make_logo)

        # Get total seqlets for this merged pattern
        total_seqlets_for_merged = merged_pattern_total_seqlets.get(merged_pattern_name, 0) if merged_pattern_total_seqlets else 0

        results.append({
            'merged_pattern': outs[0],
            'merged_logo': outs[3],
            'input_logo_fwd': outs[4],
            'input_logo_rev': outs[5],
            'input_tissue': outs[1],
            'pattern_class': pattern_class,
            'input_pattern': outs[2],
            'n_seqlets': len(node.pattern.seqlets),  # This original pattern's contribution
            'total_seqlets_in_merged': total_seqlets_for_merged  # Total seqlets in the merged pattern
        })

    return results if results else None

def get_leaf_nodes(node):
    """Get all leaf nodes from a pattern merge hierarchy"""
    leaves = []
    if not node.child_nodes:
        leaves.append(node)
    else:
        for child in node.child_nodes:
            leaves.extend(get_leaf_nodes(child))
    return leaves

def find_pattern_index(pattern, all_patterns):
    """Find the index of a pattern in the all_patterns list"""
    for i, p in enumerate(all_patterns):
        if np.allclose(pattern.contrib_scores, p.contrib_scores):
            return i
    return None

def generate_reports(cluster, merged_patterns, pattern_merge_hierarchy,
                    all_patterns, all_pattern_idxs, pattern_class, cluster_out_dir,
                    merged_pattern_total_seqlets=None):
    """Generate TSV and HTML reports"""

    # Create logos directory
    logos_dir = os.path.join(cluster_out_dir, "logos")
    os.makedirs(logos_dir, exist_ok=True)

    # Generate logos for merged patterns
    merged_logo_paths = []
    for i, pattern in enumerate(merged_patterns):
        logo_path = f"logos/{cluster['cluster_id']}_merged_{i}.png"
        full_logo_path = os.path.join(cluster_out_dir, logo_path)
        _plot_weights(pattern.contrib_scores, full_logo_path)
        merged_logo_paths.append(f"./{logo_path}")

    # Generate report data
    report_data = []

    for i, root_node in enumerate(pattern_merge_hierarchy.root_nodes):
        merged_pattern_name = f"{cluster['cluster_id']}_merged_{i}"
        leaf_nodes = get_leaf_nodes(root_node)

        for leaf in leaf_nodes:
            # Find which original pattern this leaf corresponds to
            pattern_idx = find_pattern_index(leaf.pattern, all_patterns)
            if pattern_idx is not None and pattern_idx < len(all_pattern_idxs):
                original_info = all_pattern_idxs[pattern_idx]

                # Generate logo for original pattern if needed
                input_logo_fwd = None
                input_logo_rev = None

                if len(cluster['components']) <= 50:  # Only generate logos for reasonable numbers
                    # Forward logo
                    logo_fwd_path = f"logos/{original_info['tissue']}_{original_info['pattern_name']}_fwd.png"
                    full_fwd_path = os.path.join(cluster_out_dir, logo_fwd_path)
                    _plot_weights(leaf.pattern.contrib_scores, full_fwd_path)
                    input_logo_fwd = f"./{logo_fwd_path}"

                    # Reverse logo
                    cwm_rev = leaf.pattern.contrib_scores[::-1, ::-1]
                    logo_rev_path = f"logos/{original_info['tissue']}_{original_info['pattern_name']}_rev.png"
                    full_rev_path = os.path.join(cluster_out_dir, logo_rev_path)
                    _plot_weights(cwm_rev, full_rev_path)
                    input_logo_rev = f"./{logo_rev_path}"

                # Get total seqlets for this merged pattern
                total_seqlets_for_merged = merged_pattern_total_seqlets.get(merged_pattern_name, 0) if merged_pattern_total_seqlets else 0

                report_data.append({
                    'merged_pattern': merged_pattern_name,
                    'merged_logo': merged_logo_paths[i],
                    'input_tissue': original_info['tissue'],
                    'input_pattern': original_info['pattern_name'],
                    'input_pattern_full': original_info['pattern_full'],
                    'input_logo_fwd': input_logo_fwd,
                    'input_logo_rev': input_logo_rev,
                    'n_seqlets': len(leaf.pattern.seqlets),
                    'total_seqlets_in_merged': total_seqlets_for_merged,
                    'pattern_class': pattern_class
                })

    # Create DataFrame and save reports
    if report_data:
        df = pd.DataFrame(report_data)

        # Save TSV report
        tsv_path = os.path.join(cluster_out_dir, "merge_report.tsv")
        df[['merged_pattern', 'input_tissue', 'input_pattern', 'input_pattern_full',
            'n_seqlets', 'total_seqlets_in_merged', 'pattern_class']].to_csv(tsv_path, sep='\t', index=False)

        # Save HTML report if we have logos
        if any(df['input_logo_fwd'].notna()):
            html_path = os.path.join(cluster_out_dir, "merge_report.html")
            df.to_html(html_path, escape=False,
                      formatters=dict(merged_logo=path_to_image_html,
                                    input_logo_fwd=path_to_image_html,
                                    input_logo_rev=path_to_image_html),
                      index=False)

        print(f"\tGenerated reports with {len(df)} pattern mappings")

        # Create a summary report with just the merged patterns and their total seqlets
        summary_data = []
        for merged_pattern_name, total_seqlets in merged_pattern_total_seqlets.items():
            summary_data.append({
                'merged_pattern': merged_pattern_name,
                'total_seqlets': total_seqlets
            })

        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_tsv = os.path.join(cluster_out_dir, "merged_patterns_summary.tsv")
            summary_df.to_csv(summary_tsv, sep='\t', index=False)
            print(f"\tGenerated summary report with {len(summary_df)} merged patterns")

def main():
    args = parse_args()

    print("Starting pattern merging across tissues...")
    print(f"Output directory: {args.out_dir}")
    print(f"Modisco directory: {args.modisco_dir}")
    print(f"Contribs directory: {args.contribs_dir}")

    # Load clusters to process
    clusters_to_merge = load_tissues_and_clusters(args)

    # Process each cluster
    for i, cluster in enumerate(clusters_to_merge):
        print(f"\nProcessing cluster {i+1}/{len(clusters_to_merge)}: {cluster['cluster_id']}")
        print(f"Components: {[comp['pattern_full'] for comp in cluster['components']]}")

        try:
            # Load data for this cluster
            data_tuple = load_pattern_data(cluster, args.modisco_dir, args.contribs_dir)

            # Merge patterns
            merge_patterns_for_cluster(cluster, data_tuple, args.out_dir)

            print(f"✓ Completed cluster {cluster['cluster_id']}")

        except Exception as e:
            print(f"✗ ERROR processing cluster {cluster['cluster_id']}: {str(e)}")
            continue

    print(f"\nPattern merging completed. Results saved to {args.out_dir}")

if __name__ == "__main__":
    main()
