#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
"""
Purpose: Generate report for compiled/merged patterns with CWM logos and TomTom matches
"""

import sys
import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import modiscolite.report
import types
from typing import List, Union
import tempfile
import shutil
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description='Generate report for compiled modisco patterns with TomTom matches.')
    parser.add_argument('--compiled-h5', type=str, required=True,
                       help='Path to the compiled modisco h5 file (modisco_compiled.h5).')
    parser.add_argument('--out-dir', type=str, required=True,
                       help='Path to the output directory.')
    parser.add_argument('--meme-db', type=str, required=True,
                       help='Path to the motifs database in meme format.')
    parser.add_argument('--trim-threshold', type=float, default=0.3,
                       help='Probability threshold for trimming the PPM.')
    parser.add_argument('--trim-min-length', type=int, default=3,
                       help='Minimum length of the trimmed PPM.')
    parser.add_argument('--n-matches', type=int, default=10,
                       help='Number of top matches to fetch from tomtom.')
    parser.add_argument('--top-n-matches', type=int, default=3,
                       help='Number of top matches whose logos to include in the report.')
    parser.add_argument('--pattern-groups', type=str, nargs='+', default=['pos_patterns'],
                       help='Pattern groups to process (e.g., pos_patterns neg_patterns)')
    parser.add_argument('--verbose', action='store_true', default=False,
                       help='If True, print more info.')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.isfile(args.compiled_h5):
        raise ValueError(f'Invalid path to the compiled modisco h5 file: {args.compiled_h5}')
    if not os.path.isfile(args.meme_db):
        raise ValueError(f'Invalid path to the meme database: {args.meme_db}')

    os.makedirs(args.out_dir, exist_ok=True)

    logger.info(f"Arguments: {args}")
    return args

def fetch_tomtom_matches(ppm, cwm, pattern_name, motifs_db, out_dir,
                        background=[0.25, 0.25, 0.25, 0.25],
                        tomtom_exec_path='tomtom',
                        trim_threshold=0.3, trim_min_length=3,
                        save_tomtom_output=True):
    """
    Fetches top matches from a motifs database using TomTom.

    Args:
        ppm: position probability matrix - numpy matrix of dimension (N,4)
        cwm: contribution weight matrix - numpy matrix of dimension (N,4)
        pattern_name: the name of the pattern
        motifs_db: path to motifs database in meme format
        out_dir: output directory
        background: list with ACGT background probabilities
        tomtom_exec_path: path to TomTom executable
        trim_threshold: probability threshold for trimming
        trim_min_length: minimum length after trimming
        save_tomtom_output: whether to save TomTom output files

    Returns:
        pandas.DataFrame: TomTom results
    """

    # Create temporary files
    _, fname = tempfile.mkstemp()
    _, tomtom_fname = tempfile.mkstemp()

    try:
        # Trim the pattern based on CWM scores
        score = np.sum(np.abs(cwm), axis=1)
        trim_thresh = np.max(score) * trim_threshold
        pass_inds = np.where(score >= trim_thresh)[0]

        if len(pass_inds) == 0:
            logger.warning(f"No positions pass trimming threshold for {pattern_name}")
            return pd.DataFrame()

        # Ensure minimum length
        if len(pass_inds) < trim_min_length:
            # Use the top trim_min_length positions by score
            top_inds = np.argsort(score)[-trim_min_length:]
            pass_inds = np.sort(top_inds)

        trimmed = ppm[np.min(pass_inds): np.max(pass_inds) + 1]

        # Write trimmed motif to meme file
        modiscolite.report.write_meme_file(trimmed, background, fname)

        # Check if tomtom is available
        if not shutil.which(tomtom_exec_path):
            raise ValueError(f'`tomtom` executable not found. Please install MEME suite.')

        # Run tomtom
        cmd = f'{tomtom_exec_path} -no-ssc -oc . --verbosity 1 -text -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 {fname} {motifs_db} > {tomtom_fname}'
        os.system(cmd)

        # Read tomtom results
        if os.path.getsize(tomtom_fname) > 0:
            tomtom_results = pd.read_csv(tomtom_fname, sep="\t")
            # Select relevant columns
            if not tomtom_results.empty:
                required_cols = ['Target_ID', 'p-value', 'E-value', 'q-value', 'Query_consensus']
                available_cols = [col for col in required_cols if col in tomtom_results.columns]
                tomtom_results = tomtom_results[available_cols]
        else:
            tomtom_results = pd.DataFrame()

        # Save tomtom output if requested
        if save_tomtom_output and not tomtom_results.empty:
            output_subdir = os.path.join(out_dir, "tomtom")
            os.makedirs(output_subdir, exist_ok=True)
            output_filepath = os.path.join(output_subdir, f"{pattern_name}.tomtom.tsv")
            tomtom_results.to_csv(output_filepath, sep="\t", index=False)
            logger.info(f"Saved TomTom results to {output_filepath}")

    except Exception as e:
        logger.error(f"Error running TomTom for {pattern_name}: {str(e)}")
        tomtom_results = pd.DataFrame()
    finally:
        # Clean up temporary files
        if os.path.exists(fname):
            os.remove(fname)
        if os.path.exists(tomtom_fname) and not save_tomtom_output:
            os.remove(tomtom_fname)

    return tomtom_results

def main(args):
    """Main function to generate report for compiled modisco patterns."""

    # Unpack arguments
    compiled_h5 = args.compiled_h5
    out_dir = args.out_dir
    meme_db = args.meme_db
    trim_threshold = args.trim_threshold
    n_matches = args.n_matches
    top_n_matches = args.top_n_matches
    pattern_groups = args.pattern_groups
    verbose = args.verbose

    logger.info("Starting report generation for compiled modisco patterns")

    # Create output directories
    modisco_logo_dir = os.path.join(out_dir, 'trimmed_logos')
    db_logo_dir = os.path.join(out_dir, 'db_logos')
    tomtom_dir = os.path.join(out_dir, 'tomtom')

    os.makedirs(modisco_logo_dir, exist_ok=True)
    os.makedirs(db_logo_dir, exist_ok=True)
    os.makedirs(tomtom_dir, exist_ok=True)

    # Step 1: Create trimmed logos for modisco patterns using modiscolite function
    logger.info("Creating trimmed logos for modisco patterns...")
    logger.info(f"Processing pattern groups: {pattern_groups}")

    try:
        # Use the built-in modiscolite function
        modiscolite.report.create_modisco_logos(
            compiled_h5,
            modisco_logo_dir,
            trim_threshold,
            pattern_groups
        )
        logger.info("Successfully created modisco logos")
    except Exception as e:
        logger.error(f"Error creating modisco logos: {str(e)}")
        # Fallback: try with default parameters if specific call fails
        try:
            modiscolite.report.create_modisco_logos(compiled_h5, modisco_logo_dir)
            logger.info("Successfully created modisco logos with default parameters")
        except Exception as e2:
            logger.error(f"Failed to create logos even with default parameters: {str(e2)}")
            return

    # Step 2: Initialize results dataframe
    logger.info("Initializing results dataframe...")
    results_data = {
        'pattern': [],
        'pattern_group': [],
        'modisco_cwm_fwd': [],
        'modisco_cwm_rev': [],
        'query_consensus': []
    }

    # Add columns for TomTom matches
    for i in range(n_matches):
        results_data[f'match_{i}'] = []
        results_data[f'e_value_{i}'] = []
        results_data[f'q_value_{i}'] = []
        results_data[f'p_value_{i}'] = []

    # Step 3: Load motif database
    logger.info(f"Loading motif database from {meme_db}...")
    try:
        motifs = modiscolite.report.read_meme(meme_db)
        logger.info(f"Loaded {len(motifs)} motifs from database")
    except Exception as e:
        logger.error(f"Error loading motif database: {str(e)}")
        return

    # Step 4: Process each pattern and run TomTom
    logger.info("Processing patterns and running TomTom...")

    with h5.File(compiled_h5, 'r') as f:
        for pattern_group in pattern_groups:
            if pattern_group not in f:
                logger.warning(f"Pattern group '{pattern_group}' not found in H5 file, skipping")
                continue

            metacluster = f[pattern_group]
            pattern_names = list(metacluster.keys())

            logger.info(f"Found {len(pattern_names)} patterns in group '{pattern_group}'")

            for pattern_name in pattern_names:
                if verbose:
                    logger.info(f"Processing pattern: {pattern_group}/{pattern_name}")

                pattern_data = metacluster[pattern_name]

                # Get PPM and CWM
                ppm = np.array(pattern_data['sequence'][:])
                cwm = np.array(pattern_data['contrib_scores'][:])

                # Add basic pattern info to results
                results_data['pattern'].append(pattern_name)
                results_data['pattern_group'].append(pattern_group)
                results_data['modisco_cwm_fwd'].append(
                    f"./trimmed_logos/{pattern_name}.cwm.fwd.png"
                )
                results_data['modisco_cwm_rev'].append(
                    f"./trimmed_logos/{pattern_name}.cwm.rev.png"
                )

                # Run TomTom
                tomtom_results = fetch_tomtom_matches(
                    ppm, cwm, pattern_name, meme_db, out_dir,
                    trim_threshold=trim_threshold
                )

                # Process TomTom results
                if not tomtom_results.empty:
                    # Get query consensus from first result
                    if 'Query_consensus' in tomtom_results.columns:
                        results_data['query_consensus'].append(
                            tomtom_results.iloc[0]['Query_consensus']
                        )
                    else:
                        results_data['query_consensus'].append('')

                    # Add top matches
                    for i in range(min(n_matches, len(tomtom_results))):
                        row = tomtom_results.iloc[i]
                        results_data[f'match_{i}'].append(row.get('Target_ID', ''))
                        results_data[f'e_value_{i}'].append(row.get('E-value', ''))
                        results_data[f'q_value_{i}'].append(row.get('q-value', ''))
                        results_data[f'p_value_{i}'].append(row.get('p-value', ''))

                    # Fill remaining slots with empty values
                    for i in range(len(tomtom_results), n_matches):
                        results_data[f'match_{i}'].append('')
                        results_data[f'e_value_{i}'].append('')
                        results_data[f'q_value_{i}'].append('')
                        results_data[f'p_value_{i}'].append('')
                else:
                    # No TomTom results
                    results_data['query_consensus'].append('')
                    for i in range(n_matches):
                        results_data[f'match_{i}'].append('')
                        results_data[f'e_value_{i}'].append('')
                        results_data[f'q_value_{i}'].append('')
                        results_data[f'p_value_{i}'].append('')

    # Step 5: Create dataframe
    logger.info("Creating results dataframe...")
    results_df = pd.DataFrame(results_data)

    if results_df.empty:
        logger.warning("No patterns were processed, exiting.")
        return

    # Step 6: Generate logos for top matched motifs
    logger.info(f"Generating logos for top {top_n_matches} matches...")

    for i in range(top_n_matches):
        logo_paths = []
        match_col = f'match_{i}'

        if match_col not in results_df.columns:
            continue

        for _, row in results_df.iterrows():
            motif_id = row[match_col]
            if pd.notna(motif_id) and motif_id != '':
                try:
                    # Generate logo for matched motif
                    modiscolite.report.make_logo(motif_id, db_logo_dir, motifs)
                    logo_paths.append(f"./db_logos/{motif_id}.png")
                except Exception as e:
                    logger.warning(f"Could not generate logo for {motif_id}: {str(e)}")
                    logo_paths.append('')
            else:
                logo_paths.append('')

        results_df[f'{match_col}_logo'] = logo_paths

    # Step 7: Save results
    logger.info("Saving final results...")
    output_tsv = os.path.join(out_dir, "compiled_modisco_report.tsv")
    results_df.to_csv(output_tsv, sep="\t", index=False)

    # Also save a more readable version with only top matches
    top_columns = ['pattern', 'pattern_group', 'modisco_cwm_fwd', 'modisco_cwm_rev', 'query_consensus']
    for i in range(top_n_matches):
        top_columns.extend([
            f'match_{i}', f'e_value_{i}', f'q_value_{i}',
            f'p_value_{i}', f'match_{i}_logo'
        ])

    if all(col in results_df.columns for col in top_columns):
        results_df[top_columns].to_csv(
            os.path.join(out_dir, "compiled_modisco_report_top_matches.tsv"),
            sep="\t", index=False
        )

    logger.info(f"Report generation completed! Results saved to {out_dir}")
    logger.info(f"Main report: {output_tsv}")
    logger.info(f"Number of patterns processed: {len(results_df)}")

if __name__ == '__main__':
    args = parse_args()
    main(args)
