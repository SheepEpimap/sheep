#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(file_path):
    cols = ["gene", "region", "pearson_r", "pval", "distance", "n_samples", "qval"]
    print(f"Loading {file_path}...")
    return pd.read_csv(file_path, sep="\t", skiprows=1, names=cols)

def main():
    parser = argparse.ArgumentParser(description="Calculate Jaccard similarity between H3K27ac and ATAC enhancer-target pairs at different correlation thresholds")
    parser.add_argument("--h3k27ac", required=True, help="H3K27ac correlation file")
    parser.add_argument("--atac", required=True, help="ATAC correlation file")
    parser.add_argument("--output", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    h3 = load_data(args.h3k27ac)
    atac = load_data(args.atac)

    print("Preprocessing data...")

    h3_filtered = h3[h3["qval"] <= 0.05].copy()
    atac_filtered = atac[atac["qval"] <= 0.05].copy()

    h3_filtered['id'] = h3_filtered['gene'] + '|' + h3_filtered['region']
    atac_filtered['id'] = atac_filtered['gene'] + '|' + atac_filtered['region']

    print(f"H3K27ac filtered: {len(h3_filtered)} pairs")
    print(f"ATAC filtered: {len(atac_filtered)} pairs")

    thresholds = np.arange(0.0, 1.01, 0.01)  #  0 1, 0.01
    jaccard_scores = []
    h3_sizes = []
    atac_sizes = []
    overlap_sizes = []

    print("Calculating Jaccard similarity at different thresholds...")

    for threshold in thresholds:
        h3_above_threshold = h3_filtered[h3_filtered["pearson_r"] > threshold]
        atac_above_threshold = atac_filtered[atac_filtered["pearson_r"] > threshold]

        h3_set = set(h3_above_threshold['id'])
        atac_set = set(atac_above_threshold['id'])

        intersection = h3_set.intersection(atac_set)
        union = h3_set.union(atac_set)

        if len(union) > 0:
            jaccard = len(intersection) / len(union)
        else:
            jaccard = 0

        jaccard_scores.append(jaccard)
        h3_sizes.append(len(h3_set))
        atac_sizes.append(len(atac_set))
        overlap_sizes.append(len(intersection))

    best_idx = np.argmax(jaccard_scores)
    best_threshold = thresholds[best_idx]
    best_jaccard = jaccard_scores[best_idx]

    # ------------------- 6. outputresults -------------------
    results_df = pd.DataFrame({
        'threshold': thresholds,
        'jaccard': jaccard_scores,
        'h3_size': h3_sizes,
        'atac_size': atac_sizes,
        'overlap_size': overlap_sizes
    })

    results_df.to_csv(os.path.join(args.output, "jaccard_similarity_results.csv"), index=False)

    with open(os.path.join(args.output, "best_threshold.txt"), "w") as f:
        f.write(f"Best threshold: {best_threshold:.3f}\n")
        f.write(f"Best Jaccard similarity: {best_jaccard:.4f}\n")
        f.write(f"H3K27ac pairs at best threshold: {h3_sizes[best_idx]}\n")
        f.write(f"ATAC pairs at best threshold: {atac_sizes[best_idx]}\n")
        f.write(f"Overlap pairs at best threshold: {overlap_sizes[best_idx]}\n")

    print("Generating plots...")

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    axes[0, 0].plot(thresholds, jaccard_scores, 'b-', linewidth=2)
    axes[0, 0].axvline(best_threshold, color='red', linestyle='--',
                      label=f'Best threshold: {best_threshold:.2f}\nJaccard: {best_jaccard:.3f}')
    axes[0, 0].set_xlabel('Correlation Threshold')
    axes[0, 0].set_ylabel('Jaccard Similarity')
    axes[0, 0].set_title('Jaccard Similarity vs Correlation Threshold')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()

    axes[0, 1].plot(thresholds, h3_sizes, 'b-', label='H3K27ac', linewidth=2)
    axes[0, 1].plot(thresholds, atac_sizes, 'r-', label='ATAC', linewidth=2)
    axes[0, 1].plot(thresholds, overlap_sizes, 'g-', label='Overlap', linewidth=2)
    axes[0, 1].set_xlabel('Correlation Threshold')
    axes[0, 1].set_ylabel('Number of Pairs')
    axes[0, 1].set_title('Number of Enhancer-Target Pairs vs Threshold')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()

    h3_overlap_ratio = [overlap/h3 if h3 > 0 else 0 for overlap, h3 in zip(overlap_sizes, h3_sizes)]
    atac_overlap_ratio = [overlap/atac if atac > 0 else 0 for overlap, atac in zip(overlap_sizes, atac_sizes)]

    axes[1, 0].plot(thresholds, h3_overlap_ratio, 'b-', label='Overlap/H3K27ac', linewidth=2)
    axes[1, 0].plot(thresholds, atac_overlap_ratio, 'r-', label='Overlap/ATAC', linewidth=2)
    axes[1, 0].set_xlabel('Correlation Threshold')
    axes[1, 0].set_ylabel('Overlap Ratio')
    axes[1, 0].set_title('Overlap Ratio vs Correlation Threshold')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend()

    jaccard_derivative = np.gradient(jaccard_scores, thresholds)
    axes[1, 1].plot(thresholds, jaccard_derivative, 'purple', linewidth=2)
    axes[1, 1].axhline(0, color='gray', linestyle='--')
    axes[1, 1].set_xlabel('Correlation Threshold')
    axes[1, 1].set_ylabel('d(Jaccard)/d(Threshold)')
    axes[1, 1].set_title('Rate of Change of Jaccard Similarity')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "jaccard_analysis.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, jaccard_scores, 'b-', linewidth=3)
    plt.axvline(best_threshold, color='red', linestyle='--',
               label=f'Best threshold: {best_threshold:.2f}\nJaccard: {best_jaccard:.3f}')
    plt.xlabel('Correlation Threshold')
    plt.ylabel('Jaccard Similarity')
    plt.title('Jaccard Similarity between H3K27ac and ATAC Enhancer-Target Pairs')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(args.output, "jaccard_similarity.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, h3_overlap_ratio, 'b-', label='Overlap/H3K27ac', linewidth=2)
    plt.plot(thresholds, atac_overlap_ratio, 'r-', label='Overlap/ATAC', linewidth=2)
    plt.xlabel('Correlation Threshold')
    plt.ylabel('Overlap Ratio')
    plt.title('Overlap Ratio vs Correlation Threshold')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(args.output, "overlap_ratio.png"), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, h3_sizes, 'b-', label='H3K27ac', linewidth=2)
    plt.plot(thresholds, atac_sizes, 'r-', label='ATAC', linewidth=2)
    plt.plot(thresholds, overlap_sizes, 'g-', label='Overlap', linewidth=2)
    plt.xlabel('Correlation Threshold')
    plt.ylabel('Number of Pairs')
    plt.title('Number of Enhancer-Target Pairs vs Threshold')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(args.output, "pair_counts.png"), dpi=300)
    plt.close()

    print("Analysis completed successfully!")
    print(f"Best threshold: {best_threshold:.3f}")
    print(f"Best Jaccard similarity: {best_jaccard:.4f}")

if __name__ == "__main__":
    main()
