#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
import json
import os
import pandas as pd

def extract_metrics_from_json(json_file_path):
    """ JSONfile 5 """
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)

        counts_metrics = data.get("counts_metrics", {}).get("peaks", {})
        profile_metrics = data.get("profile_metrics", {}).get("peaks", {})

        return {
            "spearmanr": counts_metrics.get("spearmanr"),
            "pearsonr": counts_metrics.get("pearsonr"),
            "mse": counts_metrics.get("mse"),
            "median_jsd": profile_metrics.get("median_jsd"),
            "median_norm_jsd": profile_metrics.get("median_norm_jsd")
        }
    except Exception as e:
        print(f"Error reading {json_file_path}: {e}")
        return None

def main():
    tissue_file = "/data/home/sczd644/run/zsw_chrombpnet/all_tissue.txt"
    output_dir = "/data/home/sczd644/run/zsw_chrombpnet/qc"
    output_file = os.path.join(output_dir, "chrombpnet_metrics_summary.tsv")

    os.makedirs(output_dir, exist_ok=True)

    try:
        with open(tissue_file, 'r') as f:
            tissues = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading tissue file: {e}")
        return

    results = []

    for tissue in tissues:
        json_file_path = f"/data/home/sczd644/run/zsw_chrombpnet/chrombpnet_model/{tissue}_chrombpnet_model/evaluation/chrombpnet_metrics.json"

        print(f"Processing {tissue}...")

        if os.path.exists(json_file_path):
            metrics = extract_metrics_from_json(json_file_path)
            if metrics:
                metrics["tissue"] = tissue
                results.append(metrics)
            else:
                print(f"  Failed to extract metrics from {json_file_path}")
        else:
            print(f"  File not found: {json_file_path}")

    if results:
        df = pd.DataFrame(results)
        columns = ["tissue", "spearmanr", "pearsonr", "mse", "median_jsd", "median_norm_jsd"]
        df = df[columns]

        df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
        print(f"\nResults saved to: {output_file}")
        print(f"Processed {len(results)} tissues out of {len(tissues)}")

        print("\nPreview of results:")
        print(df.head())
    else:
        print("No results to save.")

if __name__ == "__main__":
    main()
