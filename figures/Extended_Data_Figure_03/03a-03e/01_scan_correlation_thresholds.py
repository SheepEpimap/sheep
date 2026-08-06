#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-
"""
h3_atac_threshold_scan.py

 (0-1)  H3K27ac   ATAC  ,
 (qval>0.05 & r<0), (ROC AUC, PR AUC, F1, accuracy, precision, recall, specificity, MCC  ).
 output:
 - qvalue  (  qval & fraction q<=0.05)
 -   qval<=0.05 & r>0  , (0.05  ) -
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve,
    average_precision_score, confusion_matrix, classification_report,
    accuracy_score, precision_score, recall_score, matthews_corrcoef, f1_score
)

def load_data(file_path):
    cols = ["gene", "region", "pearson_r", "pval", "distance", "n_samples", "qval"]
    print(f"Loading {file_path} ...")
    df = pd.read_csv(file_path, sep="\t", skiprows=1, names=cols)
    return df

def compute_metrics_from_binary(y_true, y_pred):
    # binary metrics (handles cases where y_true or y_pred may be empty)
    if len(y_true) == 0:
        return dict(accuracy=np.nan, precision=np.nan, recall=np.nan, f1=np.nan, specificity=np.nan, mcc=np.nan)
    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    specificity = tn / (tn + fp + 1e-8)
    mcc = matthews_corrcoef(y_true, y_pred) if (len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1) else np.nan
    return dict(accuracy=acc, precision=prec, recall=rec, f1=f1, specificity=specificity, mcc=mcc)

def safe_auc(y_true, y_score):
    # return ROC AUC and PR AUC safely (handle degenerate cases)
    if len(np.unique(y_true)) < 2:
        return np.nan, np.nan
    try:
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
    except Exception:
        roc_auc = np.nan
    try:
        pr_auc = average_precision_score(y_true, y_score)
    except Exception:
        pr_auc = np.nan
    return roc_auc, pr_auc

def main():
    parser = argparse.ArgumentParser(description="Scan thresholds for H3K27ac vs ATAC validation")
    parser.add_argument("--h3k27ac", required=True, help="H3K27ac correlation file (tab, headerless, same cols as script expects)")
    parser.add_argument("--atac", required=True, help="ATAC correlation file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--n_thresh", type=int, default=201, help="Number of thresholds between 0 and 1 (default 201 -> step 0.005)")
    parser.add_argument("--thresh_min", type=float, default=0.0, help="minimum threshold (default 0.0)")
    parser.add_argument("--thresh_max", type=float, default=1.0, help="maximum threshold (default 1.0)")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    h3 = load_data(args.h3k27ac)
    atac = load_data(args.atac)

    h3 = h3.rename(columns={"pearson_r": "pearson_r_h3", "qval": "qval_h3"})
    atac = atac.rename(columns={"pearson_r": "pearson_r_atac", "qval": "qval_atac"})

    # -------------------------
    # -------------------------
    h3["neg_label_h3"] = np.where((h3["qval_h3"] > 0.05) & (h3["pearson_r_h3"] < 0), 1, 0)
    atac["neg_label_atac"] = np.where((atac["qval_atac"] > 0.05) & (atac["pearson_r_atac"] < 0), 1, 0)

    h3_neg = h3.loc[h3["neg_label_h3"] == 1, ["gene", "region", "pearson_r_h3", "qval_h3"]].copy()
    atac_neg = atac.loc[atac["neg_label_atac"] == 1, ["gene", "region", "pearson_r_atac", "qval_atac"]].copy()

    atac_neg_set = set(zip(atac_neg["gene"], atac_neg["region"]))
    h3_neg_set = set(zip(h3_neg["gene"], h3_neg["region"]))

    # -------------------------
    # -------------------------
    thresholds = np.linspace(args.thresh_min, args.thresh_max, args.n_thresh)
    results = []

    for t in thresholds:
        atac_pos = atac[(atac["pearson_r_atac"] > t) & (atac["qval_atac"] <= 0.05)].copy()
        atac_pos["atac_label"] = 1

        atac_fixed_neg = atac.loc[(atac["qval_atac"] > 0.05) & (atac["pearson_r_atac"] < 0), ["gene", "region"]].copy()
        atac_fixed_neg["atac_label"] = 0

        atac_labels = pd.concat([
            atac_pos[["gene", "region", "atac_label"]],
            atac_fixed_neg[["gene", "region", "atac_label"]]
        ], ignore_index=True).drop_duplicates(subset=["gene", "region"])

        h3_pos = h3[(h3["pearson_r_h3"] > t) & (h3["qval_h3"] <= 0.05)].copy()
        h3_pos["h3_label"] = 1

        h3_fixed_neg = h3.loc[(h3["qval_h3"] > 0.05) & (h3["pearson_r_h3"] < 0), ["gene", "region"]].copy()
        h3_fixed_neg["h3_label"] = 0

        h3_labels = pd.concat([
            h3_pos[["gene", "region", "h3_label", "pearson_r_h3", "qval_h3"]],
            h3_fixed_neg[["gene", "region", "h3_label"]]
        ], ignore_index=True).drop_duplicates(subset=["gene", "region"])

        # Merge on gene/region; we require atac_label exist (y_true) and h3 pearson_r exists to compute scores
        merged = pd.merge(atac_labels, h3_labels, on=["gene", "region"], how="inner")

        if merged.shape[0] == 0:
            results.append({
                "threshold": t, "n_pairs": 0,
                "roc_auc": np.nan, "pr_auc": np.nan,
                "accuracy": np.nan, "precision": np.nan, "recall": np.nan,
                "f1": np.nan, "specificity": np.nan, "mcc": np.nan
            })
            continue

        # Prepare y_true (from atac) and y_score (continuous from h3), and binary prediction from h3_label
        y_true = merged["atac_label"].astype(int).values
        # Use h3 continuous pearson_r as score (where available). If pearson_r_h3 missing for some rows (shouldn't), drop them.
        if "pearson_r_h3" not in merged.columns:
            # In rare case h3 pearson not present (e.g., row came only from fixed neg with no pearson stored), fallback to h3_label as score
            merged["pearson_r_h3"] = merged["h3_label"].astype(float)
        y_score = merged["pearson_r_h3"].astype(float).values
        y_pred = merged["h3_label"].astype(int).values

        # Compute AUCs
        roc_auc, pr_auc = safe_auc(y_true, y_score)

        # Binary metrics
        bin_metrics = compute_metrics_from_binary(y_true, y_pred)

        # Save result
        results.append({
            "threshold": t,
            "n_pairs": merged.shape[0],
            "roc_auc": roc_auc,
            "pr_auc": pr_auc,
            "accuracy": bin_metrics["accuracy"],
            "precision": bin_metrics["precision"],
            "recall": bin_metrics["recall"],
            "f1": bin_metrics["f1"],
            "specificity": bin_metrics["specificity"],
            "mcc": bin_metrics["mcc"]
        })

    df_results = pd.DataFrame(results)
    df_results.to_csv(os.path.join(args.output, "threshold_scan_metrics.tsv"), sep="\t", index=False)

    # -------------------------
    # -------------------------
    plt.figure(figsize=(8,5))
    plt.plot(df_results["threshold"], df_results["roc_auc"], label="ROC AUC")
    plt.plot(df_results["threshold"], df_results["pr_auc"], label="PR AUC")
    plt.xlabel("Correlation threshold")
    plt.ylabel("AUC")
    plt.title("AUC vs correlation threshold")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "auc_vs_threshold.pdf"), bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10,6))
    plt.plot(df_results["threshold"], df_results["f1"], label="F1")
    plt.plot(df_results["threshold"], df_results["accuracy"], label="Accuracy")
    plt.plot(df_results["threshold"], df_results["precision"], label="Precision")
    plt.plot(df_results["threshold"], df_results["recall"], label="Recall")
    plt.plot(df_results["threshold"], df_results["mcc"], label="MCC")
    plt.xlabel("Correlation threshold")
    plt.ylabel("Metric value")
    plt.title("Metrics vs correlation threshold")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "metrics_vs_threshold.pdf"), bbox_inches='tight')
    plt.close()

    # -------------------------
    # -------------------------
    q_trends = []
    thresh_for_q = np.linspace(0.0, 1.0, 201)
    for t in thresh_for_q:
        h3_sel = h3[h3["pearson_r_h3"] >= t]
        atac_sel = atac[atac["pearson_r_atac"] >= t]

        def summarize_q(df, qcol):
            if df.shape[0] == 0:
                return dict(mean_q=np.nan, frac_q_le_0_05=np.nan, n=df.shape[0])
            return dict(mean_q=df[qcol].mean(), frac_q_le_0_05=(df[qcol] <= 0.05).mean(), n=df.shape[0])

        h3_sum = summarize_q(h3_sel, "qval_h3")
        atac_sum = summarize_q(atac_sel, "qval_atac")
        q_trends.append({
            "threshold": t,
            "h3_mean_q": h3_sum["mean_q"],
            "h3_frac_q_le_0.05": h3_sum["frac_q_le_0_05"],
            "h3_n": h3_sum["n"],
            "atac_mean_q": atac_sum["mean_q"],
            "atac_frac_q_le_0.05": atac_sum["frac_q_le_0_05"],
            "atac_n": atac_sum["n"]
        })

    df_qtrends = pd.DataFrame(q_trends)
    df_qtrends.to_csv(os.path.join(args.output, "qvalue_vs_corr_threshold.tsv"), sep="\t", index=False)

    plt.figure(figsize=(8,5))
    plt.plot(df_qtrends["threshold"], df_qtrends["h3_mean_q"], label="H3 mean q")
    plt.plot(df_qtrends["threshold"], df_qtrends["atac_mean_q"], label="ATAC mean q")
    plt.xlabel("Correlation threshold")
    plt.ylabel("Mean qvalue")
    plt.title("Mean qvalue vs correlation threshold")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "mean_q_vs_threshold.pdf"), bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(df_qtrends["threshold"], df_qtrends["h3_frac_q_le_0.05"], label="H3 fraction q<=0.05")
    plt.plot(df_qtrends["threshold"], df_qtrends["atac_frac_q_le_0.05"], label="ATAC fraction q<=0.05")
    plt.xlabel("Correlation threshold")
    plt.ylabel("Fraction (q <= 0.05)")
    plt.title("Fraction of entries with q <= 0.05 vs correlation threshold")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "frac_qle0.05_vs_threshold.pdf"), bbox_inches='tight')
    plt.close()

    # -------------------------
    # -------------------------
    base_h3 = h3[(h3["qval_h3"] <= 0.05) & (h3["pearson_r_h3"] > 0)].copy()
    base_atac = atac[(atac["qval_atac"] <= 0.05) & (atac["pearson_r_atac"] > 0)].copy()

    # density plot of pearson_r under base filter (H3 and ATAC separately)
    plt.figure(figsize=(8,5))
    if base_h3.shape[0] > 0:
        sns.kdeplot(base_h3["pearson_r_h3"], label=f"H3 (n={len(base_h3)})")
    if base_atac.shape[0] > 0:
        sns.kdeplot(base_atac["pearson_r_atac"], label=f"ATAC (n={len(base_atac)})")
    plt.xlabel("Pearson r")
    plt.title("Density of Pearson r (q<=0.05 & r>0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "density_corr_qle0.05_rpos.pdf"), bbox_inches='tight')
    plt.close()

    # counts at thresholds 0, 0.05, 0.10, ... 1.0
    steps = np.arange(0.0, 1.0001, 0.05)
    counts = []
    for thr in steps:
        h3_count = base_h3[base_h3["pearson_r_h3"] >= thr].shape[0]
        atac_count = base_atac[base_atac["pearson_r_atac"] >= thr].shape[0]
        counts.append({"threshold": thr, "h3_count": h3_count, "atac_count": atac_count})

    df_counts = pd.DataFrame(counts)
    df_counts.to_csv(os.path.join(args.output, "counts_qle0.05_rpos_by_step0.05.tsv"), sep="\t", index=False)

    # bar plot of counts
    width = 0.4
    x = np.arange(len(df_counts))
    plt.figure(figsize=(12,6))
    plt.bar(x - width/2, df_counts["h3_count"], width=width, label="H3 counts")
    plt.bar(x + width/2, df_counts["atac_count"], width=width, label="ATAC counts")
    plt.xticks(x, [f"{thr:.2f}" for thr in df_counts["threshold"]], rotation=45)
    plt.xlabel("Threshold (step 0.05)")
    plt.ylabel("Number of enhancer-target pairs (q<=0.05 & r>0)")
    plt.title("Counts vs threshold (step 0.05)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "counts_bar_qle0.05_rpos_step0.05.pdf"), bbox_inches='tight')
    plt.close()

    # save merged results for best-threshold example (optional)
    # find threshold with max F1
    best_row = df_results.loc[df_results["f1"].idxmax()] if df_results["f1"].notna().any() else None
    if best_row is not None:
        best_t = best_row["threshold"]
        # recompute merged for best_t
        atac_pos = atac[(atac["pearson_r_atac"] > best_t) & (atac["qval_atac"] <= 0.05)].copy()
        atac_pos["atac_label"] = 1
        atac_fixed_neg = atac.loc[(atac["qval_atac"] > 0.05) & (atac["pearson_r_atac"] < 0), ["gene", "region"]].copy()
        atac_fixed_neg["atac_label"] = 0
        atac_labels = pd.concat([atac_pos[["gene","region","atac_label"]], atac_fixed_neg[["gene","region","atac_label"]]], ignore_index=True).drop_duplicates(subset=["gene","region"])

        h3_pos = h3[(h3["pearson_r_h3"] > best_t) & (h3["qval_h3"] <= 0.05)].copy()
        h3_pos["h3_label"] = 1
        h3_fixed_neg = h3.loc[(h3["qval_h3"] > 0.05) & (h3["pearson_r_h3"] < 0), ["gene", "region"]].copy()
        h3_fixed_neg["h3_label"] = 0
        h3_labels = pd.concat([h3_pos[["gene","region","h3_label","pearson_r_h3","qval_h3"]], h3_fixed_neg[["gene","region","h3_label"]]], ignore_index=True).drop_duplicates(subset=["gene","region"])

        merged_best = pd.merge(atac_labels, h3_labels, on=["gene","region"], how="inner")
        merged_best.to_csv(os.path.join(args.output, f"merged_best_threshold_{best_t:.3f}.tsv"), sep="\t", index=False)

    print("All done. Outputs written to:", args.output)

if __name__ == "__main__":
    main()
