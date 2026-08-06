#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
python plot_region_obs_pred_counts_sum.py \
  --region chr13:3295301-3296301    \
  --feature SNAP25 \
  --label-top-n 2 \
  --outdir /data/home/sczd644/run/zsw_chrombpnet/pred_bw/correlation_analysis/region_counts
