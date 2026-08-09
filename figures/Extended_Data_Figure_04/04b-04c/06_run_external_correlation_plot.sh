#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

python 1.py \
  --corfile /data/home/sczd644/run/zsw_chrombpnet/public/pred_bw/correlation_analysis/all_tissues_correlation_results.txt \
  --colors  /data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv \
  --outpdf  /data/home/sczd644/run/zsw_chrombpnet/public/pred_bw/correlation_analysis/tissue_cor.break_0_0.5.pdf \
  --outsvg  /data/home/sczd644/run/zsw_chrombpnet/public/pred_bw/correlation_analysis/tissue_cor.break_0_0.5.svg \
  --mode tissuepair \
  --label-col pair \
  --break-low 0.05 --break-high 0.50 \
  --sort desc
