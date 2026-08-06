#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

python plot_all_tissues_BCD_vertical_topn.py \
  --hits_root /data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/global_bg_loop_with_pval \
  --annotation_dir /data/home/sczd644/run/zsw_chrombpnet/region_ann \
  --tss_bed /data/home/sczd644/run/zsw_chrombpnet/region_ann/TSS_esemble100_colin.bed \
  --tissue_list /data/home/sczd644/run/zsw_chrombpnet/tissue.txt \
  --p_cut 0.05 \
  --out_prefix sheep_all_topn \
  --exclude_unresolved
  --top_n 50
