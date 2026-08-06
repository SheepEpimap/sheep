#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
python plot_topN_links_to_pdfs.py \
  --links_tsv /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_output_E5_confident.tsv \
  --enh_cpm  /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_E5_counts_cpms.csv \
  --expr39   /vol2/zhangshiwen/predict/Expression_tpm_39.tsv \
  --expr40   /vol2/zhangshiwen/predict/Expression_tpm_40.tsv \
  --tissue_colors /data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/summary/tissue_colors.tsv \
  --out_dir /vol2/zhangshiwen/sheep_cor/h3k27ac/E5_top8_noLOC \
  --prefix E5_top8_noLOC \
  --topN 8 \
  --exclude_loc \
  --label_top_x 2 \
  --label_top_y 2 \
  --point_size 40
