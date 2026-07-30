#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
python calc_p1t_q1t_p2_p3_ratios.py \
  --tissue-list /data/home/sczd644/run/zsw_chrombpnet/tissue.txt \
  --vcf-tar /vol2/00_panlab_rawdata/01_gwasdata/WGS_235raw/35_high_depth/35_sub_indv.vcf.tar.gz \
  --as-in-peak-dir /data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched \
  --nonas-in-peak-dir /data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_peak \
  --ann-dir /data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched/AS/02 \
  --finemo-base /data/home/sczd644/run/zsw_chrombpnet/finemo \
  --logfc-pval-threshold 0.05 \
  --motif-mode either \
  --outdir /data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched/p1t_q1t_p2_p3_ratios_strict_peak
