#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

#h3k27ac
#!/bin/bash
cd /vol2/zhangshiwen/sheep_cor/h3k27ac
for state in E6 E7 E8
do
    python /vol2/zhangshiwen/predict/egcorr_downsample_1.py \
    /vol2/zhangshiwen/predict/Expression_tpm_{39,40}.tsv \
    /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_${state}_counts_cpms.csv  \
     /vol2/mengzhu/genome/part_change_esemb100/TSS_esemble100_colin1.bed \
     /vol2/zhangshiwen/hic/tibetan_sheep_TAD.bed \
   /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_output_${state}.tsv

awk '
BEGIN{OFS="\t"}
{
    gene       = $1
    regulator  = $2
    pearson_r  = $3
    pval       = $4
    distance   = $5
    n_samples  = $6
    qval       = $7
}
$7<=0.05 && $3>0.3 {
    print regulator, pearson_r, pval, gene, distance, qval
}
' /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_output_${state}.tsv \
> /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_output_${state}_confident.tsv
done
