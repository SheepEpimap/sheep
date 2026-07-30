#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
mkdir /vol2/zhangshiwen/sheep_cor/h3k27ac
cd /vol2/zhangshiwen/sheep_cor/
for state in E5
do

multiBamSummary BED-file --BED ${state}_no_blacklist_Gs.bed --bamfiles /vol2/mengzhu/snakemake_sheep/clean/bam1/H3K27ac_*.bowtie2.mapped.filtered.sort.bam -o /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_${state}_counts.npz -p 24 --outRawCounts /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_${state}_counts.tsv -e 200

#change the table
python /vol2/mengzhu/TAD_target/prepare_table.py /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_${state}_counts.tsv > /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_${state}_counts_prepared.tsv
#nomalization
Rscript /vol2/mengzhu/TAD_target/edger_norm.R /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_${state}_counts_prepared.tsv /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_${state}_counts_cpms.csv
done
