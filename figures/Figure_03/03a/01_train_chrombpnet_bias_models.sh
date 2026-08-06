#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
module load cuda/11.8
module load cudnn/8.9.6.50_cuda11
cat /data/home/sczd644/run/zsw_chrombpnet/bia_altiss.txt |while read id
do
chrombpnet bias pipeline \
        -ibam /data/home/sczd644/run/zsw_chrombpnet/${id}/data/merged.bam \
        -d "ATAC" \
        -g /data/home/sczd644/run/zsw_chrombpnet/ref.fa \
        -c /data/home/sczd644/run/zsw_chrombpnet/sheep.chrom.sizes  \
        -p /data/home/sczd644/run/zsw_chrombpnet/${id}/data/peaks_no_blacklist.bed \
        -n /data/home/sczd644/run/zsw_chrombpnet/${id}/data/output_negatives.bed \
        -fl /data/home/sczd644/run/zsw_chrombpnet/${id}/data/splits/fold_0.json \
        -b 0.5 \
        -o /data/home/sczd644/run/zsw_chrombpnet/${id}/data/bias_model/ \
        -fp ${id}
done
