#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
mkdir -p /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_contribs
module load cuda/11.8
module load cudnn/8.9.6.50_cuda11
cat /data/home/sczd644/run/zsw_chrombpnet/tissue.txt |while read id
do
 chrombpnet contribs_bw \
    -m /data/home/sczd644/run/zsw_chrombpnet/${id}/${id}_chrombpnet_model/models/chrombpnet_nobias.h5 \
    -r /data/home/sczd644/run/zsw_chrombpnet/${id}/data/peaks_no_blacklist.bed \
    -g /data/home/sczd644/run/zsw_chrombpnet/ref.fa \
    -c /data/home/sczd644/run/zsw_chrombpnet/sheep.chrom.sizes \
    -op /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_contribs/${id}_chrombpnet_contribs

    mkdir -p /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_motifs/${id}_chrombpnet_motifs
modisco motifs -i /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_contribs/${id}_chrombpnet_contribs.counts_scores.h5 -n 1000000 -o /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_motifs/${id}_chrombpnet_motifs/${id}.modisco_results.h5


modisco report -i /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_motifs/${id}_chrombpnet_motifs/${id}.modisco_results.h5 -o /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_motifs/${id}_chrombpnet_motifs -m /data/home/sczd644/run/zsw_chrombpnet/JASPAR2024_CORE_redundant_pfms_meme.txt

done
