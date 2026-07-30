#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.

module load cuda/11.8
module load cudnn/8.9.6.50_cuda11
mkdir -p /data/home/sczd644/run/zsw_chrombpnet/pred_bw
cat  /data/home/sczd644/run/zsw_chrombpnet/tissue.txt |while read id
do
 chrombpnet pred_bw -bm /data/home/sczd644/run/zsw_chrombpnet/${id}/${id}_chrombpnet_model/models/bias_model_scaled.h5 -cm /data/home/sczd644/run/zsw_chrombpnet/${id}/${id}_chrombpnet_model/models/chrombpnet.h5 -cmb /data/home/sczd644/run/zsw_chrombpnet/${id}/${id}_chrombpnet_model/models/chrombpnet_nobias.h5 -r /data/home/sczd644/run/zsw_chrombpnet/${id}/data/peaks_no_blacklist.bed -g /data/home/sczd644/run/zsw_chrombpnet/ref.fa -c /data/home/sczd644/run/zsw_chrombpnet/sheep.chrom.sizes  -op /data/home/sczd644/run/zsw_chrombpnet/pred_bw/${id}_peaks
 done
