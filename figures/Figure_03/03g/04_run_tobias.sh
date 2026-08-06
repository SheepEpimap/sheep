#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

conda activate footprint
conda install -c bioconda tobias

TOBIAS ATACorrect \
  --bam /data/home/sczd644/run/zsw_chrombpnet/ATAC_bams/abomasum/data/merged.bam \
  --genome /data/home/sczd644/run/zsw_chrombpnet/ref.fa \
  --peaks /data/home/sczd644/run/zsw_chrombpnet/ATAC_bams/abomasum/data/peaks_no_blacklist.bed \
  --outdir /data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/tobias \
  --prefix abomasum




#!/bin/bash
tissue=${1}
TOBIAS ATACorrect \
  --bam /data/home/sczd644/run/zsw_chrombpnet/ATAC_bams/${tissue}/data/merged.bam \
  --genome /data/home/sczd644/run/zsw_chrombpnet/ref.fa \
  --peaks /data/home/sczd644/run/zsw_chrombpnet/ATAC_bams/${tissue}/data/peaks_no_blacklist.bed \
  --outdir /data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/tobias \
  --prefix ${tissue}

  cat /data/home/sczd644/run/zsw_chrombpnet/tissue.txt | while read tissue
do
    sbatch -D ./ -c 60 -p low --mem 120G 1.sh $tissue
    sleep 2
done
