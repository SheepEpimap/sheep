#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

cat /data/home/sczd644/run/zsw_chrombpnet/tissue.txt | while read tissue
do
    sbatch -D ./ -c 60 -p low --mem 120G 1.sh $tissue
    sleep 2
done
