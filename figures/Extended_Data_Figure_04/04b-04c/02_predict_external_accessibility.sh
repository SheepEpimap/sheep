#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

TISSUE_TXT="/data/home/sczd644/run/zsw_chrombpnet/public/1.txt"

ATAC_BASE="/data/home/sczd644/run/zsw_chrombpnet/public/ATAC_bams"
MODEL_BASE="/data/home/sczd644/run/zsw_chrombpnet/chrombpnet_model"

REF_FA="/data/home/sczd644/run/zsw_chrombpnet/ref.fa"
CHROM_SIZES="/data/home/sczd644/run/zsw_chrombpnet/sheep.chrom.sizes"

OUT_BASE="/data/home/sczd644/run/zsw_chrombpnet/public/pred_bw"

while read tissue; do
  [ -z "$tissue" ] && continue

  tissue_dir="${ATAC_BASE}/${tissue}"
  tissue_use="$tissue"
  if [ ! -d "$tissue_dir" ]; then
    tissue_use=$(echo "$tissue" | tr '[:upper:]' '[:lower:]')
    tissue_dir="${ATAC_BASE}/${tissue_use}"
  fi

  model="${MODEL_BASE}/${tissue_use}_chrombpnet_model/models/chrombpnet_nobias.h5"
  if [ ! -s "$model" ]; then
    model="${MODEL_BASE}/${tissue}_chrombpnet_model/models/chrombpnet_nobias.h5"
  fi
  if [ ! -s "$model" ]; then
    echo "[WARN] model not found for ${tissue} (skip)"
    continue
  fi

  for data_dir in "${tissue_dir}"/Rep*/data; do
    [ -d "$data_dir" ] || continue
    pair_tag=$(basename "$(dirname "$data_dir")")   # Rep01-02

    bed="${data_dir}/peaks_no_blacklist.bed"
    [ -s "$bed" ] || { echo "[WARN] missing bed: $bed"; continue; }

    outdir="${OUT_BASE}/${tissue_use}/${pair_tag}"
    mkdir -p "$outdir"
    prefix="${outdir}/${tissue_use}_${pair_tag}_peaks"

    ls "${prefix}"*.bw >/dev/null 2>&1 && { echo "[SKIP] ${prefix}*.bw exists"; continue; }

    echo "[RUN] ${tissue_use} ${pair_tag}"
    chrombpnet pred_bw \
      -cmb "$model" \
      -r "$bed" \
      -g "$REF_FA" \
      -c "$CHROM_SIZES" \
      -op "$prefix"
  done

done < "$TISSUE_TXT"
