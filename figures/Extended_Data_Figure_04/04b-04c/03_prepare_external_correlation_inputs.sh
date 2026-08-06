#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

TISSUE_TXT="/data/home/sczd644/run/zsw_chrombpnet/public/tissue.txt"

ATAC_BASE="/data/home/sczd644/run/zsw_chrombpnet/public/ATAC_bams"
PRED_BASE="/data/home/sczd644/run/zsw_chrombpnet/public/pred_bw"

while read -r tissue; do
  [ -z "${tissue}" ] && continue

  tuse="${tissue}"
  tdir="${ATAC_BASE}/${tuse}"
  if [ ! -d "$tdir" ]; then
    tuse=$(echo "$tissue" | tr '[:upper:]' '[:lower:]')
    tdir="${ATAC_BASE}/${tuse}"
  fi
  [ -d "$tdir" ] || { echo "[WARN] no ATAC dir for ${tissue}"; continue; }

  for data_dir in "${tdir}"/Rep*/data; do
    [ -d "$data_dir" ] || continue

    pair_tag=$(basename "$(dirname "$data_dir")")   # Rep01-02
    bed="${data_dir}/peaks_no_blacklist.observed.bed"
    obs_bw="${data_dir}/merged.bw"
    obs_out="${data_dir}/peaks_no_blacklist.observed.out"

    # ---------- observed ----------
    if [ -s "$obs_bw" ] && [ -s "$bed" ]; then
      if [ ! -s "$obs_out" ]; then
        echo "[RUN] observed  ${tuse} ${pair_tag}"
        bigWigAverageOverBed "$obs_bw" "$bed" "$obs_out"
      else
        echo "[SKIP] observed out exists: ${obs_out}"
      fi
    else
      echo "[WARN] missing observed bw/bed: ${tuse} ${pair_tag}"
      continue
    fi

    # ---------- predicted ----------
    pred_dir="${PRED_BASE}/${tuse}/${pair_tag}"
    pred_bw=""

    try1="${pred_dir}/${tuse}_${pair_tag}_peaks_chrombpnet_nobias.bw"
    [ -s "$try1" ] && pred_bw="$try1"

    if [ -z "$pred_bw" ]; then
      pred_bw=$(ls -1 "${pred_dir}"/*chrombpnet_nobias*.bw 2>/dev/null | head -n 1 || true)
    fi

    if [ -z "$pred_bw" ]; then
      pred_bw=$(ls -1 "${pred_dir}"/*.bw 2>/dev/null | head -n 1 || true)
    fi

    pred_out="${pred_dir}/${tuse}_${pair_tag}_peaks_no_blacklist.predicted.out"

    if [ -n "$pred_bw" ] && [ -s "$pred_bw" ]; then
      mkdir -p "$pred_dir"
      if [ ! -s "$pred_out" ]; then
        echo "[RUN] predicted ${tuse} ${pair_tag}"
        bigWigAverageOverBed "$pred_bw" "$bed" "$pred_out"
      else
        echo "[SKIP] predicted out exists: ${pred_out}"
      fi
    else
      echo "[WARN] predicted bw not found: ${tuse} ${pair_tag} (dir=${pred_dir})"
    fi

  done
done < "$TISSUE_TXT"
