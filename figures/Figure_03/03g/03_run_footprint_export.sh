#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

IN_DIR="/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/result"
OUT_DIR="/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/result/result_pdf_all"
PY="/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint/result/export_footprints_onefolder_i0.py"

mkdir -p "$OUT_DIR"

shopt -s nullglob
for h5 in "${IN_DIR}"/*_footprints.h5; do
  base="$(basename "$h5")"              # abomasum_footprints.h5
  tissue="${base%_footprints.h5}"       # abomasum

  echo "[INFO] $tissue : $h5"
  python3 "$PY" \
    --in_h5 "$h5" \
    --tissue "$tissue" \
    --outdir "$OUT_DIR" \
    --plot_bp 200
done
