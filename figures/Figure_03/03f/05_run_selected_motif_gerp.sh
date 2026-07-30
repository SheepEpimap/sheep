#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

python plot_pick_motif_gerp.py \
  --results_dir /data/home/sczd644/run/zsw_chrombpnet/phylop_gerp_all/results \
  --motif "CTCF#2" \
  --window 25 \
  --box_width 30 \
  --out /data/home/sczd644/run/zsw_chrombpnet/phylop_gerp_all/plots/CTCF_2.w25.box30
