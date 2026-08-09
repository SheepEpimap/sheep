#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

Rscript 2.R \
  --tissues "abomasum" \
  --region_big  "chr1:264666438-264666657" \
  --region_small "chr1:264665438-264667657" \
  --big_tracks  "obs,pred" \
  --small_tracks "logo" \
  --highlight TRUE --highlight_alpha 0.25 \
  --auto_height TRUE \
  -o multi_tissue_dual.pdf \
  --base_pt 6 \
  --highlight TRUE \
  --gene_arrow_step_bp 20 \
  --hit_arrow_len 1 \
  --hit_arrow_step 6 \
  --yscale metric \
  --track_h_map "hits=1,genes=0.7"
