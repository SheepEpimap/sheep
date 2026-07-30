#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
Rscript plot_bulk_tracks_ai_pdf.R \
 --id abomasum \
  --region_big chr1:11052498-11056712 \
  --region_small chr1:11054498-11054712 \
  --outpdf abomasum_track.pdf \
  --base_pt 6 \
  --highlight TRUE \
  --hit_arrow_len 1 \
  --hit_arrow_step 6


Rscript plot_bulk_tracks_ai_pdf.R \
  -i abomasum \
  --region_big chr1:11050598-11056812 \
  --region_small chr1:11050598-11056812 \
  --highlight FALSE \
  -o same_region_all_tracks.pdf


Rscript plot_bulk_tracks_ai_pdf.R \
  -i abomasum \
  --region_big chr1:11050598-11056812 \
  --region_small chr1:11054500-11054820 \
  --highlight TRUE \
  -o big_with_small_highlight.pdf


Rscript plot_bulk_tracks_ai_pdf.R \
 -i abomasum \
 --region_big chr1:11050598-11056812 \
 --pos_style bars \
 -o bars_style.pdf




  Rscript plot_bulk_tracks_ai_pdf.R \
  -i thyroid \
  --region_big chr9:21309559-21325598  \
  --region_small chr9:21310038-21310242  \
  -o big_with_small_highlight_thyroid.pdf \
  --base_pt 7 \
  --highlight TRUE \
  --hit_arrow_len 1 \
  --hit_arrow_step 3
