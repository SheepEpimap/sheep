#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

cd /vol2/zhangshiwen/sheep_cor/
for state in E5 E6 E7 E8
do

input_bed="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs/${state}_Gs.bed"
blacklist_bed="/vol2/zhangshiwen/blacklist/output/sheep_blacklist.bed"
output_bed="/vol2/zhangshiwen/sheep_cor/${state}_no_blacklist_Gs_new.bed"

# sort
echo "Sorting blacklist file..."
sort -k1,1 -k2,2n "$blacklist_bed" >/vol2/zhangshiwen/blacklist/output/sheep_sorted_blacklist.bed

# remove
echo "Removing blacklist regions from $input_bed..."
bedtools intersect -v -a <(sort -k1,1 -k2,2n "$input_bed") -b /vol2/zhangshiwen/blacklist/output/sheep_sorted_blacklist.bed > "$output_bed"

echo "Process completed. Output saved to $output_bed"

done
