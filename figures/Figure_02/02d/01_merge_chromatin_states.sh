#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

mkdir /vol2/zhangshiwen/sheep_cor/
cd /vol2/zhangshiwen/sheep_cor
#Gs
####get Gs: total regions of each state (Gs) across multiple tissues

#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
mkdir AAGs
for i in  E5
do
echo $i
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done
