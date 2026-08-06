
fasta-get-markov -m 1 TSR_adipose_E5_enhancer.fa TSR_adipose_E5_enhancer.bg
fimo --bgfile enhancers.bg --qv-thresh --thresh 0.01 motifs.meme enhancers.fa


#Generate sample.txt
for f in *_E5_id.bed; do
  echo "${f%_id.bed}"
done > sample.txt

########Use MEME FIMO to identify TFBSs##########
##/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5##
vim 02_fimo_motif_p0.001.sh
#!/bin/bash
set -euo pipefail

id="$1"
#fasta-get-markov -m 1 motif/${id}_enhancer.fa motif/${id}_enhancer.bg
motif_file="/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E4/motif/jaspar_sheep_core_meme/JASPAR2024_CORE_vertebrates_redundant_pfms_meme.txt"
seq_file="motif/${id}_enhancer.fa"
bg_file="motif/${id}_enhancer.bg"
outdir="motif/fimo_vertebrates_bg0.001/${id}_enhancer"

mkdir -p "$(dirname "$outdir")"

# To use q-value (FDR) as the threshold (recommended for enhancer scans), uncomment the following line:
# Common q-value thresholds are 0.01 or 0.05; use 0.001 for greater stringency
fimo --oc "$outdir" \
     --max-stored-scores 5000000 \
     --thresh 0.001 \
     --bfile "$bg_file" \
     "$motif_file" \
     "$seq_file"
#for i in `cat sample1.txt`; do sbatch -c 40  -p low --mem 80G 02_fimo_motif_p0.001.sh $i; done
sbatch -c 5  -p smp --mem 40G 02_fimo_motif_p0.001.sh TSR_adipose_E5z
