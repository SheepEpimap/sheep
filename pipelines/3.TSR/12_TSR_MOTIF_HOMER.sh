###########Final code#############
###########Use HOMER to identify motifs#############
###########Remove only sequences from the corresponding tissue from the background###########
##/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/Hormer_motif_remove1##
vim 03_TSR_remove1.sh
#!/bin/bash
set -e
id=$1
findMotifsGenome.pl ${id}_id.bed \
/vol2/mengzhu/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna \
Hormer_motif_remove1/${id}_remove1 \
-bg background_$(basename $id "_E5").bed \
-len 8,10,12 -size 200 -mask -p 5
# The script processes one ID supplied as $1; the unmatched source `done` was removed.
#for i in `cat sample.txt`; do sbatch -c 5  -p low --mem 10G 03_TSR_remove1.sh $i; done
