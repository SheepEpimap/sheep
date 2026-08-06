#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

annoFile="/data/home/sczd644/run/zsw_chrombpnet/uniquemotif_result/simplified_motif_anno.tsv"   # pattern<TAB>TF
tissueList="/data/home/sczd644/run/zsw_chrombpnet/tissue.txt"
baseDir="/data/home/sczd644/run/zsw_chrombpnet/finemo"
outDir="/data/home/sczd644/run/zsw_chrombpnet/03-syntax/04/footprint"
pwmFile="$outDir/motif_pwm.tsv"   # TF<TAB>

awk 'NR>1{print $1"\t"$2}' "$annoFile" > pattern2tf.tmp

while read -r id; do
    echo "[$id]   pattern→TF ..."
    inBed="$baseDir/${id}_finemo/hits.bed"
    outBed="$baseDir/${id}_finemo/hits_tf.bed"
    uniqTF="$baseDir/${id}_finemo/motif_unique_tf.tsv"

    [[ -f "$inBed" ]] || { echo "WARN: $inBed  , "; continue; }

    awk -F'\t' -v OFS='\t' '
        BEGIN{while((getline<"pattern2tf.tmp")>0){pat[$1]=$2}}
        {
            split($4,a,".");
            pattern=a[2];
            if(pattern in pat) $4=pat[pattern];
            print
        }
    ' "$inBed" > "$outBed"

    cut -f4 "$outBed" | sort -u | awk -v org="$id" '{print org"\t"$0}' > "$uniqTF"

done < "$tissueList"

cat "$baseDir"/*_finemo/motif_unique_tf.tsv > "$outDir/all_motif_unique_tf.tsv"

awk 'NR>1{print $1"\t"$2}' "$pwmFile" > tf2seq.tmp
awk -F'\t' 'NR==FNR{seq[$1]=$2;next}
            {print $0 "\t" ($2 in seq ? seq[$2] : "NA")}' \
    tf2seq.tmp "$outDir/all_motif_unique_tf.tsv" \
    > "$outDir/all_motif_unique_tf_withSeq.tsv"

while read -r id; do
    awk -F'\t' -v org="$id" '$1==org{print $2"\t"$3}' \
        "$outDir/all_motif_unique_tf_withSeq.tsv" \
        > "$outDir/${id}_motif_pwm.tsv"
done < "$tissueList"

rm -f pattern2tf.tmp tf2seq.tmp
echo " !hits_tf.bed  , file  $outDir"

mkdir -p /data/home/sczd644/run/zsw_chrombpnet/footprint/result
#!/bin/bash
module load cuda/11.8
module load cudnn/8.9.6.50_cuda11

 cat /data/home/sczd644/run/zsw_chrombpnet/1.txt  |while read id
do

chrombpnet footprints -m /data/home/sczd644/run/zsw_chrombpnet/chrombpnet_model/${id}_chrombpnet_model/models/chrombpnet_nobias.h5 -r /data/home/sczd644/run/zsw_chrombpnet/nopeak/${id}_negatives.bed -g /data/home/sczd644/run/zsw_chrombpnet/ref.fa -fl /data/home/sczd644/run/zsw_chrombpnet/fold_0.json -op /data/home/sczd644/run/zsw_chrombpnet/footprint/result/${id} -pwm_f  /data/home/sczd644/run/zsw_chrombpnet/footprint/${id}_motif_pwm.tsv
done
