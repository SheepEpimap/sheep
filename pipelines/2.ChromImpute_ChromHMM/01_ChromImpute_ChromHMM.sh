#!/bin/bash
#This step follows peak calling and produces *_treat_pileup.bdg and *_control_lambda.bdg files for each sample
#ls *_treat_pileup.bdg | cut -d '_' -f 1-3 |while read id; do
id=$1
echo $id
chipReads=$(cat /vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/${id}.bed | wc -l | awk '{printf "%f", $1/1000000}')
macs2 bdgcmp -t ${id}_treat_pileup.bdg -c ${id}_control_lambda.bdg -o ${id}_ppois.bdg -m ppois -S ${chipReads}
slopBed -i ${id}_ppois.bdg -g /vol2/mengzhu/genome/reference/sheep.size -b 0 | bedClip stdin /vol2/mengzhu/genome/reference/sheep.size ${id}.pval.signal.bedgraph
awk 'BEGIN{OFS="\t"}
    NR==1 {
      last_chrom=$1; last_end=$3
      print $0
      next
    }
    {
      if ($1 == last_chrom && $2 < last_end) {
        $2 = last_end
      }
      if ($2 < $3) {
        print $0
        last_chrom = $1; last_end = $3
      }
    }' ${id}.pval.signal.bedgraph > ${id}.pval.signal_1.bedgraph
sort -k1,1 -k2,2n ${id}.pval.signal_1.bedgraph > ${id}.pval.signal.bedgraph.tmp
bedGraphToBigWig ${id}.pval.signal.bedgraph.tmp /vol2/mengzhu/genome/reference/sheep.size ${id}.pval.signal.bigwig
gzip -c -f ${id}.pval.signal.bedgraph.tmp > ${id}.pval.signal.bedgraph.gz
#done
#sbatch -p low -c 8 --mem=32G -t 10-0 01_pval_signal_bedgraph.sh
#for i in cat sample.tab; do sbatch -p low -w comput3 -c 4 --mem=8G -t 10-0 01_pval_signal_bedgraph.sh $i ; done

#0.Prepare Imputed_samples.tab
mkdir CHROMHMMDIR CONVERTEDDIR DISTANCEDIR IMPUTED INPUTDATADIR PREDICTORDIR TRAINDATA
for bdg in *.pval.signal.bedgraph.gz; do
  # Extract the tissue name
  tissue=$(echo "$bdg" | cut -d'_' -f 2-3 | cut -d '.' -f 1)
  # Extract the mark name (the prefix, such as ATAC or H3K27ac)
  mark=$(echo "$bdg" | cut -d'_' -f1)
  # Get the full path
  echo -e "${tissue}\t${mark}\t${bdg}" >> Imputed_samples.tab
done


#1.Convert data to the format required by ChromImpute
vim ChromImpute_Convert1.sh 
#!/bin/bash
java -jar -Xmx40G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar Convert INPUTDATADIR  Imputed_samples.tab /vol2/mengzhu/genome/reference/sheep1.size CONVERTEDDIR


#2.Generate the training dataset
vim ChromImpute_ComputeGlobalDist.sh
#!/bin/bash
java -jar -Xmx32G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ComputeGlobalDist -m ${mark} /vol2/mengzhu/soft/ChromImpute/sheep/CONVERTEDDIR ${sample_tab} /vol2/mengzhu/genome/reference/sheep.size /vol2/mengzhu/soft/ChromImpute/sheep/DISTANCEDIR
java -jar -Xmx32G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar GenerateTrainData /vol2/mengzhu/soft/ChromImpute/sheep/CONVERTEDDIR /vol2/mengzhu/soft/ChromImpute/sheep/DISTANCEDIR ${sample_tab} /vol2/mengzhu/genome/reference/sheep.size /vol2/mengzhu/soft/ChromImpute/sheep/TRAINDATA ${mark}
for i in H3K4me1 H3K4me3 H3K27ac H3K27me3 ATAC RNASeq; do sbatch -p smp -c 4 --mem=32G -t 10-0 ChromImpute_ComputeGlobalDist.sh $i imputed_samples.tab; done

#3.Start training and generate imputed results (cannot run in parallel)
vim ChromImpute_TrainApply.sh
#!/bin/bash
set -e
module load openjdk/16.0.2
EID=$1
assay=$2
sample_tab=$3
java -jar -Xmx8G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar Train TRAINDATA ${sample_tab} PREDICTORDIR ${EID} ${assay}
java -jar -Xmx8G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar Apply CONVERTEDDIR DISTANCEDIR PREDICTORDIR ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size IMPUTED ${EID} ${assay}
#for i in `cat Imputed_samples_39.tab | awk '{print $1}' | sort | uniq`; do for j in H3K4me3 H3K27ac H3K27me3; do sbatch -p smp -c 4 --mem=8G -t 10-0 -o Logs/GeneratePredictors.${i}_${j}.%j.out ChromImpute_TrainApply.sh $i $j Imputed_samples.tab; done; done
#for i in `cat Imputed_samples_40.tab | awk '{print $1}' | sort | uniq`; do for j in ATAC H3K4me1 H3K4me3 H3K27ac H3K27me3; do sbatch -p low -w comput3 -c 4 --mem=8G -t 10-0 -o Logs/GeneratePredictors.${i}_${j}.%j.out ChromImpute_TrainApply.sh $i $j Imputed_samples_40.tab; done; done

#4.Convert predictions to the input format required for ChromHMM learning
#The -g 2 option enables binarization; the exported files are already binarized and can be passed directly to ChromHMM LearnModel
mkdir ATAC H3K27ac H3K27me3 H3K4me1 H3K4me3
vim 04_ExportToChromHMM_ATAC.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/ATAC
#for i in {1..23}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_ATAC.${i}.out 04_ExportToChromHMM_ATAC.sh Temp/export2chromhmm_ATAC.${i}.tab; done
vim 04_ExportToChromHMM_H3K4me1.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K4me1
#for i in {1..3}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K4me1.${i}.out 04_ExportToChromHMM_H3K4me1.sh Temp/export2chromhmm_H3K4me1.${i}.tab; done
vim 04_ExportToChromHMM_H3K4me3.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K4me3
#for i in {1}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K4me3.1.out 04_ExportToChromHMM_H3K4me3.sh Temp/export2chromhmm_H3K4me3.1.tab; done
vim 04_ExportToChromHMM_H3K27ac.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K27ac
#for i in {1..20}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K27ac.${i}.out 04_ExportToChromHMM_H3K27ac.sh Temp/export2chromhmm_H3K27ac.${i}.tab; done
vim 04_ExportToChromHMM_H3K27me3.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K27me3
#for i in {1..26}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K27me3.${i}.out 04_ExportToChromHMM_H3K27me3.sh Temp/export2chromhmm_H3K27me3.${i}.tab; done


#5.Merge imputed data into the binarized files generated by ChromHMM
5.1 ATAC
#!/usr/bin/env bash
set -euo pipefail

ATAC_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/ATAC"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${ATAC_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    ATAC_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $ATAC_TXT"

    awk 'NR>2{print $1}' "$ATAC_TXT" > atac_col.tmp
    tail -n +3 "$MATRIX" | cut -f2- > other_marks.tmp
    paste atac_col.tmp other_marks.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm atac_col.tmp other_marks.tmp body.tmp
done

5.2 H3K27ac
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K27ac"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${H3_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    awk 'NR>2{print $1}' "$H3_TXT" > h27ac_col.tmp
    tail -n +3 "$MATRIX" > matrix_body.tmp

    cut -f1 matrix_body.tmp > col1.tmp
    cut -f3- matrix_body.tmp > others.tmp

    paste col1.tmp h27ac_col.tmp others.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h27ac_col.tmp matrix_body.tmp col1.tmp others.tmp body.tmp
done

5.3 H3K27me3
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K27me3"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${H3_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    # Extract the 0/1 column starting at row 3 from the single-column H3K27me3 file
    awk 'NR>2{print $1}' "$H3_TXT" > h27me3_col.tmp

    # Matrix body starting at row 3
    tail -n +3 "$MATRIX" > matrix_body.tmp

    # Retain columns 1-2 of the original matrix (ATAC, H3K27ac)
    cut -f1-2 matrix_body.tmp > col12.tmp
    # Retain column 4 onward from the original matrix (H3K4me1, H3K4me3)
    cut -f4-  matrix_body.tmp > others.tmp

    paste col12.tmp h27me3_col.tmp others.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h27me3_col.tmp matrix_body.tmp col12.tmp others.tmp body.tmp
done

5.4 H3K4me1
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K4me1"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${H3_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    # Extract the 0/1 column starting at row 3 from the single-column H3K4me1 file
    awk 'NR>2{print $1}' "$H3_TXT" > h4me1_col.tmp

    # Matrix body starting at row 3
    tail -n +3 "$MATRIX" > matrix_body.tmp

    # Retain columns 1-3 of the original matrix (ATAC, H3K27ac, H3K27me3)
    cut -f1-3 matrix_body.tmp > col123.tmp
    # Retain column 5 onward from the original matrix (H3K4me3 etc.)
    cut -f5-  matrix_body.tmp > others.tmp

    # Reassemble: ATAC, H3K27ac, H3K27me3, new H3K4me1, H3K4me3...
    paste col123.tmp h4me1_col.tmp others.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h4me1_col.tmp matrix_body.tmp col123.tmp others.tmp body.tmp
done

5.5 H3K4me3
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K4me3"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_39"

for f in "${H3_DIR}"/*39*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    awk 'NR>2{print $1}' "$H3_TXT" > h4me3_col.tmp
    tail -n +3 "$MATRIX" > matrix_body.tmp

    cut -f1-4 matrix_body.tmp > col1234.tmp

    paste col1234.tmp h4me3_col.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h4me3_col.tmp matrix_body.tmp col1234.tmp body.tmp
done



#5.Use ChromHMM to predict chromatin states
vim ChromHMM_ModelOptm.sh
#!/bin/bash
# Script: ChromHMM_ModelOptm.sh
mkdir -p ChroHMM/LearnModel_${num_model}
num_model=$1
java -jar -Xmx64G /vol2/mengzhu/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /vol2/mengzhu/genome/reference/sheep1.size /vol2/mengzhu/soft/ChromImpute/sheep/CHROMHMMDIR ChroHMM/LearnModel_${num_model} ${num_model} Ramb2
# for i in {10..18}; do sbatch -p low -c 12 --mem=48G -t 10-0 -o LearnModel_${i}.%j.out ChromHMM_ModelOptm.sh $i; done

#6.Generate 2-20 chromatin-state solutions for Rep1 and Rep2
vim 02_LearnModel.sh
#!/bin/bash
num_model=$1
mkdir Rep1/LearnModel_${num_model}
mkdir Rep2/LearnModel_${num_model}
java -jar -Xmx22G /vol2/mengzhu/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /vol2/mengzhu/genome/reference/sheep1.size SAMPLEDATA_sheep_five_maker_39 Rep1/LearnModel_${num_model} ${num_model} Ramb2
java -jar -Xmx22G /vol2/mengzhu/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /vol2/mengzhu/genome/reference/sheep1.size SAMPLEDATA_sheep_five_maker_40 Rep2/LearnModel_${num_model} ${num_model} Ramb2
# for i in {2..20}; do sbatch -p low -c 11 --mem=22G -t 10-0 -o Logs/LearnModel_${i}.out 02_CompareModels.sh ${i}; done


#Plot correlations among the 2-20 chromatin-state solutions and select the optimal number of states
mkdir emissions_Rep1
cp LearnModel_*/emissions_*.txt emissions_Rep1
mkdir emissions_Rep2
cp LearnModel_*/emissions_*.txt emissions_Rep2
vim 03_CompareModels_Rep2.sh
#!/bin/bash
#four maker and 16 chromatin state 39
java -mx80000M -jar /vol2/mengzhu/ChromHMM/ChromHMM.jar CompareModels \
    -color 255,0,0 \
    /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/LearnModel_20/emissions_20.txt \
    /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/emissions_Rep2 \
    /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/CompareModels/Rep2
#sbatch -c 40  -p low --mem 80G 03_CompareModels_Rep2.sh
# Extract the first column, remove rows with the value "1", transpose, and save as the first row
cut -f1 Rep1.txt | grep -v "^1$" | tr '\n' '\t' > Rep1_average.txt
echo "" >> Rep1_average.txt  # Append a newline

# Calculate the mean of each column and insert "average" before the first column of the second row
awk '
NR > 1 {
    for(i=2; i<=NF; i++) {
        sum[i] += $i;
    }
    count++;
}
END {
    printf "average\t";  # Insert "average" before the first column of the second output row
    for(i=2; i<=NF; i++) {
        printf "%.2f\t", sum[i]/count;
    }
    print "";  # Newline
}' Rep1.txt >> Rep1_average.txt
