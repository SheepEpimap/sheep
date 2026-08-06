vim tss.sh
#!/bin/bash
cd /vol2/mengzhu/snakemake_sheep
for i in protein_coding lncRNA pseudogene miRNA snRNA snoRNA 
do
echo $i
grep $i /vol2/mengzhu/genome/Gene_esemble100_colin.bed > 2.bed
    for tissue in cerebellum_39;
    do
    computeMatrix scale-regions \
    -R 2.bed \
    -S "DeepTools/H3K27ac_"$tissue"_ZScores.bw" "DeepTools/H3K27me3_"$tissue"_ZScores.bw" "DeepTools/H3K4me1_"$tissue"_ZScores.bw" "DeepTools/H3K4me3_"$tissue"_ZScores.bw" "DeepTools/ATAC_"$tissue"_ZScores.bw" \
    -a 2500 -b 2500 \
    -out ComputeMatrix/${tissue}_${i}_region.mat.gz \
    -p 24 \
    --skipZeros

    plotProfile -m ComputeMatrix/${tissue}_${i}_region.mat.gz \
    -out Figures/${tissue}_${i}_region.PDF \
    --perGroup \
    --colors red gray yellow green pink \
    --samplesLabel H3K27ac H3K27me3 H3K4me1 H3K4me3 ATAC \
    -z ""
    done
done
