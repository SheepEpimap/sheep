#!/bin/bash
cd /vol2/mengzhu/snakemake_sheep/clean/bw
rm  labels.txt
ls ATAC*ZScores.bw H3*ZScores.bw RNASeq*ZScores.bw | while read id;
do
   echo $(basename $id "_ZScores.bw") >> labels.txt
done

#!/bin/bash
# Read the label file and combine its contents into one line
ma=$(cat labels.txt | perl -p -e 's/\n/ /g')

# Check the number of matching BigWig files
bw_files=$(ls ATAC_*ZScores.bw H3K27ac_*ZScores.bw H3K27me3*ZScores.bw H3K4me1*ZScores.bw H3K4me3*ZScores.bw RNASeq*ZScores.bw)
num_bw_files=$(echo "$bw_files" | wc -w)
num_labels=$(echo "$ma" | wc -w)

echo "Number of BigWig files: $num_bw_files"
echo "Number of labels: $num_labels"

# Check whether the numbers of labels and files match
if [ "$num_bw_files" -ne "$num_labels" ]; then
    echo "Error: The number of labels does not match the number of BigWig files."
    exit 1
fi

# Print matching files and labels for debugging
echo "BigWig files: $bw_files"
echo "Labels: $ma"

# Run multiBigwigSummary
multiBigwigSummary bins -p=110 \
    -b $bw_files \
    -out all_MultiBigwigSummary.npz \
    --labels $ma \
    --binSize=1000 \
    --outRawCounts all_MultiBigwigSummary.txt

plotCorrelation -in /vol2/mengzhu/snakemake_sheep/clean/bw/all_MultiBigwigSummary.npz \
--corMethod pearson \
--skipZeros --plotTitle "Pearson Correlation of all Read Counts" \
--whatToPlot heatmap \
-o /vol2/mengzhu/snakemake_sheep/clean/bw/all_Pearson_Correlation.PDF \
--removeOutliers \
--colorMap bwr \
--outFileCorMatrix /vol2/mengzhu/snakemake_sheep/clean/bw/pearsonCorr_readCounts.tab

sed "s/'//g" pearsonCorr_readCounts.tab > 1.txt
cut -f 1 1.txt | cut -d "_" -f 1  > 2.txt
cut -f 1 1.txt | cut -d "_" -f 2 > 3.txt
cut -f 1 1.txt | cut -d "_" -f 3 > 4.txt
paste 1.txt 2.txt 3.txt 4.txt > pearsonCorr_readCounts_last.csv
sed -i '1d' pearsonCorr_readCounts_last.csv
cut -f 1 pearsonCorr_readCounts.tab | cut -d "_" -f 1 | sed "s/'//g"    #Remove single quotes

#!/bin/bash
cd /vol2/mengzhu/snakemake_sheep/clean/bw
sed '1d' pearsonCorr_readCounts.tab | sed "s/'//g" > 1.txt
cut -f 1 1.txt | cut -d "_" -f 1  > 2.txt
cut -f 1 1.txt | cut -d "_" -f 2 > 3.txt
cut -f 1 1.txt | cut -d "_" -f 3 > 4.txt

cd /vol2/mengzhu/snakemake_sheep/clean/bw
cp 3.txt 3_1.txt
cat /vol2/mengzhu/snakemake_sheep/Figures/tissue_list_sheep.txt |sed 's/\t/_/g' > 6.txt
paste /vol2/mengzhu/snakemake_sheep/Figures/tissue_list_sheep.txt 6.txt >tissue_list_1.txt 
sed 's/\r/\t/g' tissue_list_1.txt > tissue_list_1_modified.txt
cat tissue_list_1_modified.txt | while read line
do 
echo $line
arr=($line)
tissue=${arr[0]}
layer1=${arr[2]}
layer2=${arr[3]}
echo $tissue
sed -i "s/\<$tissue\>/$layer2/g" 3_1.txt 
done
cat 3_1.txt |sed 's/_/\t/g' > 3_2.txt
sed '1c Tissues Layer Bigtype' 3_2.txt |  sed 's/ /\t/g' > 3_3.txt
paste 1.txt 2.txt 3_3.txt 4.txt > PearsonCorr_readCounts_last.csv


