#!/bin/bash
for assay in H3K27ac H3K27me3 H3K4me1 H3K4me3 RNASeq
do

cd /vol2/mengzhu/snakemake_sheep/clean/bw

ma=$(grep ${assay} labels.txt|perl -p -e 's/\n/ /g')
echo ${ma}
multiBigwigSummary bins -p=5 \
-b ${assay}_*ZScores.bw \
-out ${assay}_MultiBigwigSummary.npz \
--labels ${ma} \
--binSize=1000

cd /vol2/mengzhu/snakemake_sheep
plotCorrelation -in ${assay}_MultiBigwigSummary.npz \
--corMethod pearson \
--skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap \
-o ${assay}_Pearson_Correlation.PDF \
--removeOutliers \
--colorMap bwr \
--outFileCorMatrix ${assay}_pearsonCorr_readCounts.tab


plotPCA -in ${assay}_MultiBigwigSummary.npz \
-o ${assay}_PCA_1_vs_2.pdf \
--plotTitle "PCA" \
--transpose --PCs 1 2 \
--outFileNameData ${assay}_PCA.tab
done


#!/bin/bash
for assay in ATAC
do
cd /vol2/mengzhu/snakemake_sheep/clean/bw

ma=$(grep ${assay} labels.txt|perl -p -e 's/\n/ /g')
echo ${ma}
multiBigwigSummary bins -p=24 \
-b  ATAC_*ZScores.bw \
-out ${assay}_MultiBigwigSummary.npz \
--labels ${ma} \
--binSize=1000 \
--outRawCounts ${assay}_MultiBigwigSummary.txt


plotCorrelation -in ${assay}_MultiBigwigSummary.npz \
--corMethod pearson \
--skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap \
-o ${assay}_Pearson_Correlation.PDF \
--removeOutliers \
--colorMap bwr \
--outFileCorMatrix ${assay}_pearsonCorr_readCounts.tab


plotPCA -in ${assay}_MultiBigwigSummary.npz \
-o ${assay}_PCA_1_vs_2.pdf \
--plotTitle "PCA" \
--transpose --PCs 1 2 \
--outFileNameData ${assay}_PCA.tab

done
