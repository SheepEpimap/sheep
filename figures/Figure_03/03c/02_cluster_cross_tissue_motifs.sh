#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
module load cuda/11.8
module load cudnn/8.9.6.50_cuda11

set -euo pipefail

echo "=================================================="
echo " motif "
echo " : $(date)"
echo "=================================================="



INPUT_PFM="/data/home/sczd644/run/zsw_chrombpnet/cluster/all_tissues_motifs.pfm"
INPUT_METADATA="/data/home/sczd644/run/zsw_chrombpnet/cluster/all_tissues_motifs_metadata.tsv"
CLUSTER_THRESHOLD=0.9 #
OUTPUT_DIR="/data/home/sczd644/run/zsw_chrombpnet/cluster/cross_tissue_cluster_analysis"

mkdir -p $OUTPUT_DIR

echo " :"
echo "  - inputPFM: $INPUT_PFM"
echo "  - input : $INPUT_METADATA"
echo "  -  : $CLUSTER_THRESHOLD"
echo "  - Output directory: $OUTPUT_DIR"
echo ""

# checkinputfile
if [[ ! -f $INPUT_PFM ]]; then
    echo " : Not foundPFMfile $INPUT_PFM"
    exit 1
fi

if [[ ! -f $INPUT_METADATA ]]; then
    echo " : Not found file $INPUT_METADATA"
    exit 1
fi

MOTIF_COUNT=$(grep -c '>' $INPUT_PFM)
TISSUE_COUNT=$(awk -F'\t' 'NR>1 {print $2}' $INPUT_METADATA | sort -u | wc -l)

echo "input statistics:"
echo "  - Motif : $MOTIF_COUNT"
echo "  -  : $TISSUE_COUNT"
echo ""

echo " 1:  motif ..."
CLUSTER_RESULTS="$OUTPUT_DIR/cluster_results"
mkdir -p $CLUSTER_RESULTS

gimme cluster $INPUT_PFM $CLUSTER_RESULTS -t $CLUSTER_THRESHOLD

if [[ $? -eq 0 ]]; then
    echo "✓  "
else
    echo "✗  "
    exit 1
fi
