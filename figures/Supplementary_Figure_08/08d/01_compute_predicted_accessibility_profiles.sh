#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# ============================================================
# bash compute_chrombpnet_one_tissue_per_tissue_bed.sh abomasum
# bash compute_chrombpnet_one_tissue_per_tissue_bed.sh hippocampus
# ============================================================

tissue=${1}

# ============================================================
# ============================================================


# ============================================================
# ============================================================

# ChromBPNet predicted accessibility bigWig
BW_ROOT="/data/home/sczd644/run/zsw_chrombpnet/pred_bw"

# hippocampus_E6.bed
# abomasum_E1.bed
BED_DIR="/vol2/zhangshiwen/GWAS/GWAS_enrichment/chromatin"

BW_FILE="${BW_ROOT}/${tissue}/${tissue}_peaks_chrombpnet_nobias.bw"

# ============================================================
# ============================================================

OUT_ROOT="/data/home/sczd644/run/zsw_chrombpnet/matrix/per_tissue_bed_chrombpnet_11states"

PER_TISSUE_DIR="${OUT_ROOT}/per_tissue_by_state"
LOG_DIR="${OUT_ROOT}/logs"

mkdir -p "${PER_TISSUE_DIR}" "${LOG_DIR}"

# ============================================================
# ============================================================

UPSTREAM=10000
REGION_BODY=2000
DOWNSTREAM=10000
BIN_SIZE=50
THREADS=8

STATES=(
E1
E2
E3
E4
E5
E6
E7
E8
E9
E10
E11
)

# ============================================================
# 5. check bigWig
# ============================================================

if [ ! -s "${BW_FILE}" ]; then
    echo "[ERROR] ChromBPNet bigWig not found or empty:"
    echo "${BW_FILE}"
    exit 1
fi

echo "============================================================"
echo "[INFO] Tissue : ${tissue}"
echo "[INFO] BW     : ${BW_FILE}"
echo "[INFO] BEDDIR : ${BED_DIR}"
echo "[INFO] OUT    : ${OUT_ROOT}"
echo "============================================================"

# ============================================================
# ============================================================
#
#
# ============================================================

for state in "${STATES[@]}"
do
    BED_FILE="${BED_DIR}/${tissue}_${state}.bed"

    STATE_OUT_DIR="${PER_TISSUE_DIR}/${state}"
    mkdir -p "${STATE_OUT_DIR}"

    OUT_PREFIX="${STATE_OUT_DIR}/${tissue}.${state}.ChromBPNet"

    echo "------------------------------------------------------------"
    echo "[INFO] Processing tissue=${tissue}, state=${state}"
    echo "[INFO] BED: ${BED_FILE}"
    echo "[INFO] BW : ${BW_FILE}"
    echo "[INFO] OUT: ${OUT_PREFIX}.mat.gz"
    echo "------------------------------------------------------------"

    if [ ! -s "${BED_FILE}" ]; then
        echo "[WARN] BED file not found or empty, skip:"
        echo "${BED_FILE}"
        continue
    fi

    computeMatrix scale-regions \
        -S "${BW_FILE}" \
        -R "${BED_FILE}" \
        --beforeRegionStartLength "${UPSTREAM}" \
        --regionBodyLength "${REGION_BODY}" \
        --afterRegionStartLength "${DOWNSTREAM}" \
        --binSize "${BIN_SIZE}" \
        --missingDataAsZero \
        --samplesLabel "${tissue}" \
        -p "${THREADS}" \
        -o "${OUT_PREFIX}.mat.gz" \
        --outFileNameMatrix "${OUT_PREFIX}.matrix.tsv" \
        --outFileSortedRegions "${OUT_PREFIX}.sorted.bed"

    if [ ! -s "${OUT_PREFIX}.mat.gz" ]; then
        echo "[ERROR] computeMatrix failed:"
        echo "${OUT_PREFIX}.mat.gz"
        exit 1
    fi

    echo "[DONE] ${OUT_PREFIX}.mat.gz"
done

echo "============================================================"
echo "[DONE] Finished ChromBPNet matrices for tissue: ${tissue}"
echo "============================================================"
