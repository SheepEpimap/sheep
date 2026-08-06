#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# ============================================================
#
# ${MATRIX_ROOT}/per_tissue_by_state/E1/abomasum.E1.ChromBPNet.mat.gz
# ${MATRIX_ROOT}/per_tissue_by_state/E2/abomasum.E2.ChromBPNet.mat.gz
# ...
#
# output:
# ============================================================

# ============================================================
# ============================================================

source /vol2/lixinyi/miniconda3_legacy/etc/profile.d/conda.sh
conda activate deeptools

# ============================================================
# ============================================================

TISSUE_LIST="/vol2/zhangshiwen/GWAS/GWAS_enrichment/chromatin/tissue.txt"

MATRIX_ROOT="/data/home/sczd644/run/zsw_chrombpnet/matrix/per_tissue_bed_chrombpnet_11states"

MATRIX_DIR="${MATRIX_ROOT}/per_tissue_by_state"

PROFILE_DIR="${MATRIX_ROOT}/profile_by_tissue_state"

FINAL_DIR="${MATRIX_ROOT}/final_profile_by_state"

mkdir -p "${PROFILE_DIR}" "${FINAL_DIR}"

# ============================================================
# ============================================================

UPSTREAM=10000
REGION_BODY=2000
DOWNSTREAM=10000
BIN_SIZE=50

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
# 4. check tissue list
# ============================================================

if [ ! -s "${TISSUE_LIST}" ]; then
    echo "[ERROR] tissue list not found or empty:"
    echo "${TISSUE_LIST}"
    exit 1
fi

echo "============================================================"
echo "[INFO] MATRIX_DIR  : ${MATRIX_DIR}"
echo "[INFO] PROFILE_DIR : ${PROFILE_DIR}"
echo "[INFO] FINAL_DIR   : ${FINAL_DIR}"
echo "============================================================"

# ============================================================
# ============================================================

for state in "${STATES[@]}"
do
    echo "============================================================"
    echo "[INFO] Extracting profiles for state=${state}"
    echo "============================================================"

    STATE_MATRIX_DIR="${MATRIX_DIR}/${state}"
    STATE_PROFILE_DIR="${PROFILE_DIR}/${state}"

    mkdir -p "${STATE_PROFILE_DIR}"

    while read tissue
    do
        [ -z "${tissue}" ] && continue
        [[ "${tissue}" =~ ^# ]] && continue
        [ "${tissue}" = "tissue" ] && continue

        MAT_FILE="${STATE_MATRIX_DIR}/${tissue}.${state}.ChromBPNet.mat.gz"
        PROFILE_TSV="${STATE_PROFILE_DIR}/${tissue}.${state}.ChromBPNet.profile.tsv"

        if [ ! -s "${MAT_FILE}" ]; then
            echo "[WARN] Missing matrix, skip: ${MAT_FILE}"
            continue
        fi

        if [ -s "${PROFILE_TSV}" ]; then
            echo "[INFO] Profile exists, skip: ${PROFILE_TSV}"
            continue
        fi

        echo "[INFO] Extracting profile: tissue=${tissue}, state=${state}"

        python - "${MAT_FILE}" "${PROFILE_TSV}" "${tissue}" "${state}" "${UPSTREAM}" "${REGION_BODY}" "${DOWNSTREAM}" <<'PY'
import sys
import numpy as np
from deeptools import heatmapper

mat_file = sys.argv[1]
out_tsv = sys.argv[2]
tissue = sys.argv[3]
state = sys.argv[4]

upstream = int(sys.argv[5])
region_body = int(sys.argv[6])
downstream = int(sys.argv[7])

hm = heatmapper.heatmapper()
hm.read_matrix_file(mat_file)

mat = hm.matrix.matrix

# rows = regions, columns = bins
profile = np.nanmean(mat, axis=0)

n_bins = len(profile)
total_len = upstream + region_body + downstream

relative_bp = np.linspace(
    -upstream,
    region_body + downstream,
    n_bins,
    endpoint=False
) + (total_len / n_bins) / 2.0

with open(out_tsv, "w") as out:
    out.write("tissue\tstate\tbin\trelative_bp\tregion_part\tsignal\n")

    for i, (x, y) in enumerate(zip(relative_bp, profile)):
        if x < 0:
            part = "upstream"
        elif x <= region_body:
            part = "body"
        else:
            part = "downstream"

        if np.isnan(y):
            y_str = "NA"
        else:
            y_str = f"{y:.8g}"

        out.write(
            f"{tissue}\t{state}\t{i}\t{x:.3f}\t{part}\t{y_str}\n"
        )
PY

        if [ ! -s "${PROFILE_TSV}" ]; then
            echo "[ERROR] Failed to generate profile:"
            echo "${PROFILE_TSV}"
            exit 1
        fi

        echo "[DONE] ${PROFILE_TSV}"

    done < "${TISSUE_LIST}"
done

# ============================================================
# ============================================================

for state in "${STATES[@]}"
do
    echo "============================================================"
    echo "[INFO] Merging profiles for state=${state}"
    echo "============================================================"

    STATE_PROFILE_DIR="${PROFILE_DIR}/${state}"

    OUT_LONG="${FINAL_DIR}/All_tissues.${state}.ChromBPNet.profile.long.tsv"
    OUT_MEAN="${FINAL_DIR}/All_tissues.${state}.ChromBPNet.profile.mean.tsv"

    rm -f "${OUT_LONG}" "${OUT_MEAN}"

    first=1

    while read tissue
    do
        [ -z "${tissue}" ] && continue
        [[ "${tissue}" =~ ^# ]] && continue
        [ "${tissue}" = "tissue" ] && continue

        PROFILE_TSV="${STATE_PROFILE_DIR}/${tissue}.${state}.ChromBPNet.profile.tsv"

        if [ ! -s "${PROFILE_TSV}" ]; then
            echo "[WARN] Missing profile, skip: ${PROFILE_TSV}"
            continue
        fi

        if [ "${first}" -eq 1 ]; then
            cat "${PROFILE_TSV}" > "${OUT_LONG}"
            first=0
        else
            tail -n +2 "${PROFILE_TSV}" >> "${OUT_LONG}"
        fi

    done < "${TISSUE_LIST}"

    if [ ! -s "${OUT_LONG}" ]; then
        echo "[WARN] No long profile generated for ${state}"
        continue
    fi

    # ========================================================
    # output:
    # state bin relative_bp region_part n_tissue mean_signal se_signal
    # ========================================================

    python - "${OUT_LONG}" "${OUT_MEAN}" <<'PY'
import sys
import numpy as np
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep="\t")
df["signal"] = pd.to_numeric(df["signal"], errors="coerce")

def se(x):
    x = x.dropna()
    if len(x) <= 1:
        return np.nan
    return x.std(ddof=1) / np.sqrt(len(x))

summary = (
    df.groupby(["state", "bin", "relative_bp", "region_part"], as_index=False)
      .agg(
          n_tissue=("signal", lambda x: x.dropna().shape[0]),
          mean_signal=("signal", "mean"),
          se_signal=("signal", se)
      )
)

summary.to_csv(outfile, sep="\t", index=False)
PY

    if [ ! -s "${OUT_MEAN}" ]; then
        echo "[ERROR] Failed to generate mean profile:"
        echo "${OUT_MEAN}"
        exit 1
    fi

    echo "[DONE] Long profile: ${OUT_LONG}"
    echo "[DONE] Mean profile: ${OUT_MEAN}"
done

# ============================================================
# ============================================================

ALL_MEAN="${FINAL_DIR}/All_tissues.E1_E11.ChromBPNet.profile.mean.tsv"

rm -f "${ALL_MEAN}"

first=1

for state in "${STATES[@]}"
do
    f="${FINAL_DIR}/All_tissues.${state}.ChromBPNet.profile.mean.tsv"

    if [ ! -s "${f}" ]; then
        echo "[WARN] Missing mean profile, skip: ${f}"
        continue
    fi

    if [ "${first}" -eq 1 ]; then
        cat "${f}" > "${ALL_MEAN}"
        first=0
    else
        tail -n +2 "${f}" >> "${ALL_MEAN}"
    fi
done

if [ -s "${ALL_MEAN}" ]; then
    echo "============================================================"
    echo "[DONE] Final E1-E11 mean profile:"
    echo "${ALL_MEAN}"
    echo "============================================================"
else
    echo "[ERROR] Final E1-E11 mean profile was not generated."
    exit 1
fi

echo "[INFO] All profile merging finished."
