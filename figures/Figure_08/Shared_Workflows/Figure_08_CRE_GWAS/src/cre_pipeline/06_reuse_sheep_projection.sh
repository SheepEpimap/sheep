#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] line $LINENO, exit code $?" >&2' ERR

# The sheep-projected E5 sfCRE path depends on E5 sfCRE and sheep E5 reference.
# With the preserved overlap rule, E10/E11 blank changes sd/so evidence, not
# same-state E5 sfCRE. Validate the E5 sfCRE files, then reuse the projection.

: "${SRC_WORKDIR:?Set SRC_WORKDIR in the paths configuration}"
: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR in the paths configuration}"
USE_SYMLINKS="${USE_SYMLINKS:-1}"
FORCE_RELINK="${FORCE_RELINK:-0}"
VALIDATE_E5_SF="${VALIDATE_E5_SF:-1}"
REQUIRE_IDENTICAL_E5_SF="${REQUIRE_IDENTICAL_E5_SF:-1}"

HUMAN_TISSUES=(Adipose Colon Cortex Heart Liver Lung Muscle Ovary Sintest Spleen Stomach Testis)
SHEEP_TISSUES=(adipose colon cerebral-cortex heart liver lung muscle ovary jejunum spleen abomasum testis)

SRC_DIR="${SRC_WORKDIR}/sf_sheepE5_TSRlabeled"
DST_DIR="${WORKDIR}/sf_sheepE5_TSRlabeled"

log() { echo "[$(date +%F' '%T)] $*" >&2; }

[[ -d "${SRC_DIR}" ]] || { echo "ERROR: missing source projection directory: ${SRC_DIR}" >&2; exit 1; }
[[ -d "${SRC_DIR}/hg19" ]] || { echo "ERROR: missing source projection hg19 directory: ${SRC_DIR}/hg19" >&2; exit 1; }
mkdir -p "${WORKDIR}"

if [[ "${VALIDATE_E5_SF}" == "1" ]]; then
    log "validate E5 sfCRE files before reusing sheep projection"
    VALIDATION="${WORKDIR}/E5_sfCRE_reuse_validation.tsv"
    {
        echo -e "human_tissue\tsheep_tissue\told_lines\tnew_lines\tstatus"
    } > "${VALIDATION}"

    mismatches=0
    for i in "${!HUMAN_TISSUES[@]}"; do
        ht="${HUMAN_TISSUES[$i]}"
        st="${SHEEP_TISSUES[$i]}"
        old="${SRC_WORKDIR}/classified/${ht}_${st}_E5.sfCRE.bed"
        new="${WORKDIR}/classified/${ht}_${st}_E5.sfCRE.bed"
        [[ -e "${old}" ]] || { echo "ERROR: missing old sfCRE: ${old}" >&2; exit 1; }
        [[ -e "${new}" ]] || { echo "ERROR: missing new sfCRE: ${new}" >&2; exit 1; }

        old_n=$(wc -l < "${old}")
        new_n=$(wc -l < "${new}")
        if cmp -s "${old}" "${new}"; then
            status="identical"
        else
            status="different"
            mismatches=$((mismatches + 1))
        fi
        echo -e "${ht}\t${st}\t${old_n}\t${new_n}\t${status}" >> "${VALIDATION}"
    done

    log "validation table: ${VALIDATION}"
    if [[ "${mismatches}" -gt 0 && "${REQUIRE_IDENTICAL_E5_SF}" == "1" ]]; then
        echo "ERROR: ${mismatches} E5 sfCRE files differ; do not reuse sheep projection. Run the full projection step instead." >&2
        exit 1
    fi
fi

if [[ -e "${DST_DIR}" || -L "${DST_DIR}" ]]; then
    if [[ "${FORCE_RELINK}" == "1" ]]; then
        rm -rf "${DST_DIR}"
    else
        log "keep existing ${DST_DIR}"
        exit 0
    fi
fi

if [[ "${USE_SYMLINKS}" == "1" ]]; then
    ln -s "${SRC_DIR}" "${DST_DIR}"
    log "symlinked ${DST_DIR} -> ${SRC_DIR}"
else
    cp -a "${SRC_DIR}" "${DST_DIR}"
    log "copied ${SRC_DIR} -> ${DST_DIR}"
fi
