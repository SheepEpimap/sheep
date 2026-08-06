#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] line $LINENO, exit code $?" >&2' ERR

# Human E5 TSR labels depend only on human E5 hg38 intervals, not on the
# E10/E11 blank rule. Reuse the final existing TSR directory.

: "${SRC_WORKDIR:?Set SRC_WORKDIR in the paths configuration}"
: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR in the paths configuration}"
STATE="${STATE:-E5}"
USE_SYMLINKS="${USE_SYMLINKS:-1}"
FORCE_RELINK="${FORCE_RELINK:-0}"

SRC_DIR="${SRC_WORKDIR}/human_TSR_${STATE}_hg38"
DST_DIR="${WORKDIR}/human_TSR_${STATE}_hg38"

log() { echo "[$(date +%F' '%T)] $*" >&2; }

[[ -s "${SRC_DIR}/labeled.bed" ]] || { echo "ERROR: missing source TSR labels: ${SRC_DIR}/labeled.bed" >&2; exit 1; }
mkdir -p "${WORKDIR}"

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
