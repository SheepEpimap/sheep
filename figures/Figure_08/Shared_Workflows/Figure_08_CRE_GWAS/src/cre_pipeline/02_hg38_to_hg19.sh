#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] line $LINENO, exit code $?" >&2' ERR
shopt -s nullglob

NCPU="${SLURM_CPUS_PER_TASK:-8}"

log() { echo "[$(date +%F' '%T)] $*" >&2; }

# =========================================================
# G2 classified results: recover human hg38 from col4
# then liftover hg38 -> hg19 for LDSC
#
# Inputs:
#   classified/*.sfCRE.bed / *.sdCRE.bed / *.soCRE.bed
#     col1-3 = sheep V3
#     col4   = original human hg38 key: chr:start-end
#
#   lifted/*.unmapped
#     used as ssCRE source (original human hg38 intervals)
#
# Outputs:
#   classified_hg19/*.hg19.bed
#   ssCRE_hg19/*.ssCRE.hg19.bed
# =========================================================

# ---- tools ----
: "${LIFTOVER_BIN:?Set LIFTOVER_BIN in the paths configuration}"
: "${CHAIN_HG38_HG19:?Set CHAIN_HG38_HG19 in the paths configuration}"

# ---- params ----
MINMATCH_HG38_HG19="0.8"

# ---- dirs ----
: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR in the paths configuration}"
CLASS_DIR="${WORKDIR}/classified"
LIFTED_DIR="${WORKDIR}/lifted"

OUT_CLASS_DIR="${WORKDIR}/classified_hg19"
OUT_SS_DIR="${WORKDIR}/ssCRE_hg19"
TMP_DIR="${WORKDIR}/tmp_hg19"

mkdir -p "${OUT_CLASS_DIR}" "${OUT_SS_DIR}" "${TMP_DIR}"

[[ -x "${LIFTOVER_BIN}" ]] || { echo "ERROR: liftOver not found: ${LIFTOVER_BIN}" >&2; exit 1; }
[[ -f "${CHAIN_HG38_HG19}" ]] || { echo "ERROR: chain not found: ${CHAIN_HG38_HG19}" >&2; exit 1; }

export LIFTOVER_BIN CHAIN_HG38_HG19 MINMATCH_HG38_HG19
export OUT_CLASS_DIR OUT_SS_DIR TMP_DIR

# ---------------------------------------------------------
# class files: sf/sd/so
# ---------------------------------------------------------
process_class_file() {
    local f="$1"
    local base
    base=$(basename "$f" .bed)

    local tmp_hg38="${TMP_DIR}/${base}.hg38.bed"
    local out_hg19="${OUT_CLASS_DIR}/${base}.hg19.bed"
    local out_unmap="${OUT_CLASS_DIR}/${base}.hg38_to_hg19.unmapped"

    # Recover the original hg38 coordinates from column 4.
    awk 'BEGIN{FS=OFS="\t"}
         NF>=4 {
             split($4, a, ":")
             if (length(a) < 2) next
             split(a[2], b, "-")
             if (length(b) < 2) next
             print a[1], b[1], b[2], $4
         }' "$f" \
      | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u > "${tmp_hg38}"

    if [[ -s "${tmp_hg38}" ]]; then
        "${LIFTOVER_BIN}" -minMatch="${MINMATCH_HG38_HG19}" -bedPlus=4 \
            "${tmp_hg38}" "${CHAIN_HG38_HG19}" \
            "${out_hg19}" "${out_unmap}" || true
    else
        : > "${out_hg19}"
        : > "${out_unmap}"
    fi

    local n_in n_out n_un
    n_in=$(wc -l < "${tmp_hg38}")
    n_out=$(wc -l < "${out_hg19}")
    n_un=$(grep -cv '^#' "${out_unmap}" 2>/dev/null || echo 0)

    echo -e "class\t${base}\t${n_in}\t${n_out}\t${n_un}"
}
export -f process_class_file

# ---------------------------------------------------------
# ssCRE: from lifted/*.unmapped
# ---------------------------------------------------------
process_ss_file() {
    local f="$1"
    local base
    base=$(basename "$f" .unmapped)

    local tmp_hg38="${TMP_DIR}/${base}.ssCRE.hg38.bed"
    local out_hg19="${OUT_SS_DIR}/${base}.ssCRE.hg19.bed"
    local out_unmap="${OUT_SS_DIR}/${base}.ssCRE.hg38_to_hg19.unmapped"

    grep -v '^#' "$f" \
      | awk 'BEGIN{FS=OFS="\t"} NF>=4 {print $1,$2,$3,$4}' \
      | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u > "${tmp_hg38}"

    if [[ -s "${tmp_hg38}" ]]; then
        "${LIFTOVER_BIN}" -minMatch="${MINMATCH_HG38_HG19}" -bedPlus=4 \
            "${tmp_hg38}" "${CHAIN_HG38_HG19}" \
            "${out_hg19}" "${out_unmap}" || true
    else
        : > "${out_hg19}"
        : > "${out_unmap}"
    fi

    local n_in n_out n_un
    n_in=$(wc -l < "${tmp_hg38}")
    n_out=$(wc -l < "${out_hg19}")
    n_un=$(grep -cv '^#' "${out_unmap}" 2>/dev/null || echo 0)

    echo -e "ssCRE\t${base}\t${n_in}\t${n_out}\t${n_un}"
}
export -f process_ss_file

# =========================================================
# run class files in parallel
# =========================================================
log "===== Step 1: class files (sf/sd/so) ====="

CLASS_TASKS="${TMP_DIR}/class_tasks.txt"
: > "${CLASS_TASKS}"
for f in "${CLASS_DIR}"/*.sfCRE.bed "${CLASS_DIR}"/*.sdCRE.bed "${CLASS_DIR}"/*.soCRE.bed; do
    [[ -e "$f" ]] && echo "$f" >> "${CLASS_TASKS}"
done

log "class tasks: $(wc -l < "${CLASS_TASKS}")"
CLASS_SUMMARY="${TMP_DIR}/class_summary.tsv"
: > "${CLASS_SUMMARY}"

if [[ -s "${CLASS_TASKS}" ]]; then
    xargs -P "${NCPU}" -I {} bash -c 'process_class_file "$1"' _ {} \
        < "${CLASS_TASKS}" > "${CLASS_SUMMARY}"
fi

# =========================================================
# run ssCRE files in parallel
# =========================================================
log "===== Step 2: ssCRE files ====="

SS_TASKS="${TMP_DIR}/ss_tasks.txt"
: > "${SS_TASKS}"
for f in "${LIFTED_DIR}"/*.unmapped; do
    [[ -e "$f" ]] && echo "$f" >> "${SS_TASKS}"
done

log "ssCRE tasks: $(wc -l < "${SS_TASKS}")"
SS_SUMMARY="${TMP_DIR}/ss_summary.tsv"
: > "${SS_SUMMARY}"

if [[ -s "${SS_TASKS}" ]]; then
    xargs -P "${NCPU}" -I {} bash -c 'process_ss_file "$1"' _ {} \
        < "${SS_TASKS}" > "${SS_SUMMARY}"
fi

# =========================================================
# merge summary
# =========================================================
log "===== Step 3: merge summary ====="

FINAL_SUMMARY="${WORKDIR}/hg38_to_hg19_summary.tsv"
{
    echo -e "type\tname\thg38_input\thg19_output\thg19_unmapped"
    [[ -s "${CLASS_SUMMARY}" ]] && cat "${CLASS_SUMMARY}"
    [[ -s "${SS_SUMMARY}" ]] && cat "${SS_SUMMARY}"
} > "${FINAL_SUMMARY}"

log "Summary written to: ${FINAL_SUMMARY}"
log "classified hg19 dir: ${OUT_CLASS_DIR}"
log "ssCRE hg19 dir: ${OUT_SS_DIR}"
log "===== ALL DONE ====="
