#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] line $LINENO, exit code $?" >&2' ERR
shopt -s nullglob

# Collect the E5 hg19 CRE files used by G2_TSR_label.
# This makes the previously manual E5_hg19_all step explicit and reproducible.

: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR in the paths configuration}"
CLASS_HG19_DIR="${CLASS_HG19_DIR:-${WORKDIR}/classified_hg19}"
SS_HG19_DIR="${SS_HG19_DIR:-${WORKDIR}/ssCRE_hg19}"
OUT_DIR="${OUT_DIR:-${WORKDIR}/E5_hg19_all}"

mkdir -p "${OUT_DIR}"

rm -f "${OUT_DIR}"/*.sfCRE.hg19.bed \
      "${OUT_DIR}"/*.sdCRE.hg19.bed \
      "${OUT_DIR}"/*.soCRE.hg19.bed \
      "${OUT_DIR}"/*.ssCRE.hg19.bed

copied=0

for f in "${CLASS_HG19_DIR}"/*_E5.sfCRE.hg19.bed \
         "${CLASS_HG19_DIR}"/*_E5.sdCRE.hg19.bed \
         "${CLASS_HG19_DIR}"/*_E5.soCRE.hg19.bed; do
    [[ -e "${f}" ]] || continue
    cp "${f}" "${OUT_DIR}/"
    copied=$((copied + 1))
done

for f in "${SS_HG19_DIR}"/*_E5.ssCRE.hg19.bed; do
    [[ -e "${f}" ]] || continue
    cp "${f}" "${OUT_DIR}/"
    copied=$((copied + 1))
done

echo "[INFO] class_hg19_dir=${CLASS_HG19_DIR}"
echo "[INFO] ss_hg19_dir=${SS_HG19_DIR}"
echo "[INFO] out_dir=${OUT_DIR}"
echo "[INFO] copied_files=${copied}"

if (( copied == 0 )); then
    echo "[ERROR] no E5 hg19 CRE files copied; check upstream classified_hg19/ssCRE_hg19 outputs" >&2
    exit 1
fi

{
    echo -e "cre_type\tfiles\trows"
    for cre in sfCRE sdCRE soCRE ssCRE; do
        files=( "${OUT_DIR}"/*.${cre}.hg19.bed )
        rows=0
        for f in "${files[@]}"; do
            [[ -s "${f}" ]] && rows=$((rows + $(wc -l < "${f}")))
        done
        echo -e "${cre}\t${#files[@]}\t${rows}"
    done
} > "${OUT_DIR}/E5_hg19_all_summary.tsv"

echo "[DONE] summary=${OUT_DIR}/E5_hg19_all_summary.tsv"
