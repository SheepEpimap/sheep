#!/usr/bin/env bash

set -euo pipefail

# Figure 2j-2k, step 1: prepare the sheep Juicer reference and MboI map.

JUICER_ROOT="${JUICER_ROOT:-/vol2/zhangshiwen/juicer}"
REFERENCE_FASTA="${REFERENCE_FASTA:-/vol2/mengzhu/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna}"
SHEEP_FASTA="${SHEEP_FASTA:-${JUICER_ROOT}/references/sheep.fa}"
ENZYME="${ENZYME:-MboI}"
GENOME_ID="${GENOME_ID:-sheep}"
SITE_GENERATOR="${SITE_GENERATOR:-/vol2/zhangshiwen/miniconda3/envs/hic/share/juicer-1.6-0/scripts/misc/generate_site_positions.py}"
THREADS="${THREADS:-16}"

REFERENCE_DIR="${JUICER_ROOT}/references"
SITE_DIR="${JUICER_ROOT}/restriction_sites"
SITE_FILE="${SITE_DIR}/${GENOME_ID}_${ENZYME}.txt"
CHROM_SIZES="${SITE_DIR}/${GENOME_ID}.chrom.sizes"

[[ -s "${REFERENCE_FASTA}" ]] || {
  echo "[ERROR] Sheep reference FASTA was not found: ${REFERENCE_FASTA}" >&2
  exit 1
}
[[ -s "${SITE_GENERATOR}" ]] || {
  echo "[ERROR] Juicer restriction-site generator was not found: ${SITE_GENERATOR}" >&2
  exit 1
}

for program in bwa samtools python; do
  command -v "${program}" >/dev/null 2>&1 || {
    echo "[ERROR] Required program was not found: ${program}" >&2
    exit 1
  }
done

mkdir -p "${REFERENCE_DIR}" "${SITE_DIR}" "${JUICER_ROOT}/work"
ln -sfn "${REFERENCE_FASTA}" "${SHEEP_FASTA}"

if [[ ! -s "${SHEEP_FASTA}.bwt" ]]; then
  bwa index -p "${SHEEP_FASTA}" "${SHEEP_FASTA}"
fi

samtools faidx "${SHEEP_FASTA}"
cut -f1,2 "${SHEEP_FASTA}.fai" > "${CHROM_SIZES}"

(
  cd "${SITE_DIR}"
  python "${SITE_GENERATOR}" \
    "${ENZYME}" \
    "${GENOME_ID}" \
    "${SHEEP_FASTA}"
)

[[ -s "${SITE_FILE}" ]] || {
  echo "[ERROR] Restriction-site map was not generated: ${SITE_FILE}" >&2
  exit 1
}

echo "[DONE] Reference FASTA: ${SHEEP_FASTA}"
echo "[DONE] Chromosome sizes: ${CHROM_SIZES}"
echo "[DONE] Restriction-site map: ${SITE_FILE}"
