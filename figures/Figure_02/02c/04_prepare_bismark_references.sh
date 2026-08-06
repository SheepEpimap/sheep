#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 4: prepare sheep and lambda Bismark references.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
SHEEP_FASTA="${SHEEP_FASTA:-/vol2/mengzhu/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna}"
LAMBDA_FASTA="${LAMBDA_FASTA:-${BASE}/extra_genome/lambda.fa.fasta}"
REFERENCE_ROOT="${REFERENCE_ROOT:-${BASE}/wgbs_bismark/ref}"

SHEEP_REFERENCE="${REFERENCE_ROOT}/sheep"
LAMBDA_REFERENCE="${REFERENCE_ROOT}/lambda"

[[ -s "${SHEEP_FASTA}" ]] || {
  echo "[ERROR] Sheep reference FASTA was not found: ${SHEEP_FASTA}" >&2
  exit 1
}
[[ -s "${LAMBDA_FASTA}" ]] || {
  echo "[ERROR] Lambda reference FASTA was not found: ${LAMBDA_FASTA}" >&2
  exit 1
}

for program in samtools bismark_genome_preparation; do
  command -v "${program}" >/dev/null 2>&1 || {
    echo "[ERROR] Required program was not found: ${program}" >&2
    exit 1
  }
done

mkdir -p "${SHEEP_REFERENCE}" "${LAMBDA_REFERENCE}"
ln -sfn "${SHEEP_FASTA}" "${SHEEP_REFERENCE}/genome.fa"

awk '
  /^>/ {
    if (!header_written) {
      print ">chrL"
      header_written = 1
    }
    next
  }
  { print }
' "${LAMBDA_FASTA}" > "${LAMBDA_REFERENCE}/genome.fa"

for reference_dir in "${SHEEP_REFERENCE}" "${LAMBDA_REFERENCE}"; do
  samtools faidx "${reference_dir}/genome.fa"
  cut -f1,2 "${reference_dir}/genome.fa.fai" \
    > "${reference_dir}/chrom.sizes"
  bismark_genome_preparation \
    --bowtie2 \
    "${reference_dir}"
done

echo "[DONE] Sheep Bismark reference: ${SHEEP_REFERENCE}"
echo "[DONE] Lambda Bismark reference: ${LAMBDA_REFERENCE}"
