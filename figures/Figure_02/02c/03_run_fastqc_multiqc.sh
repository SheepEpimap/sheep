#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 3: run raw-read FastQC and aggregate reports with MultiQC.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
FASTQDIR="${FASTQDIR:-${BASE}/fastq}"
OUTDIR="${OUTDIR:-${BASE}/qc}"
THREADS="${THREADS:-8}"

command -v fastqc >/dev/null 2>&1 || {
  echo "[ERROR] fastqc was not found." >&2
  exit 1
}
command -v multiqc >/dev/null 2>&1 || {
  echo "[ERROR] multiqc was not found." >&2
  exit 1
}

mkdir -p "${OUTDIR}/fastqc"

shopt -s nullglob
r1_files=("${FASTQDIR}"/*_1.fastq.gz)
if (( ${#r1_files[@]} == 0 )); then
  echo "[ERROR] No *_1.fastq.gz files were found in ${FASTQDIR}." >&2
  exit 1
fi

for r1 in "${r1_files[@]}"; do
  r2="${r1%_1.fastq.gz}_2.fastq.gz"
  if [[ ! -s "${r2}" ]]; then
    echo "[WARN] Paired R2 file was not found for: ${r1}" >&2
    continue
  fi
  fastqc -t "${THREADS}" -o "${OUTDIR}/fastqc" "${r1}" "${r2}"
done

multiqc --force -o "${OUTDIR}" "${OUTDIR}/fastqc"
echo "[DONE] MultiQC report: ${OUTDIR}/multiqc_report.html"
