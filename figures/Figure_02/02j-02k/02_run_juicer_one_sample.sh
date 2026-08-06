#!/usr/bin/env bash

set -euo pipefail

# Figure 2j-2k, step 2: run the CPU Juicer pipeline for one sheep Hi-C sample.

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 <sample_id>" >&2
  exit 1
fi

SAMPLE_ID="$1"
JUICER_ROOT="${JUICER_ROOT:-/vol2/zhangshiwen/juicer}"
WORK_ROOT="${WORK_ROOT:-${JUICER_ROOT}/work}"
JUICER_SCRIPT="${JUICER_SCRIPT:-${JUICER_ROOT}/scripts/juicer.sh}"
SHEEP_FASTA="${SHEEP_FASTA:-${JUICER_ROOT}/references/sheep.fa}"
CHROM_SIZES="${CHROM_SIZES:-${JUICER_ROOT}/restriction_sites/sheep.chrom.sizes}"
SITE_FILE="${SITE_FILE:-${JUICER_ROOT}/restriction_sites/sheep_MboI.txt}"
ENZYME="${ENZYME:-MboI}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-16}}"

SAMPLE_DIR="${WORK_ROOT}/${SAMPLE_ID}"
FASTQ_DIR="${SAMPLE_DIR}/fastq"
LOG_DIR="${SAMPLE_DIR}/logs"

[[ -s "${JUICER_SCRIPT}" ]] || {
  echo "[ERROR] Juicer script was not found: ${JUICER_SCRIPT}" >&2
  exit 1
}
[[ -s "${SHEEP_FASTA}" ]] || {
  echo "[ERROR] Sheep reference was not found: ${SHEEP_FASTA}" >&2
  exit 1
}
[[ -s "${CHROM_SIZES}" && -s "${SITE_FILE}" ]] || {
  echo "[ERROR] Juicer chromosome sizes or restriction-site map is missing." >&2
  exit 1
}
[[ -d "${FASTQ_DIR}" ]] || {
  echo "[ERROR] Sample FASTQ directory was not found: ${FASTQ_DIR}" >&2
  exit 1
}

fastq_count="$(
  find "${FASTQ_DIR}" -maxdepth 1 -type f \
    \( -name '*.fastq' -o -name '*.fastq.gz' \) | wc -l
)"
(( fastq_count >= 2 )) || {
  echo "[ERROR] Fewer than two FASTQ files were found for ${SAMPLE_ID}." >&2
  exit 1
}

mkdir -p "${LOG_DIR}"

bash "${JUICER_SCRIPT}" \
  -z "${SHEEP_FASTA}" \
  -p "${CHROM_SIZES}" \
  -y "${SITE_FILE}" \
  -s "${ENZYME}" \
  -d "${SAMPLE_DIR}" \
  -D "${JUICER_ROOT}" \
  -t "${THREADS}" \
  > "${LOG_DIR}/juicer.log" 2>&1

[[ -s "${SAMPLE_DIR}/aligned/merged_nodups.txt" ]] || {
  echo "[ERROR] merged_nodups.txt was not generated for ${SAMPLE_ID}." >&2
  exit 1
}
[[ -s "${SAMPLE_DIR}/aligned/inter_30.hic" ]] || {
  echo "[ERROR] MAPQ30 Hi-C map was not generated for ${SAMPLE_ID}." >&2
  exit 1
}

touch "${SAMPLE_DIR}/JUICER_COMPLETE"
echo "[DONE] ${SAMPLE_ID}"
