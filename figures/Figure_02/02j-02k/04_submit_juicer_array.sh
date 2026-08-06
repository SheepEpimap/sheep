#!/usr/bin/env bash

set -euo pipefail

# Figure 2j-2k, step 4: submit the sheep Hi-C samples as a SLURM array.

JUICER_ROOT="${JUICER_ROOT:-/vol2/zhangshiwen/juicer}"
SAMPLE_FILE="${SAMPLE_FILE:-${JUICER_ROOT}/tissue.txt}"
SBATCH_SCRIPT="${SBATCH_SCRIPT:-${JUICER_ROOT}/figure_code/03_run_juicer_array.sbatch}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"

[[ -s "${SAMPLE_FILE}" ]] || {
  echo "[ERROR] Sample file was not found: ${SAMPLE_FILE}" >&2
  exit 1
}
[[ -s "${SBATCH_SCRIPT}" ]] || {
  echo "[ERROR] SLURM worker was not found: ${SBATCH_SCRIPT}" >&2
  exit 1
}

sample_count="$(awk 'NF { n++ } END { print n + 0 }' "${SAMPLE_FILE}")"
(( sample_count > 0 )) || {
  echo "[ERROR] The sample file contains no sample IDs." >&2
  exit 1
}

mkdir -p "${JUICER_ROOT}/log"
sbatch \
  --array="1-${sample_count}%${MAX_CONCURRENT}" \
  "${SBATCH_SCRIPT}"

echo "[DONE] Submitted array: 1-${sample_count}%${MAX_CONCURRENT}"
