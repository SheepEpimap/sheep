#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 8: submit all sample rows as a bounded SLURM array.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
SAMPLESHEET="${SAMPLESHEET:-${BASE}/wgbs_bismark/samples.tsv}"
SBATCH_SCRIPT="${SBATCH_SCRIPT:-${BASE}/wgbs_bismark/scripts/07_run_bismark_array.sbatch}"
MAX_CONCURRENT="${MAX_CONCURRENT:-10}"

[[ -s "${SAMPLESHEET}" ]] || {
  echo "[ERROR] Sample sheet was not found: ${SAMPLESHEET}" >&2
  exit 1
}
[[ -s "${SBATCH_SCRIPT}" ]] || {
  echo "[ERROR] SLURM worker script was not found: ${SBATCH_SCRIPT}" >&2
  exit 1
}

sample_count="$(( $(wc -l < "${SAMPLESHEET}") - 1 ))"
(( sample_count > 0 )) || {
  echo "[ERROR] The sample sheet contains no samples." >&2
  exit 1
}

sbatch \
  --array="1-${sample_count}%${MAX_CONCURRENT}" \
  "${SBATCH_SCRIPT}"

echo "[DONE] Submitted array: 1-${sample_count}%${MAX_CONCURRENT}"
