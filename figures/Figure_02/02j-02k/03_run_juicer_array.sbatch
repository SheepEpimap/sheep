#!/usr/bin/env bash
#SBATCH -J sheepHiC
#SBATCH -D /vol2/zhangshiwen/juicer
#SBATCH -p low
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 4-00:00:00
#SBATCH -o /vol2/zhangshiwen/juicer/log/%x.%A_%a.out
#SBATCH -e /vol2/zhangshiwen/juicer/log/%x.%A_%a.err
#SBATCH --array=1-1

set -euo pipefail

# Figure 2j-2k, step 3: SLURM worker for one Hi-C sample ID.

JUICER_ROOT="${JUICER_ROOT:-/vol2/zhangshiwen/juicer}"
SAMPLE_FILE="${SAMPLE_FILE:-${JUICER_ROOT}/tissue.txt}"
RUN_ONE="${RUN_ONE:-${JUICER_ROOT}/figure_code/02_run_juicer_one_sample.sh}"

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
  echo "[ERROR] SLURM_ARRAY_TASK_ID is not set." >&2
  exit 1
}
[[ -s "${SAMPLE_FILE}" ]] || {
  echo "[ERROR] Sample file was not found: ${SAMPLE_FILE}" >&2
  exit 1
}
[[ -x "${RUN_ONE}" ]] || {
  echo "[ERROR] Per-sample Juicer script is not executable: ${RUN_ONE}" >&2
  exit 1
}

sample_id="$(
  awk -v n="${SLURM_ARRAY_TASK_ID}" 'NF && ++i == n { print $1; exit }' \
    "${SAMPLE_FILE}"
)"
[[ -n "${sample_id}" ]] || {
  echo "[ERROR] No sample was found for array index ${SLURM_ARRAY_TASK_ID}." >&2
  exit 1
}

"${RUN_ONE}" "${sample_id}"
