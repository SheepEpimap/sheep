#!/usr/bin/env bash
#SBATCH -J BSSEQ
#SBATCH -D /vol2/zhangshiwen/rrbs
#SBATCH -p low
#SBATCH -c 14
#SBATCH --mem=28G
#SBATCH -t 2-00:00:00
#SBATCH -o /vol2/zhangshiwen/rrbs/wgbs_bismark/log/%x.%A_%a.out
#SBATCH -e /vol2/zhangshiwen/rrbs/wgbs_bismark/log/%x.%A_%a.err
#SBATCH --array=1-1

set -euo pipefail

# Figure 2c, step 7: SLURM worker for one row of samples.tsv.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
SAMPLESHEET="${SAMPLESHEET:-${BASE}/wgbs_bismark/samples.tsv}"
RUN_ONE="${RUN_ONE:-${BASE}/wgbs_bismark/scripts/06_run_bismark_one_sample.sh}"

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
  echo "[ERROR] SLURM_ARRAY_TASK_ID is not set." >&2
  exit 1
}
[[ -s "${SAMPLESHEET}" ]] || {
  echo "[ERROR] Sample sheet was not found: ${SAMPLESHEET}" >&2
  exit 1
}
[[ -x "${RUN_ONE}" ]] || {
  echo "[ERROR] Per-sample script is not executable: ${RUN_ONE}" >&2
  exit 1
}

mkdir -p "${BASE}/wgbs_bismark/log"
line_number="$((SLURM_ARRAY_TASK_ID + 1))"
IFS=$'\t' read -r sample_id library_type r1 r2 < <(
  awk -v n="${line_number}" 'NR == n { print; exit }' "${SAMPLESHEET}"
)

[[ -n "${sample_id:-}" ]] || {
  echo "[ERROR] No sample was found for array index ${SLURM_ARRAY_TASK_ID}." >&2
  exit 1
}

"${RUN_ONE}" "${sample_id}" "${library_type}" "${r1}" "${r2}"
