#!/usr/bin/env bash
#SBATCH -p low
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 1-00:00:00

set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_config.sh"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_lib.sh"

as_require_cmds find sbatch
as_require_dir "${INPUT_BAM_DIR}" "input BAM directory"

max_submit="${MAX_SUBMIT:-80}"
submitted=0
worker="${SCRIPT_DIR}/03_process_mappability_and_peaks.sh"

submit_one() {
  local sample="$1" suffix="$2" kind="$3"
  local final_peak="${PEAK_DIR}/${sample}_Peaks.bed"
  local out_log="${LOG_DIR}/03_${sample}_%j.out"
  local err_log="${LOG_DIR}/03_${sample}_%j.err"

  if [[ -s "${final_peak}" ]]; then
    as_info "Skip ${sample}: final peak exists"
    return 0
  fi
  if (( submitted >= max_submit )); then
    return 1
  fi

  as_info "Submit ${kind}: ${sample}"
  sbatch \
    --job-name="${sample}_peak" \
    --output="${out_log}" \
    --error="${err_log}" \
    "${worker}" "${sample}" "${suffix}"
  submitted=$((submitted + 1))
  sleep 1
}

# Non-ATAC peak assays. RNASeq is explicitly excluded because it has no peak stage.
while IFS= read -r -d '' bam; do
  basename_file="$(basename "${bam}")"
  sample="${basename_file%.${NON_ATAC_BAM_SUFFIX}}"
  [[ "${sample}" == "${basename_file}" ]] && continue
  [[ "${sample}" == RNASeq* ]] && continue
  [[ "${sample}" == ATAC* ]] && continue
  submit_one "${sample}" "${NON_ATAC_BAM_SUFFIX}" "non-ATAC peak assay" || break
done < <(find "${INPUT_BAM_DIR}" -maxdepth 1 -type f \
  -name "*.${NON_ATAC_BAM_SUFFIX}" \
  ! -name "*raw_sorted.bam" -print0 | sort -z)

if (( submitted < max_submit )); then
  while IFS= read -r -d '' bam; do
    basename_file="$(basename "${bam}")"
    sample="${basename_file%.${ATAC_BAM_SUFFIX}}"
    [[ "${sample}" == ATAC* ]] || continue
    submit_one "${sample}" "${ATAC_BAM_SUFFIX}" "ATAC" || break
  done < <(find "${INPUT_BAM_DIR}" -maxdepth 1 -type f \
    -name "ATAC*.${ATAC_BAM_SUFFIX}" -print0 | sort -z)
fi

as_info "Submitted ${submitted} peak-pipeline jobs (limit=${max_submit})"
