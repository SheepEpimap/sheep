#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
sample_list="${1:-${SAMPLE_MAP}}"; max_submit="${MAX_JOBS:-100}"; as_require_file "${sample_list}" "sample list"; count=0
while read -r sample _individual _rest; do
  [[ -n "${sample}" && "${sample}" != \#* && "${sample}" != sample_id ]] || continue
  [[ -s "${WASP_FINAL_DIR}/${sample}.keep.merged.sorted.bam" ]] && { as_info "Skipping completed ${sample}"; continue; }
  (( count < max_submit )) || break
  sbatch --job-name="wasp_${sample}" --output="${LOG_DIR}/wasp_${sample}_%j.out" --error="${LOG_DIR}/wasp_${sample}_%j.err" "${SCRIPT_DIR}/10_run_wasp.sh" "${sample}"; count=$((count+1))
done < "${sample_list}"
as_info "Submitted ${count} WASP jobs"
