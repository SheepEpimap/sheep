#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
sample_list="${1:-${SAMPLE_MAP}}"; mode="${2:-adjust}"; as_require_file "${sample_list}" "sample list"
while read -r sample _individual _rest; do [[ -n "${sample}" && "${sample}" != \#* && "${sample}" != sample_id ]] || continue; [[ "${mode}" == adjust && -s "${BINOMIAL_DIR}/${sample}_binomial_results.txt" ]] && continue; sbatch --job-name="fdr_${sample}" --output="${LOG_DIR}/fdr_${sample}_%j.out" --error="${LOG_DIR}/fdr_${sample}_%j.err" "${SCRIPT_DIR}/19_run_fdr.sh" "${sample}" "${mode}"; done < "${sample_list}"
