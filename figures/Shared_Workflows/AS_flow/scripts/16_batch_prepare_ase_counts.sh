#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
sample_list="${1:-${SAMPLE_MAP}}"; as_require_file "${sample_list}" "sample list"
while read -r sample _individual _rest; do [[ -n "${sample}" && "${sample}" != \#* && "${sample}" != sample_id ]] || continue; [[ -s "${ASE_COUNT_DIR}/${sample}.ase_counts.tsv" ]] && continue; sbatch --job-name="counts_${sample}" --output="${LOG_DIR}/counts_${sample}_%j.out" --error="${LOG_DIR}/counts_${sample}_%j.err" "${SCRIPT_DIR}/15_prepare_ase_counts.sh" "${sample}"; done < "${sample_list}"
