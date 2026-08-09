#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
sample_list="${1:-${SAMPLE_MAP}}"; as_require_file "${sample_list}" "sample list"
while read -r sample _individual _rest; do [[ -n "${sample}" && "${sample}" != \#* && "${sample}" != sample_id ]] || continue; [[ -s "${ASE_TABLE_DIR}/${sample}.output.table" ]] && continue; sbatch --job-name="ase_${sample}" --output="${LOG_DIR}/ase_${sample}_%j.out" --error="${LOG_DIR}/ase_${sample}_%j.err" "${SCRIPT_DIR}/13_run_gatk_ase.sh" "${sample}"; done < "${sample_list}"
