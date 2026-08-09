#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 6
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
sample_id="${1:-}"; mode="${2:-adjust}"; [[ -n "${sample_id}" ]] || as_die "Usage: 19_run_fdr.sh SAMPLE_ID [adjust|simulation]"
as_activate_conda "${ASE_ENV}"; as_require_cmds Rscript; input="${ASE_COUNT_DIR}/${sample_id}.ase_counts.tsv"; as_require_file "${input}" "ASE count table"
case "${mode}" in adjust) Rscript "${SCRIPT_DIR}/17_fdr_adjust.R" "${sample_id}" "${input}" "${BINOMIAL_DIR}" "${FDR_DETAIL_DIR}";; simulation) Rscript "${SCRIPT_DIR}/18_fdr_simulation.R" "${sample_id}" "${input}" "${BINOMIAL_DIR}" "${FDR_SIMULATIONS:-1000}" "${FDR_CORES:-6}";; *) as_die "Unknown mode: ${mode}";; esac
