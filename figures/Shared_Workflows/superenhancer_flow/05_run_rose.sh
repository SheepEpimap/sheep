#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

ROSE_MAIN="${ROSE_HOME}/bin/ROSE_main.py"
require_file "${ROSE_MAIN}"
mkdir -p "${RESULT_DIR}"

export PYTHONPATH="${ROSE_HOME}/lib${PYTHONPATH:+:${PYTHONPATH}}"
export PATH="${ROSE_HOME}/bin:${PATH}"

# Set SKIP_EXISTING=0 to force regeneration of existing ROSE output folders.
SKIP_EXISTING="${SKIP_EXISTING:-1}"

for tissue in "${TISSUES[@]}"; do
  gff="${WORK_ROOT}/${tissue}_enhancer.gff"
  require_file "${gff}"

  for replicate in "${REPLICATES[@]}"; do
    bam="${WORK_ROOT}/H3K27ac_${tissue}_${replicate}.bowtie2.mapped.filtered.sort.bam"
    outdir="${RESULT_DIR}/${tissue}_enhancer_${replicate}"
    final_table="${outdir}/${tissue}_enhancer_SuperStitched.table.txt"
    require_file "${bam}"

    if [[ "${SKIP_EXISTING}" == 1 && -s "${final_table}" ]]; then
      log "Skipping existing ROSE result: ${final_table}"
      continue
    fi

    log "Running ROSE: tissue=${tissue}, replicate=${replicate}"
    "${ROSE_MAIN}" \
      -g "${ROSE_GENOME_KEY}" \
      -i "${gff}" \
      -r "${bam}" \
      -o "${outdir}" \
      -s "${ROSE_STITCH_DISTANCE}" \
      -t "${ROSE_TSS_EXCLUSION}"
  done
done

log "ROSE identification completed."
