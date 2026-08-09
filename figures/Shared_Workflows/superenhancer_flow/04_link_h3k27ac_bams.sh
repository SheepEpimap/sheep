#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

mkdir -p "${WORK_ROOT}"
missing=0

for tissue in "${TISSUES[@]}"; do
  for replicate in "${REPLICATES[@]}"; do
    bam="${H3K27AC_BAM_DIR}/H3K27ac_${tissue}_${replicate}.bowtie2.mapped.filtered.sort.bam"
    bai="${bam}.bai"

    if [[ ! -s "${bam}" ]]; then
      log "Missing BAM: ${bam}"
      missing=1
      continue
    fi

    ln -sfn "${bam}" "${WORK_ROOT}/$(basename "${bam}")"
    if [[ -s "${bai}" ]]; then
      ln -sfn "${bai}" "${WORK_ROOT}/$(basename "${bai}")"
    else
      log "BAM index not found (ROSE may still run if the BAM is indexed elsewhere): ${bai}"
    fi
  done
done

if (( missing != 0 )); then
  die "One or more expected H3K27ac BAM files are absent."
fi
log "H3K27ac BAM links created below ${WORK_ROOT}."
