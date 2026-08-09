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

as_require_dir "${RNASEQ_RAW_READS_DIR}" "RNA-seq raw-read directory"
as_require_cmds find sbatch

worker="${SCRIPT_DIR}/05_process_rnaseq.sh"
max_submit="${MAX_SUBMIT:-80}"
submitted=0

while IFS= read -r -d '' r1; do
  base="$(basename "${r1}")"
  sample="${base%_R1.fq.gz}"
  [[ "${sample}" == RNASeq* ]] || continue
  [[ -f "${RNASEQ_RAW_READS_DIR}/${sample}_R2.fq.gz" ]] || {
    as_warn "R2 missing; skip ${sample}"
    continue
  }
  if [[ -s "${RNASEQ_QUANT_DIR}/${sample}_salmon/quant.sf" &&
        -s "${RNASEQ_QUANT_DIR}/${sample}_gene_abundance.tsv" ]]; then
    as_info "Skip ${sample}: quantification outputs exist"
    continue
  fi
  (( submitted < max_submit )) || break
  sbatch \
    --job-name="${sample}_rnaseq" \
    --output="${LOG_DIR}/05_${sample}_%j.out" \
    --error="${LOG_DIR}/05_${sample}_%j.err" \
    "${worker}" "${sample}"
  submitted=$((submitted + 1))
  sleep 1
done < <(find "${RNASEQ_RAW_READS_DIR}" -maxdepth 1 -type f -name 'RNASeq*_R1.fq.gz' -print0 | sort -z)

as_info "Submitted ${submitted} RNA-seq jobs"
