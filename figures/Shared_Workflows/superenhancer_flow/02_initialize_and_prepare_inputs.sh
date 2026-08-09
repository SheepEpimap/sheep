#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd awk
require_cmd sort

log "Creating the anchored super-enhancer workspace: ${WORK_ROOT}"
mkdir -p \
  "${WORK_ROOT}" \
  "${WORK_ROOT}/annotation" \
  "${RESULT_DIR}" \
  "${STATE_ID_DIR}" \
  "${ENHANCER_COMBINE_DIR}" \
  "${RESULT_DIR}/AA_normal_enhancer/AA_one_count" \
  "${RESULT_DIR}/AA_normal_enhancer/AA_separate/AA_one_count" \
  "${RESULT_DIR}/AA_one_count" \
  "${RESULT_DIR}/Target_gene_independent" \
  "${RESULT_DIR}/Target_gene_cloest/AA_cluster/AA_Motif" \
  "${RESULT_DIR}/Target_gene_cloest/AA_cluster/human" \
  "${RESULT_DIR}/AA_Motif"

missing=0
for tissue in "${TISSUES[@]}"; do
  combined_gff="${WORK_ROOT}/${tissue}_enhancer.gff"
  combined_bed="${ENHANCER_COMBINE_DIR}/${tissue}_enhancer.bed"
  : > "${combined_gff}"
  : > "${combined_bed}"

  for state in "${STATES[@]}"; do
    input_bed="${STATE_DIR}/${tissue}_${state}.bed"
    output_gff="${WORK_ROOT}/${tissue}_${state}.gff"
    output_id="${STATE_ID_DIR}/${tissue}_${state}_id.bed"

    if [[ ! -s "${input_bed}" ]]; then
      log "Missing state BED: ${input_bed}"
      missing=1
      continue
    fi

    # Generate a stable genomic-region identifier while preserving the original
    # BED coordinates and optional fourth column.
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
      id=$1 ":" $2 "-" $3;
      original=(NF>=4 ? $4 : ".");
      print $1,$2,$3,id,original
    }' "${input_bed}" > "${output_id}"

    # ROSE GFF layout follows the source Markdown:
    # chrom, region_id, region_id, start, end, length, ., -, region_id
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
      id=$1 ":" $2 "-" $3;
      print $1,id,id,$2,$3,$3-$2,".","-",id
    }' "${input_bed}" > "${output_gff}"

    cat "${output_gff}" >> "${combined_gff}"
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {print $1,$2,$3,(NF>=4?$4:".")}' \
      "${input_bed}" >> "${combined_bed}"
  done

  if [[ -s "${combined_bed}" ]]; then
    sort -k1,1 -k2,2n -k3,3n "${combined_bed}" -o "${combined_bed}"
  fi
  log "Prepared ROSE input: ${combined_gff}"
done

if (( missing != 0 )); then
  die "One or more E6-E10 state BED files are absent. Review the messages above."
fi

log "Enhancer input preparation completed."
