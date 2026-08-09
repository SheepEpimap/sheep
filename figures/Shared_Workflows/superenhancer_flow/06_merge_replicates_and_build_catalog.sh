#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd bedtools
require_cmd awk
require_cmd sort
mkdir -p "${RESULT_DIR}"
cd "${RESULT_DIR}"

summary="${RESULT_DIR}/all_super_enhancer_summary.txt"
all_regions="${RESULT_DIR}/all_super_enhancer.csv"
printf 'Tissue\tNumber\tMean_length_bp\tMean_ROSE_signal\tTotal_length_bp\tGenome_fraction\n' > "${summary}"
: > "${all_regions}"

for tissue in "${TISSUES[@]}"; do
  tmp="$(mktemp "${RESULT_DIR}/.${tissue}.replicates.XXXXXX")"
  trap 'rm -f "${tmp:-}"' EXIT

  for replicate in "${REPLICATES[@]}"; do
    table="${RESULT_DIR}/${tissue}_enhancer_${replicate}/${tissue}_enhancer_SuperStitched.table.txt"
    require_file "${table}"
    # The first six ROSE lines are metadata/header lines in the source workflow.
    sed '1,6d' "${table}" | awk 'BEGIN{OFS="\t"} NF>=7 {print $2,$3,$4,$7}' >> "${tmp}"
  done

  output="${RESULT_DIR}/${tissue}_super_enhancer.bed"
  sort -k1,1 -k2,2n -k3,3n "${tmp}" \
    | bedtools merge -i stdin -c 4 -o mean \
    | awk -v t="${tissue}" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,t}' > "${output}"
  rm -f "${tmp}"
  trap - EXIT

  n=$(awk 'END{print NR+0}' "${output}")
  mean_len=$(awk '{s+=($3-$2)} END{if(NR) printf "%.6f",s/NR; else print "NA"}' "${output}")
  mean_signal=$(awk '{s+=$4} END{if(NR) printf "%.6f",s/NR; else print "NA"}' "${output}")
  total_len=$(awk '{s+=($3-$2)} END{print s+0}' "${output}")
  fraction=$(safe_ratio "${total_len}" "${GENOME_SIZE_BP}" 10)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tissue}" "${n}" "${mean_len}" "${mean_signal}" "${total_len}" "${fraction}" >> "${summary}"
  cat "${output}" >> "${all_regions}"
done

sort -k1,1 -k2,2n -k3,3n "${all_regions}" | bedtools merge -i stdin \
  > "${RESULT_DIR}/all_super_enhancer_combine.csv"

# Negative -d values in bedtools merge require a minimum amount of overlap.
# These names are retained because downstream commands in the source Markdown
# reference them directly.
sort -k1,1 -k2,2n -k3,3n "${all_regions}" | bedtools merge -i stdin -d -100 \
  > "${RESULT_DIR}/all_super_enhancer_combine_1.csv"
sort -k1,1 -k2,2n -k3,3n "${all_regions}" | bedtools merge -i stdin -d -10000 -c 1 -o count \
  > "${RESULT_DIR}/all_super_enhancer_combine_10000.csv"

log "Merged replicate-specific ROSE calls and built the cross-tissue catalog."
