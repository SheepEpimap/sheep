#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd bedtools
normal_dir="${RESULT_DIR}/AA_normal_enhancer"
mkdir -p "${normal_dir}"

summary_split="${normal_dir}/all_super_enhancer_summary_2.txt"
summary_normal_state="${normal_dir}/all_super_enhancer_summary_normal.txt"
summary_super_state="${RESULT_DIR}/all_super_enhancer_summary_1.txt"
printf 'Tissue\tAll\tNormal_number\tSuper_number\tNormal_percent\tSuper_percent\n' > "${summary_split}"
printf 'Tissue\tAll_number\tE6n\tE7n\tE8n\tE9n\tE10n\tE6\tE7\tE8\tE9\tE10\n' > "${summary_normal_state}"
printf 'Tissue\tNumber\tSuper_size\tOrigin_size\tAll_number\tE6n\tE7n\tE8n\tE9n\tE10n\tSuper_average_size\tE6\tE7\tE8\tE9\tE10\n' > "${summary_super_state}"

for tissue in "${TISSUES[@]}"; do
  enhancer="${ENHANCER_COMBINE_DIR}/${tissue}_enhancer.bed"
  super="${RESULT_DIR}/${tissue}_super_enhancer.bed"
  require_file "${enhancer}"
  require_file "${super}"

  normal="${normal_dir}/${tissue}_enhancer_normal.bed"
  inside="${normal_dir}/${tissue}_enhancer_super.bed"
  bedtools intersect -v -a "${enhancer}" -b "${super}" > "${normal}"
  bedtools intersect -wo -a "${enhancer}" -b "${super}" > "${inside}"

  all_n=$(awk 'END{print NR+0}' "${enhancer}")
  normal_n=$(awk 'END{print NR+0}' "${normal}")
  super_components_n=$(awk 'END{print NR+0}' "${inside}")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tissue}" "${all_n}" "${normal_n}" "${super_components_n}" \
    "$(safe_ratio "${normal_n}" "${all_n}" 4)" \
    "$(safe_ratio "${super_components_n}" "${all_n}" 4)" >> "${summary_split}"

  normal_counts=()
  super_counts=()
  for state in "${STATES[@]}"; do
    state_bed="${STATE_DIR}/${tissue}_${state}.bed"
    require_file "${state_bed}"

    normal_count_bed="${normal_dir}/${tissue}_enhancer_normal_${state}_number.bed"
    super_count_bed="${normal_dir}/${tissue}_enhancer_super_${state}_number.bed"
    bedtools intersect -c -b "${state_bed}" -a "${normal}" > "${normal_count_bed}"
    bedtools intersect -c -b "${state_bed}" -a "${super}" > "${super_count_bed}"
    normal_counts+=("$(awk '{s+=$NF}END{print s+0}' "${normal_count_bed}")")
    super_counts+=("$(awk '{s+=$NF}END{print s+0}' "${super_count_bed}")")
  done

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tissue}" "${normal_n}" \
    "${normal_counts[0]}" "${normal_counts[1]}" "${normal_counts[2]}" "${normal_counts[3]}" "${normal_counts[4]}" \
    "$(safe_ratio "${normal_counts[0]}" "${normal_n}" 4)" \
    "$(safe_ratio "${normal_counts[1]}" "${normal_n}" 4)" \
    "$(safe_ratio "${normal_counts[2]}" "${normal_n}" 4)" \
    "$(safe_ratio "${normal_counts[3]}" "${normal_n}" 4)" \
    "$(safe_ratio "${normal_counts[4]}" "${normal_n}" 4)" >> "${summary_normal_state}"

  super_n=$(awk 'END{print NR+0}' "${super}")
  mean_super_len=$(awk '{s+=($3-$2)}END{if(NR)printf "%.6f",s/NR;else print "NA"}' "${super}")
  mean_origin_len=$(awk '{s+=($3-$2)}END{if(NR)printf "%.6f",s/NR;else print "NA"}' "${inside}")
  component_count_file="${normal_dir}/${tissue}_enhancer_super_number.bed"
  bedtools intersect -c -b "${enhancer}" -a "${super}" > "${component_count_file}"
  component_total=$(awk '{s+=$NF}END{print s+0}' "${component_count_file}")
  mean_component_count=$(awk '{s+=$NF}END{if(NR)printf "%.10f",s/NR;else print 0}' "${component_count_file}")
  super_average_size=$(safe_ratio "${mean_super_len}" "${mean_component_count}" 4)

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tissue}" "${super_n}" "${mean_super_len}" "${mean_origin_len}" "${component_total}" \
    "${super_counts[0]}" "${super_counts[1]}" "${super_counts[2]}" "${super_counts[3]}" "${super_counts[4]}" \
    "${super_average_size}" \
    "$(safe_ratio "${super_counts[0]}" "${component_total}" 4)" \
    "$(safe_ratio "${super_counts[1]}" "${component_total}" 4)" \
    "$(safe_ratio "${super_counts[2]}" "${component_total}" 4)" \
    "$(safe_ratio "${super_counts[3]}" "${component_total}" 4)" \
    "$(safe_ratio "${super_counts[4]}" "${component_total}" 4)" >> "${summary_super_state}"
done

log "Normal-versus-super-enhancer data tables completed."
