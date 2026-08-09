#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd bedtools
require_cmd paste
require_cmd awk

build_matrix() {
  local catalog="$1" output_counts="$2" output_binary="$3" bed_pattern="$4" count_dir="$5"
  require_file "${catalog}"
  mkdir -p "${count_dir}"

  local tmpdir
  tmpdir=$(mktemp -d "${RESULT_DIR}/.matrix.XXXXXX")
  trap 'rm -rf "${tmpdir:-}"' RETURN

  cut -f1-3 "${catalog}" > "${tmpdir}/coords.tsv"
  printf 'chr\tstart\tend' > "${tmpdir}/header.tsv"
  local cols=("${tmpdir}/coords.tsv")

  for tissue in "${TISSUES[@]}"; do
    local bed count_bed column
    printf -v bed "${bed_pattern}" "${tissue}"
    require_file "${bed}"
    count_bed="${count_dir}/$(basename "${bed}" .bed)_Gs.bed"
    column="${tmpdir}/${tissue}.count"
    bedtools intersect \
      -a <(sort -k1,1 -k2,2n "${catalog}") \
      -b <(sort -k1,1 -k2,2n "${bed}") \
      -c -sorted > "${count_bed}"
    awk '{print $NF}' "${count_bed}" > "${column}"
    cols+=("${column}")
    printf '\t%s' "${tissue}" >> "${tmpdir}/header.tsv"
  done
  printf '\n' >> "${tmpdir}/header.tsv"

  paste "${cols[@]}" > "${tmpdir}/matrix.counts.tsv"
  cat "${tmpdir}/header.tsv" "${tmpdir}/matrix.counts.tsv" > "${output_counts}"

  {
    awk 'BEGIN{OFS="\t"} NR==1 {print $0,"tissue_count"; next}
      {
        sum=0;
        printf "%s\t%s\t%s",$1,$2,$3;
        for(i=4;i<=NF;i++) {v=($i>0?1:0); sum+=v; printf "\t%d",v}
        printf "\t%d\n",sum
      }' "${output_counts}"
  } > "${output_binary}"

  rm -rf "${tmpdir}"
  trap - RETURN
}

# ---- Super-enhancer matrix --------------------------------------------------
super_catalog="${RESULT_DIR}/all_super_enhancer_combine_10000.csv"
build_matrix \
  "${super_catalog}" \
  "${RESULT_DIR}/AA_one_count/all_super_enhancer_Gs.csv" \
  "${RESULT_DIR}/AA_one_count/all_super_enhancer_Gs_one_count.csv" \
  "${RESULT_DIR}/%s_super_enhancer.bed" \
  "${RESULT_DIR}/AA_one_count"

# ---- Normal-enhancer matrix -------------------------------------------------
normal_dir="${RESULT_DIR}/AA_normal_enhancer"
normal_catalog="${normal_dir}/all_normal_enhancer_combine.csv"
cat "${normal_dir}"/*_enhancer_normal.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > "${normal_catalog}"
build_matrix \
  "${normal_catalog}" \
  "${normal_dir}/AA_one_count/all_normal_enhancer_Gs.csv" \
  "${normal_dir}/AA_one_count/all_normal_enhancer_Gs_one_count.csv" \
  "${normal_dir}/%s_enhancer_normal.bed" \
  "${normal_dir}/AA_one_count"

# ---- State-separated normal enhancers --------------------------------------
separate_dir="${normal_dir}/AA_separate"
mkdir -p "${separate_dir}/AA_one_count"
for state in "${STATES[@]}"; do
  for tissue in "${TISSUES[@]}"; do
    bedtools intersect -v \
      -a "${STATE_DIR}/${tissue}_${state}.bed" \
      -b "${RESULT_DIR}/${tissue}_super_enhancer.bed" \
      > "${separate_dir}/${tissue}_${state}.bed"
  done

  catalog="${separate_dir}/all_normal_${state}_combine.csv"
  cat "${separate_dir}"/*_"${state}".bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin -c 4 -o distinct > "${catalog}"

  build_matrix \
    "${catalog}" \
    "${separate_dir}/AA_one_count/all_normal_${state}_Gs.csv" \
    "${separate_dir}/AA_one_count/all_normal_${state}_Gs_one_count.csv" \
    "${separate_dir}/%s_${state}.bed" \
    "${separate_dir}/AA_one_count"
done

log "Super-, normal-, and E6-E10 presence matrices completed."
