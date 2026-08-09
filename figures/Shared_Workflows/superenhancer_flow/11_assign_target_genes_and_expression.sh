#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd bedtools
require_file "${PROTEIN_TSS_BED}"
require_file "${EXPRESSION_TPM}"
require_file "${EXPRESSION_TAU}"
require_file "${SHEEP_HUMAN_CONSERVATION}"

independent_dir="${RESULT_DIR}/Target_gene_independent"
global_dir="${RESULT_DIR}/Target_gene_cloest"
mkdir -p "${independent_dir}" "${global_dir}"

join_by_first_field() {
  local left="$1" right="$2" output="$3"
  join -t $'\t' -1 1 -2 1 \
    <(sed 's/ /\t/g' "${left}" | sort -t $'\t' -k1,1) \
    <(sed 's/ /\t/g' "${right}" | sort -t $'\t' -k1,1) \
    > "${output}" || true
}

add_expression_summary() {
  local input="$1" group="$2" output="$3"
  awk -v first="${EXPRESSION_FIRST_COLUMN}" -v last="${EXPRESSION_LAST_COLUMN}" -v group="${group}" \
    'BEGIN{OFS="\t"}
     {
       sum=0; n=0;
       for(i=first;i<=last && i<=NF;i++){sum+=$i; n++}
       printf "%s",$1;
       for(i=2;i<=6 && i<=NF;i++) printf OFS "%s",$i;
       printf OFS "%.10g" OFS "%.10g" OFS "%s\n",sum,(n?sum/n:0),group
     }' "${input}" > "${output}"
}

all_tissue_expression="${independent_dir}/all_super_enhancer_gene_expression.txt"
: > "${all_tissue_expression}"

for tissue in "${TISSUES[@]}"; do
  super="${RESULT_DIR}/${tissue}_super_enhancer.bed"
  closest="${independent_dir}/all_${tissue}_super_enhancer_gene.bed"
  gene_with="${independent_dir}/gene_with_${tissue}_super_enhancer.txt"
  gene_without="${independent_dir}/gene_without_${tissue}_super_enhancer.txt"
  non_bed="${independent_dir}/non_${tissue}_super_enhancer_gene.bed"
  require_file "${super}"

  bedtools closest \
    -a <(sort -k1,1 -k2,2n "${super}") \
    -b <(sort -k1,1 -k2,2n "${PROTEIN_TSS_BED}") \
    -D b > "${closest}"

  # Super-enhancer BED has five columns; the gene identifier is field 4 of the
  # TSS BED, therefore field 9 in the combined closest output.
  awk 'BEGIN{OFS="\t"} NF>=9{print $9}' "${closest}" | awk 'NF&&!seen[$1]++' > "${gene_with}"
  awk 'NR==FNR{seen[$1]=1;next}!seen[$4]' "${gene_with}" "${PROTEIN_TSS_BED}" > "${non_bed}"
  awk 'NF>=4{print $4}' "${non_bed}" | awk 'NF&&!seen[$1]++' > "${gene_without}"

  with_expr="${independent_dir}/.${tissue}.with_expr.tsv"
  without_expr="${independent_dir}/.${tissue}.without_expr.tsv"
  join_by_first_field "${gene_with}" "${EXPRESSION_TPM}" "${with_expr}"
  join_by_first_field "${gene_without}" "${EXPRESSION_TPM}" "${without_expr}"

  add_expression_summary "${with_expr}" "WithSuperEnhancer" "${independent_dir}/.${tissue}.with_summary.tsv"
  add_expression_summary "${without_expr}" "WithoutSuperEnhancer" "${independent_dir}/.${tissue}.without_summary.tsv"
  cat "${independent_dir}/.${tissue}.with_summary.tsv" \
      "${independent_dir}/.${tissue}.without_summary.tsv" \
      > "${independent_dir}/AA_with_without_${tissue}_super_enhancer_expression.txt"

  awk -v t="${tissue}" 'BEGIN{OFS="\t"}{$(NF+1)=t;print}' "${with_expr}" >> "${all_tissue_expression}"

  join_by_first_field "${gene_with}" "${EXPRESSION_TAU}" "${independent_dir}/.${tissue}.with_tau.tsv"
  join_by_first_field "${gene_without}" "${EXPRESSION_TAU}" "${independent_dir}/.${tissue}.without_tau.tsv"
  join_by_first_field "${independent_dir}/.${tissue}.with_tau.tsv" "${SHEEP_HUMAN_CONSERVATION}" "${independent_dir}/.${tissue}.with_human.tsv"
  join_by_first_field "${independent_dir}/.${tissue}.without_tau.tsv" "${SHEEP_HUMAN_CONSERVATION}" "${independent_dir}/.${tissue}.without_human.tsv"
  cat "${independent_dir}/.${tissue}.with_human.tsv" \
      "${independent_dir}/.${tissue}.without_human.tsv" \
      > "${independent_dir}/AA_with_without_${tissue}_super_enhancer_compare_to_human.txt"
done

{
  header=$(head -n 1 "${EXPRESSION_TPM}")
  printf '%s\ttissue\n' "${header}"
  cat "${all_tissue_expression}"
} > "${independent_dir}/all_super_enhancer_gene_expression_last.txt"

# ---- Global with/without-super-enhancer comparison --------------------------
global_catalog="${RESULT_DIR}/all_super_enhancer_combine_10000.csv"
require_file "${global_catalog}"
bedtools closest \
  -a <(sort -k1,1 -k2,2n "${global_catalog}") \
  -b <(sort -k1,1 -k2,2n "${PROTEIN_TSS_BED}") \
  -D b > "${global_dir}/all_super_enhancer_gene.bed"

# The global catalog has four columns, so the TSS gene identifier is field 8.
awk 'NF>=8{print $8}' "${global_dir}/all_super_enhancer_gene.bed" | awk 'NF&&!seen[$1]++' \
  > "${global_dir}/gene_with_super_enhancer.txt"
awk 'NR==FNR{seen[$1]=1;next}!seen[$4]' "${global_dir}/gene_with_super_enhancer.txt" "${PROTEIN_TSS_BED}" \
  > "${global_dir}/non_super_enhancer_gene.bed"
awk 'NF>=4{print $4}' "${global_dir}/non_super_enhancer_gene.bed" | awk 'NF&&!seen[$1]++' \
  > "${global_dir}/gene_without_super_enhancer.txt"

join_by_first_field "${global_dir}/gene_with_super_enhancer.txt" "${EXPRESSION_TPM}" "${global_dir}/.with_expr.tsv"
join_by_first_field "${global_dir}/gene_without_super_enhancer.txt" "${EXPRESSION_TPM}" "${global_dir}/.without_expr.tsv"
add_expression_summary "${global_dir}/.with_expr.tsv" "WithSuperEnhancer" "${global_dir}/output2.txt"
add_expression_summary "${global_dir}/.without_expr.tsv" "WithoutSuperEnhancer" "${global_dir}/output1.txt"
cat "${global_dir}/output1.txt" "${global_dir}/output2.txt" > "${global_dir}/AA_with_without_super_enhancer_expression.txt"

join_by_first_field "${global_dir}/gene_with_super_enhancer.txt" "${EXPRESSION_TAU}" "${global_dir}/2_tau.txt"
join_by_first_field "${global_dir}/gene_without_super_enhancer.txt" "${EXPRESSION_TAU}" "${global_dir}/1_tau.txt"
join_by_first_field "${global_dir}/output1.txt" "${global_dir}/1_tau.txt" "${global_dir}/output1_tau.txt"
join_by_first_field "${global_dir}/output2.txt" "${global_dir}/2_tau.txt" "${global_dir}/output2_tau.txt"
cat "${global_dir}/output1_tau.txt" "${global_dir}/output2_tau.txt" > "${global_dir}/AA_with_without_super_enhancer_tau.txt"
join_by_first_field "${global_dir}/output1_tau.txt" "${SHEEP_HUMAN_CONSERVATION}" "${global_dir}/output1_human.txt"
join_by_first_field "${global_dir}/output2_tau.txt" "${SHEEP_HUMAN_CONSERVATION}" "${global_dir}/output2_human.txt"
cat "${global_dir}/output1_human.txt" "${global_dir}/output2_human.txt" \
  > "${global_dir}/AA_with_without_super_enhancer_compare_to_human.txt"

log "Target-gene and expression comparison tables completed."
