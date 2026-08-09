#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Core processing order. Signal quantification and motif enrichment are included
# because they are non-visual analyses; set the variables below to 0 to skip them.
RUN_SIGNAL_QUANTIFICATION="${RUN_SIGNAL_QUANTIFICATION:-1}"
RUN_MOTIF_ENRICHMENT="${RUN_MOTIF_ENRICHMENT:-1}"

steps=(
  02_initialize_and_prepare_inputs.sh
  03_build_and_install_rose_annotation.sh
  04_link_h3k27ac_bams.sh
  05_run_rose.sh
  06_merge_replicates_and_build_catalog.sh
  07_genomic_annotation_summary.sh
  08_compare_normal_and_super_enhancers.sh
  09_build_enhancer_presence_matrices.sh
  11_assign_target_genes_and_expression.sh
)

for step in "${steps[@]}"; do
  printf '[pipeline] running %s\n' "${step}" >&2
  bash "${SCRIPT_DIR}/${step}"
done

if [[ "${RUN_SIGNAL_QUANTIFICATION}" == 1 ]]; then
  bash "${SCRIPT_DIR}/10_quantify_h3k27ac_signal.sh"
fi
if [[ "${RUN_MOTIF_ENRICHMENT}" == 1 ]]; then
  bash "${SCRIPT_DIR}/12_run_motif_enrichment.sh"
fi

printf '[pipeline] all requested data-processing steps completed\n' >&2
