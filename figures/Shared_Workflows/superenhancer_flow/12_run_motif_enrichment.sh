#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd findMotifsGenome.pl
require_file "${GENOME_FASTA}"
mkdir -p "${RESULT_DIR}/AA_Motif"

# Per-tissue super-enhancer motif enrichment. HOMER output includes its own HTML
# summaries, but this script does not create any custom figure or plotting code.
for tissue in "${TISSUES[@]}"; do
  bed="${RESULT_DIR}/${tissue}_super_enhancer.bed"
  require_file "${bed}"
  out="${RESULT_DIR}/AA_Motif/${tissue}_super_enhancer_motif"
  log "HOMER motif enrichment: ${tissue}"
  findMotifsGenome.pl "${bed}" "${GENOME_FASTA}" "${out}" \
    -len 8,10,12 -size 200 -mask -p "${HOMER_THREADS}"
done

# Optional cluster-specific analysis is executed only when cluster BED-like files
# already exist. Cluster generation was embedded in visualization code in the
# source Markdown and is therefore intentionally not recreated here.
cluster_dir="${RESULT_DIR}/Target_gene_cloest/AA_cluster"
shopt -s nullglob
cluster_files=("${cluster_dir}"/*_cluster.txt)
for cluster in "${cluster_files[@]}"; do
  id_bed="${cluster%.txt}_id.txt"
  awk 'BEGIN{OFS="\t"} NF>=3{id=$1":"$2"-"$3;print $1,$2,$3,id}' "${cluster}" > "${id_bed}"
  out="${cluster_dir}/AA_Motif/$(basename "${cluster}" .txt)_motif"
  findMotifsGenome.pl "${id_bed}" "${GENOME_FASTA}" "${out}" \
    -len 8,10,12 -size 200 -mask -p "${HOMER_THREADS}"
done

log "Motif enrichment completed."
