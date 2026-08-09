#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00

set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_config.sh"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_lib.sh"

as_activate_conda "${ATAC_ENV}"
as_require_cmds samtools cut awk genmap
as_require_file "${REFERENCE_FASTA}" "reference FASTA"

as_info "Preparing chromosome sizes and GenMap mappability resources"

samtools faidx "${REFERENCE_FASTA}"

if [[ -f "${REFERENCE_CHROM_SIZES_SOURCE}" ]]; then
  cp -f "${REFERENCE_CHROM_SIZES_SOURCE}" "${CHROM_SIZES}"
else
  cut -f1,2 "${REFERENCE_FASTA}.fai" > "${CHROM_SIZES}"
fi

# Keep canonical chromosomes by default. Edit this regex in 00_config.sh workflow
# only if scaffolds/unplaced contigs are intentionally required.
awk '$1 ~ /^chr([0-9]+|X|Y|M|MT)$/ {print $1 "\t" $2}' "${CHROM_SIZES}" > "${CHROM_LIST}.tmp"
if [[ -s "${CHROM_LIST}.tmp" ]]; then
  mv -f "${CHROM_LIST}.tmp" "${CHROM_LIST}"
else
  as_warn "No chr-prefixed canonical names detected; retaining all reference sequences"
  cp -f "${CHROM_SIZES}" "${CHROM_LIST}"
  rm -f "${CHROM_LIST}.tmp"
fi

calculated_size="$(awk '{sum += $2} END {printf "%.0f", sum}' "${CHROM_SIZES}")"
as_info "Genome size from chromosome-size file: ${calculated_size}; configured MACS2 size: ${GENOME_SIZE}"

mkdir -p "${MAPPABILITY_RESOURCE_DIR}"
index_prefix="${MAPPABILITY_RESOURCE_DIR}/Ovis_aries_v2.0.genmap.index"

if [[ ! -d "${index_prefix}" ]]; then
  genmap index -F "${REFERENCE_FASTA}" -I "${index_prefix}"
else
  as_info "GenMap index exists; skipping: ${index_prefix}"
fi

for k in 36 40 72; do
  output_dir="${MAPPABILITY_RESOURCE_DIR}/k${k}"
  final_bed="${output_dir}/Ovis_aries_v2.0.k${k}.bed"
  mkdir -p "${output_dir}"
  if [[ -s "${final_bed}" ]]; then
    as_info "Mappability BED exists; skipping k=${k}: ${final_bed}"
    continue
  fi

  genmap map \
    -K "${k}" \
    -I "${index_prefix}" \
    -O "${output_dir}" \
    --fast \
    --min-mapq 1 \
    --bed

  find "${output_dir}" -maxdepth 1 -type f -name '*.bed' ! -name "$(basename "${final_bed}")" -print0 \
    | sort -z \
    | xargs -0 cat > "${final_bed}"

  [[ -s "${final_bed}" ]] || as_die "Failed to create mappability BED for k=${k}"
  as_info "Created ${final_bed}"
done

as_info "Reference and mappability preparation completed"
