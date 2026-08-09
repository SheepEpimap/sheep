#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/00_lib.sh"
as_activate_conda "${WASP_ENV}"
as_require_cmds samtools bcftools awk
as_require_file "${REFERENCE_FASTA}" "reference FASTA"
as_require_file "${PHASING_VCF}" "unphased VCF"
# Derived resources stay inside this package. External defaults are replaceable in 00_config.sh.
if [[ -s "${REFERENCE_CHROM_SIZES_SOURCE}" ]]; then
  cp -f "${REFERENCE_CHROM_SIZES_SOURCE}" "${CHROM_SIZES}"
else
  [[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "${REFERENCE_FASTA}"
  cut -f1,2 "${REFERENCE_FASTA}.fai" > "${CHROM_SIZES}"
fi
awk '$1 ~ /^chr([0-9]+|X|Y|M|MT)$/ {print $1"\t"$2}' "${CHROM_SIZES}" > "${CHROM_LIST}.tmp"
if [[ -s "${CHROM_LIST}.tmp" ]]; then mv -f "${CHROM_LIST}.tmp" "${CHROM_LIST}"; else cp -f "${CHROM_SIZES}" "${CHROM_LIST}"; rm -f "${CHROM_LIST}.tmp"; fi
mkdir -p "${PHASING_MAP_DIR}"
while IFS=$'\t' read -r chrom _length; do
  [[ -n "${chrom}" ]] || continue
  out="${PHASING_MAP_DIR}/beagle_${chrom}.map"
  as_info "Preparing Beagle map: ${chrom}"
  # Fallback retained from the source workflow: physical position/1e6 is used as approximate cM.
  # Replace these maps with an authoritative recombination map when available.
  bcftools query -r "${chrom}" -f '%CHROM\t%POS\t%ID\n' "${PHASING_VCF}" \
    | awk -F'\t' -v OFS='\t' '{id=($3=="." || $3=="" ? $1"_"$2 : $3); print $1,id,$2/1000000,$2}' > "${out}"
  [[ -s "${out}" ]] || as_warn "No VCF records found for ${chrom}; map is empty"
done < "${CHROM_LIST}"
as_info "Phasing resources prepared under ${RESOURCE_DIR}"
