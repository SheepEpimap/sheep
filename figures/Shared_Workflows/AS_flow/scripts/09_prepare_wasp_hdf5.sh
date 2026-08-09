#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
as_activate_conda "${WASP_ENV}"
as_require_file "${CHROM_SIZES}" "chromosome-size file"
snp2h5="${WASP_PATH}/snp2h5/snp2h5"; as_require_file "${snp2h5}" "compiled WASP snp2h5 executable"
[[ -x "${snp2h5}" ]] || as_die "snp2h5 is not executable; compile it under ${WASP_PATH}/snp2h5"
shopt -s nullglob; vcfs=("${PHASING_DIR}"/phased_*.vcf.gz); (( ${#vcfs[@]} > 0 )) || as_die "No phased VCF files found"
mkdir -p "${WASP_HDF5_DIR}"
if [[ -s "${WASP_HAPS_H5}" && -s "${WASP_SNP_TAB_H5}" && -s "${WASP_SNP_INDEX_H5}" ]]; then as_info "WASP HDF5 files exist; skipping"; exit 0; fi
rm -f "${WASP_HAPS_H5}" "${WASP_SNP_TAB_H5}" "${WASP_SNP_INDEX_H5}"
"${snp2h5}" --chrom "${CHROM_SIZES}" --format vcf --haplotype "${WASP_HAPS_H5}" --snp_index "${WASP_SNP_INDEX_H5}" --snp_tab "${WASP_SNP_TAB_H5}" "${vcfs[@]}"
for f in "${WASP_HAPS_H5}" "${WASP_SNP_TAB_H5}" "${WASP_SNP_INDEX_H5}"; do [[ -s "${f}" ]] || as_die "Missing HDF5 output: ${f}"; done
