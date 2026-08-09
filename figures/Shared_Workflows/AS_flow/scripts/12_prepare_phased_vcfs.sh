#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
as_activate_conda "${ASE_ENV}"; as_require_cmds bgzip tabix bcftools
shopt -s nullglob; vcfs=("${PHASING_DIR}"/phased_*.vcf.gz); (( ${#vcfs[@]} > 0 )) || as_die "No phased VCF files"
for vcf in "${vcfs[@]}"; do bgzip -t "${vcf}" || as_die "Invalid bgzip VCF: ${vcf}"; tabix -f -p vcf "${vcf}"; bcftools view -h "${vcf}" >/dev/null; done
