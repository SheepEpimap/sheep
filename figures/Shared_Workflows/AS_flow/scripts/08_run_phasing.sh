#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 28-00:00:00
#SBATCH -J as_phasing
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${SCRIPT_DIR}/00_config.sh"
source "${SCRIPT_DIR}/00_lib.sh"
as_activate_conda "${WASP_ENV}"
as_require_cmds beagle bcftools
as_require_file "${PHASING_VCF}" "unphased VCF"
as_require_file "${CHROM_LIST}" "chromosome list"
while IFS=$'\t' read -r chrom _length; do
  [[ -n "${chrom}" ]] || continue
  map_file="${PHASING_MAP_DIR}/beagle_${chrom}.map"; out_prefix="${PHASING_DIR}/phased_${chrom}"; out_vcf="${out_prefix}.vcf.gz"
  as_require_file "${map_file}" "Beagle map for ${chrom}"
  if [[ -s "${out_vcf}" ]]; then as_info "Phased VCF exists; skipping ${chrom}"; else
    beagle -Xmx140g gt="${PHASING_VCF}" out="${out_prefix}" chrom="${chrom}" nthreads="${PHASING_THREADS}" impute=false window=100 map="${map_file}"
  fi
  bcftools index -f -t "${out_vcf}"
done < "${CHROM_LIST}"
as_info "Phasing completed: ${PHASING_DIR}"
