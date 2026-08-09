#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 10-15:00:00
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
sample_id="${1:-}"; [[ -n "${sample_id}" ]] || as_die "Usage: 13_run_gatk_ase.sh SAMPLE_ID"
as_activate_conda "${ASE_ENV}"; as_require_cmds gatk samtools bgzip; as_require_file "${REFERENCE_FASTA}" "reference FASTA"
input_bam="${WASP_FINAL_DIR}/${sample_id}.keep.merged.sorted.bam"; as_require_file "${input_bam}" "WASP-corrected BAM"; [[ -s "${input_bam}.bai" ]] || samtools index "${input_bam}"
ref_dict="${REFERENCE_FASTA%.*}.dict"; [[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "${REFERENCE_FASTA}"; [[ -s "${ref_dict}" ]] || as_die "Missing GATK dictionary ${ref_dict}; create with gatk CreateSequenceDictionary"
shopt -s nullglob; vcfs=("${PHASING_DIR}"/phased_*.vcf.gz); (( ${#vcfs[@]} > 0 )) || as_die "No phased VCF files"; vcf_args=()
for vcf in "${vcfs[@]}"; do [[ -s "${vcf}.tbi" || -s "${vcf}.csi" ]] || as_die "VCF index missing: ${vcf}"; bgzip -t "${vcf}"; vcf_args+=( -V "${vcf}" ); done
output="${ASE_TABLE_DIR}/${sample_id}.output.table"; [[ -s "${output}" ]] && { as_info "Output exists; skipping"; exit 0; }
gatk --java-options '-Xmx56g' ASEReadCounter -R "${REFERENCE_FASTA}" -I "${input_bam}" "${vcf_args[@]}" -O "${output}"; [[ -s "${output}" ]] || as_die "GATK output missing"
