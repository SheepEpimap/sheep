#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 10-15:00:00
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
sample_id="${1:-}"; input_override="${2:-}"; [[ -n "${sample_id}" ]] || as_die "Usage: 10_run_wasp.sh SAMPLE_ID [INPUT_BAM]"
as_activate_conda "${WASP_ENV}"; as_require_cmds python samtools
as_require_file "${REFERENCE_FASTA}" "reference FASTA"; as_require_file "${SAMPLE_MAP}" "sample map"
for f in "${WASP_HAPS_H5}" "${WASP_SNP_TAB_H5}" "${WASP_SNP_INDEX_H5}" "${WASP_PATH}/mapping/find_intersecting_snps.py" "${WASP_PATH}/mapping/filter_remapped_reads.py"; do as_require_file "${f}" "WASP resource"; done
individual="$(awk -v s="${sample_id}" '$1==s {print $2; exit}' "${SAMPLE_MAP}" | tr -d '\r')"; [[ -n "${individual}" ]] || as_die "No individual ID for ${sample_id}"
if [[ -n "${input_override}" ]]; then input_bam="${input_override}"; elif as_sample_is_rnaseq "${sample_id}"; then input_bam="${RNASEQ_ALIGN_DIR}/${sample_id}.filtered.bam"; else input_bam="${MAPPABILITY_DIR}/${sample_id}.72.bam"; fi
as_require_file "${input_bam}" "WASP input BAM"
sample_find_dir="${WASP_FIND_DIR}/${sample_id}"; sample_remap_dir="${WASP_REMAP_DIR}/${sample_id}"; mkdir -p "${sample_find_dir}" "${sample_remap_dir}" "${WASP_FINAL_DIR}"
exec > >(tee -a "${LOG_DIR}/10_wasp_${sample_id}.log") 2>&1
rg_bam="${WASP_TEMP_DIR}/${sample_id}.RG.sorted.bam"
if [[ ! -s "${rg_bam}" ]]; then
  samtools addreplacerg -@ "${WASP_THREADS}" -r "ID:${sample_id}" -r "SM:${individual}" -r "LB:${individual}" -r "PL:ILLUMINA" -o "${WASP_TEMP_DIR}/${sample_id}.RG.bam" "${input_bam}"
  samtools sort -@ "${WASP_THREADS}" -o "${rg_bam}" "${WASP_TEMP_DIR}/${sample_id}.RG.bam"; samtools index -f "${rg_bam}"; rm -f "${WASP_TEMP_DIR}/${sample_id}.RG.bam"
fi
paired_flag=(); if (( $(samtools view -c -f 1 "${rg_bam}") > 0 )); then paired_flag=(--is_paired_end); fi
python "${WASP_PATH}/mapping/find_intersecting_snps.py" --is_sorted "${paired_flag[@]}" --output_dir "${sample_find_dir}" --snp_tab "${WASP_SNP_TAB_H5}" --snp_index "${WASP_SNP_INDEX_H5}" --haplotype "${WASP_HAPS_H5}" --samples "${individual}" "${rg_bam}"
base="$(basename "${rg_bam}" .bam)"; keep_bam="${sample_find_dir}/${base}.keep.bam"; to_remap_bam="${sample_find_dir}/${base}.to.remap.bam"; as_require_file "${keep_bam}" "WASP keep BAM"
final_bam="${WASP_FINAL_DIR}/${sample_id}.keep.merged.sorted.bam"; remap_keep="${sample_remap_dir}/${sample_id}.remap.keep.bam"
if (( ${#paired_flag[@]} > 0 )); then fq1="${sample_find_dir}/${base}.remap.fq1.gz"; fq2="${sample_find_dir}/${base}.remap.fq2.gz"; [[ -s "${fq1}" && -s "${fq2}" ]] && remap_has_reads=true || remap_has_reads=false; else fq="${sample_find_dir}/${base}.remap.fq.gz"; [[ -s "${fq}" ]] && remap_has_reads=true || remap_has_reads=false; fi
if [[ "${remap_has_reads}" == true ]]; then
  remap_bam="${sample_remap_dir}/${sample_id}.remap.sorted.bam"
  if as_sample_is_rnaseq "${sample_id}"; then
    as_require_cmds STAR; as_require_dir "${STAR_INDEX}" "STAR index"; prefix="${sample_remap_dir}/${sample_id}.STAR."
    if (( ${#paired_flag[@]} > 0 )); then STAR --genomeDir "${STAR_INDEX}" --readFilesIn "${fq1}" "${fq2}" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN "${WASP_THREADS}" --outFileNamePrefix "${prefix}"; else STAR --genomeDir "${STAR_INDEX}" --readFilesIn "${fq}" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN "${WASP_THREADS}" --outFileNamePrefix "${prefix}"; fi
    mv -f "${prefix}Aligned.sortedByCoord.out.bam" "${remap_bam}"
  else
    as_require_cmds bwa; [[ -s "${REFERENCE_FASTA}.bwt" ]] || bwa index "${REFERENCE_FASTA}"
    if (( ${#paired_flag[@]} > 0 )); then bwa mem -M -t "${WASP_THREADS}" "${REFERENCE_FASTA}" "${fq1}" "${fq2}" | samtools view -b - | samtools sort -@ "${WASP_THREADS}" -o "${remap_bam}"; else bwa mem -M -t "${WASP_THREADS}" "${REFERENCE_FASTA}" "${fq}" | samtools view -b - | samtools sort -@ "${WASP_THREADS}" -o "${remap_bam}"; fi
  fi
  samtools index -f "${remap_bam}"; python "${WASP_PATH}/mapping/filter_remapped_reads.py" "${to_remap_bam}" "${remap_bam}" "${remap_keep}"; as_require_file "${remap_keep}" "filtered remapped BAM"
  merged="${WASP_TEMP_DIR}/${sample_id}.keep.merged.bam"; samtools merge -f "${merged}" "${keep_bam}" "${remap_keep}"; samtools sort -@ "${WASP_THREADS}" -o "${final_bam}" "${merged}"; rm -f "${merged}"
else
  as_warn "No reads require remapping; using keep set"; samtools sort -@ "${WASP_THREADS}" -o "${final_bam}" "${keep_bam}"
fi
samtools index -f "${final_bam}"; as_info "WASP-corrected BAM: ${final_bam}"
