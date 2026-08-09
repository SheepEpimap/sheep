#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 7-00:00:00

set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_config.sh"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_lib.sh"

usage() {
  echo "Usage: sbatch scripts/05_process_rnaseq.sh RNASeq_SAMPLE_ID" >&2
}

sample_id="${1:-${sample_id:-}}"
[[ -n "${sample_id}" ]] || { usage; exit 2; }
as_sample_is_rnaseq "${sample_id}" || as_die "RNA-seq sample ID must start with RNASeq"

as_activate_conda "${RNASEQ_ENV}"
as_require_cmds STAR samtools salmon trim_galore stringtie
as_require_file "${REFERENCE_FASTA}" "reference FASTA"
as_require_file "${REFERENCE_GTF}" "reference GTF"
as_require_dir "${STAR_INDEX}" "STAR index"
as_require_dir "${SALMON_INDEX}" "Salmon index"

step_log="${RNASEQ_METRICS_DIR}/${sample_id}_step_log.txt"
{
  echo "Sample: ${sample_id}"
  echo "Started: $(date)"
  echo "Raw reads directory: ${RNASEQ_RAW_READS_DIR}"
  echo "Output root: ${RNASEQ_WORK_DIR}"
  echo "Important: this RNA-seq branch does not call peaks."
} > "${step_log}"

raw_r1="${RNASEQ_RAW_READS_DIR}/${sample_id}_R1.fq.gz"
raw_r2="${RNASEQ_RAW_READS_DIR}/${sample_id}_R2.fq.gz"
trimmed_r1="${RNASEQ_TRIM_DIR}/${sample_id}_R1_val_1.fq.gz"
trimmed_r2="${RNASEQ_TRIM_DIR}/${sample_id}_R2_val_2.fq.gz"

if [[ ! -s "${trimmed_r1}" || ! -s "${trimmed_r2}" ]]; then
  as_require_file "${raw_r1}" "RNA-seq R1"
  as_require_file "${raw_r2}" "RNA-seq R2"
  as_info "Running Trim Galore for ${sample_id}"
  trim_galore \
    -q "${MAPQ}" \
    --cores "${RNASEQ_THREADS}" \
    --paired \
    "${raw_r1}" "${raw_r2}" \
    -o "${RNASEQ_TRIM_DIR}" \
    >> "${step_log}" 2>&1
else
  as_info "Trimmed reads exist; skipping Trim Galore"
fi

aligned_bam="${RNASEQ_ALIGN_DIR}/${sample_id}.aligned.bam"
filtered_bam="${RNASEQ_ALIGN_DIR}/${sample_id}.filtered.bam"
star_prefix="${RNASEQ_ALIGN_DIR}/${sample_id}."
star_bam="${star_prefix}Aligned.sortedByCoord.out.bam"
star_log="${RNASEQ_TEMP_DIR}/star_${sample_id}.log"

if [[ ! -s "${aligned_bam}" ]]; then
  if [[ -s "${star_bam}" ]]; then
    mv -f "${star_bam}" "${aligned_bam}"
  else
    as_info "Running STAR for ${sample_id}"
    STAR \
      --genomeDir "${STAR_INDEX}" \
      --readFilesIn "${trimmed_r1}" "${trimmed_r2}" \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --runThreadN "${RNASEQ_THREADS}" \
      --outFileNamePrefix "${star_prefix}" \
      --outSAMattrRGline "ID:${sample_id}" "SM:${sample_id}" "PL:ILLUMINA" \
      --quantMode GeneCounts \
      --sjdbGTFfile "${REFERENCE_GTF}" \
      > "${star_log}" 2>&1
    [[ -s "${star_bam}" ]] || as_die "STAR did not create ${star_bam}; see ${star_log}"
    mv -f "${star_bam}" "${aligned_bam}"
  fi
else
  as_info "Aligned BAM exists; skipping STAR"
fi

[[ -s "${aligned_bam}.bai" ]] || samtools index -@ "${RNASEQ_THREADS}" "${aligned_bam}"

if [[ ! -s "${filtered_bam}" ]]; then
  as_info "Filtering RNA-seq BAM at MAPQ >= ${MAPQ}"
  samtools view -@ "${RNASEQ_THREADS}" -b -q "${MAPQ}" "${aligned_bam}" -o "${filtered_bam}"
  samtools index -@ "${RNASEQ_THREADS}" "${filtered_bam}"
else
  as_info "Filtered BAM exists; skipping BAM filtering"
fi

gene_gtf="${RNASEQ_QUANT_DIR}/${sample_id}.gtf"
gene_abundance="${RNASEQ_QUANT_DIR}/${sample_id}_gene_abundance.tsv"
if [[ ! -s "${gene_abundance}" ]]; then
  as_info "Running StringTie"
  stringtie \
    -p "${RNASEQ_THREADS}" \
    -G "${REFERENCE_GTF}" \
    -o "${gene_gtf}" \
    -A "${gene_abundance}" \
    "${filtered_bam}" \
    > "${RNASEQ_TEMP_DIR}/stringtie_${sample_id}.log" 2>&1
else
  as_info "StringTie result exists; skipping"
fi

salmon_out="${RNASEQ_QUANT_DIR}/${sample_id}_salmon"
if [[ ! -s "${salmon_out}/quant.sf" ]]; then
  as_info "Running Salmon"
  salmon quant \
    -i "${SALMON_INDEX}" \
    -l A \
    -1 "${trimmed_r1}" \
    -2 "${trimmed_r2}" \
    -p "${RNASEQ_THREADS}" \
    -o "${salmon_out}" \
    --gcBias \
    > "${RNASEQ_TEMP_DIR}/salmon_${sample_id}.log" 2>&1
else
  as_info "Salmon result exists; skipping"
fi

{
  echo "Completed: $(date)"
  echo "Filtered BAM: ${filtered_bam}"
  echo "StringTie abundance: ${gene_abundance}"
  echo "Salmon output: ${salmon_out}"
} >> "${step_log}"

as_info "RNA-seq processing completed for ${sample_id}"
