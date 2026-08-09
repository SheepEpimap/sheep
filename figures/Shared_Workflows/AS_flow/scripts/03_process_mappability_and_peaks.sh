#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 12
#SBATCH --mem=96G
#SBATCH -t 5-00:00:00

set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_config.sh"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_lib.sh"

usage() {
  cat <<'EOF'
Usage:
  sbatch scripts/03_process_mappability_and_peaks.sh SAMPLE_ID BAM_SUFFIX
  bash   scripts/03_process_mappability_and_peaks.sh SAMPLE_ID BAM_SUFFIX

Examples:
  ... H3K27ac_liver_01 bowtie2.mapped.filtered.sort.bam
  ... ATAC_liver_01 last.bam

RNASeq samples are intentionally rejected. Use 05_process_rnaseq.sh instead.
EOF
}

sample_id="${1:-${sample_id:-}}"
bam_suffix="${2:-${BAM_SUFFIX:-}}"
[[ -n "${sample_id}" && -n "${bam_suffix}" ]] || { usage; exit 2; }
as_sample_is_rnaseq "${sample_id}" && as_die "RNASeq does not use MACS2 peak calling in this package; run 05_process_rnaseq.sh"

as_activate_conda "${ATAC_ENV}"
as_require_cmds samtools bedtools macs2 Rscript gzip awk
as_require_file "${REFERENCE_FASTA}" "reference FASTA"
as_require_file "${CHROM_SIZES}" "chromosome-size file"

input_bam="${INPUT_BAM_DIR}/${sample_id}.${bam_suffix}"
as_require_file "${input_bam}" "input BAM"

for k in 36 40 72; do
  as_require_file "${MAPPABILITY_RESOURCE_DIR}/k${k}/Ovis_aries_v2.0.k${k}.bed" "k=${k} mappability BED"
done

as_info "Sample: ${sample_id}"
as_info "Input BAM: ${input_bam}"

for k in 40 36 72; do
  map_bed="${MAPPABILITY_RESOURCE_DIR}/k${k}/Ovis_aries_v2.0.k${k}.bed"
  out_bam="${MAPPABILITY_DIR}/${sample_id}.${k}.bam"
  out_tag="${MAPPABILITY_DIR}/${sample_id}.${k}.tagAlign.gz"
  stat_file="${PEAK_METRICS_DIR}/${sample_id}_${k}.mappability_read_count.txt"

  if [[ ! -s "${out_bam}" ]]; then
    as_info "Applying k=${k} mappability filter"
    bedtools intersect -abam "${input_bam}" -b "${map_bed}" -wa > "${out_bam}"
    samtools index -f "${out_bam}"
  else
    as_info "Filtered BAM exists; skipping k=${k}: ${out_bam}"
  fi

  if [[ ! -s "${out_tag}" ]]; then
    bedtools bamtobed -i "${out_bam}" \
      | awk 'BEGIN{OFS="\t"} {$4="N"; $5="1000"; print}' \
      | gzip -nc > "${out_tag}"
  fi
  samtools view -c "${out_bam}" > "${stat_file}"
done

[[ "${sample_id}" == Input* ]] && {
  as_info "Input/control sample detected; peak calling skipped by design"
  exit 0
}

bam72="${MAPPABILITY_DIR}/${sample_id}.72.bam"
spp_stats="${PEAK_METRICS_DIR}/${sample_id}.spp_stats.txt"
spp_pdf="${PEAK_METRICS_DIR}/${sample_id}.Cross_Correlation.pdf"
fraglen=100

if [[ -f "${SPP_RUNNER}" ]]; then
  as_info "Running SPP cross-correlation"
  export R_LIBS_USER="${SPP_R_LIB}"
  if Rscript "${SPP_RUNNER}" \
      -c="${bam72}" \
      -rf \
      -out="${spp_stats}" \
      -p=8 \
      -s=0:2:400 \
      -savp="${spp_pdf}" \
      -tmpdir="${PEAK_TEMP_DIR}"; then
    candidate="$(awk -v target="${sample_id}.72.bam" '
      $1 == target {
        split($3, a, ",");
        print a[1];
        exit
      }' "${spp_stats}")"
    if [[ "${candidate}" =~ ^[0-9]+$ ]] && (( candidate >= 20 && candidate <= 1000 )); then
      fraglen="${candidate}"
    else
      as_warn "SPP fragment length unavailable/invalid; using ${fraglen} bp"
    fi
  else
    as_warn "SPP failed; using ${fraglen} bp"
    : > "${spp_stats}"
  fi
  unset R_LIBS_USER
else
  as_warn "SPP runner not found; using ${fraglen} bp: ${SPP_RUNNER}"
fi

if as_sample_is_atac "${sample_id}"; then
  qvalue="${ATAC_QVALUE}"
else
  qvalue="${OTHER_PEAK_QVALUE}"
fi
as_info "MACS2 q-value cutoff for ${sample_id}: ${qvalue}"

peak_file="${PEAK_DIR}/${sample_id}_peaks.broadPeak"
final_peak="${PEAK_DIR}/${sample_id}_Peaks.bed"

if [[ ! -s "${final_peak}" ]]; then
  macs2 callpeak \
    -t "${bam72}" \
    -n "${sample_id}" \
    --outdir "${PEAK_DIR}" \
    -g "${GENOME_SIZE}" \
    -q "${qvalue}" \
    -f BAMPE \
    --nomodel \
    --extsize "${fraglen}" \
    --shift -100 \
    --keep-dup all \
    --broad

  [[ -s "${peak_file}" ]] || as_die "MACS2 did not create ${peak_file}"
  mv -f "${peak_file}" "${final_peak}"
else
  as_info "Peak file exists; skipping MACS2: ${final_peak}"
fi

wc -l < "${final_peak}" > "${PEAK_METRICS_DIR}/${sample_id}_peak_number.txt"
bedtools genomecov -i "${final_peak}" -g "${CHROM_SIZES}" \
  | tail -n 1 \
  | awk -v sample="${sample_id}" '{print sample "\t" $5}' \
  > "${PEAK_METRICS_DIR}/${sample_id}_peak_coverage.txt"

as_info "Completed mappability filtering and peak calling for ${sample_id}"
