#!/usr/bin/env bash

set -euo pipefail

# Figure 2c: summarize tissue-matched DNA methylation across E1-E11.

RESULTS_BASE="/vol2/zhangshiwen/rrbs/wgbs_bismark/results"
SAMPLE_LIST="/vol2/zhangshiwen/rrbs/wgbs_bismark/samples.tsv"
META_TSV="/vol2/zhangshiwen/rrbs/md5sum.txt"
STATE_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/All_chromatin_state"
OUTDIR="/vol2/zhangshiwen/rrbs/wgbs_bismark/state_methylation_by_sample_tissue"

REUSE_VALID_RESULTS="${REUSE_VALID_RESULTS:-1}"
STATES=(E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11)

STATE_BED_DIR="${OUTDIR}/tmp/tissue_state_beds"
SAMPLE_TMP_DIR="${OUTDIR}/tmp/per_sample"
PER_SAMPLE_DIR="${OUTDIR}/per_sample"
LOG_DIR="${OUTDIR}/log"
PLOT_DIR="${OUTDIR}/plot"

mkdir -p \
  "${STATE_BED_DIR}" \
  "${SAMPLE_TMP_DIR}" \
  "${PER_SAMPLE_DIR}" \
  "${LOG_DIR}" \
  "${PLOT_DIR}"

command -v awk >/dev/null 2>&1 || {
  echo "[ERROR] awk was not found." >&2
  exit 1
}
command -v bedtools >/dev/null 2>&1 || {
  echo "[ERROR] bedtools was not found." >&2
  exit 1
}
command -v zcat >/dev/null 2>&1 || {
  echo "[ERROR] zcat was not found." >&2
  exit 1
}

[[ -d "${RESULTS_BASE}" ]] || {
  echo "[ERROR] Results directory was not found: ${RESULTS_BASE}" >&2
  exit 1
}
[[ -s "${SAMPLE_LIST}" ]] || {
  echo "[ERROR] Sample table was not found or is empty: ${SAMPLE_LIST}" >&2
  exit 1
}
[[ -s "${META_TSV}" ]] || {
  echo "[ERROR] Metadata table was not found or is empty: ${META_TSV}" >&2
  exit 1
}
[[ -d "${STATE_DIR}" ]] || {
  echo "[ERROR] Chromatin-state directory was not found: ${STATE_DIR}" >&2
  exit 1
}

: > "${LOG_DIR}/reused_valid_results.log"
: > "${LOG_DIR}/processed_results.log"
: > "${LOG_DIR}/missing_cov.log"
: > "${LOG_DIR}/missing_tissue_in_metadata.log"
: > "${LOG_DIR}/missing_state_bed.log"
: > "${LOG_DIR}/no_state_overlap.log"
: > "${LOG_DIR}/missing_result_when_merging.log"

is_valid_per_sample_result() {
  local result_file="$1"

  [[ -s "${result_file}" ]] || return 1

  awk -F'\t' '
    NF != 9 {
      exit 1
    }
    $4 !~ /^E([1-9]|10|11)$/ {
      exit 1
    }
    {
      count[$4]++
    }
    END {
      for (i = 1; i <= 11; i++) {
        state = "E" i
        if (count[state] != 1) {
          exit 1
        }
      }
    }
  ' "${result_file}"
}

resolve_state_bed() {
  local tissue="$1"
  local candidate
  local candidates=(
    "${tissue}"
    "${tissue//_/-}"
    "${tissue// /-}"
    "${tissue//-/_}"
  )

  for candidate in "${candidates[@]}"; do
    if [[ -s "${STATE_BED_DIR}/${candidate}.states.sorted.bed" ]]; then
      printf '%s\t%s\n' \
        "${candidate}" \
        "${STATE_BED_DIR}/${candidate}.states.sorted.bed"
      return 0
    fi
  done

  return 1
}

echo "[INFO] Building tissue-level E1-E11 chromatin-state BED files."

shopt -s nullglob
for old_file in \
  "${STATE_BED_DIR}"/*.unsorted.bed \
  "${STATE_BED_DIR}"/*.states.sorted.bed; do
  rm -f -- "${old_file}"
done

state_file_count=0
for state in "${STATES[@]}"; do
  state_files=("${STATE_DIR}"/*_"${state}".bed)

  for state_file in "${state_files[@]}"; do
    base_name="$(basename "${state_file}")"
    tissue="${base_name%_"${state}".bed}"

    awk -v OFS='\t' -v state="${state}" '
      BEGIN {
        FS = OFS = "\t"
      }
      $0 !~ /^#/ &&
      NF >= 3 &&
      $2 ~ /^[0-9]+$/ &&
      $3 ~ /^[0-9]+$/ {
        start = $2
        end = $3
        if (start > end) {
          temporary = start
          start = end
          end = temporary
        }
        if (end > start) {
          print $1, start, end, state
        }
      }
    ' "${state_file}" >> "${STATE_BED_DIR}/${tissue}.unsorted.bed"

    state_file_count=$((state_file_count + 1))
  done
done

if (( state_file_count == 0 )); then
  echo "[ERROR] No files matching <tissue>_E1.bed through <tissue>_E11.bed were found in ${STATE_DIR}." >&2
  exit 1
fi

combined_bed_count=0
for unsorted_bed in "${STATE_BED_DIR}"/*.unsorted.bed; do
  sorted_bed="${unsorted_bed%.unsorted.bed}.states.sorted.bed"
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "${unsorted_bed}" > "${sorted_bed}"
  rm -f -- "${unsorted_bed}"
  combined_bed_count=$((combined_bed_count + 1))
done

if (( combined_bed_count == 0 )); then
  echo "[ERROR] No tissue-level chromatin-state BED files were created." >&2
  exit 1
fi

mapfile -t SAMPLES < <(
  awk -F'\t' '
    NR > 1 && $1 != "" {
      sub(/\r$/, "", $1)
      print $1
    }
  ' "${SAMPLE_LIST}" | LC_ALL=C sort -u
)

if (( ${#SAMPLES[@]} == 0 )); then
  echo "[ERROR] No sample IDs were found in the first column of ${SAMPLE_LIST}." >&2
  exit 1
fi

declare -A SAMPLE_TO_TISSUE
while IFS=$'\t' read -r run_accession fastq_aspera fastq_md5 tissue_type extra_fields; do
  run_accession="${run_accession%$'\r'}"
  tissue_type="${tissue_type%$'\r'}"

  [[ "${run_accession}" == "run_accession" ]] && continue
  [[ -z "${run_accession}" ]] && continue
  [[ -z "${tissue_type}" ]] && continue

  SAMPLE_TO_TISSUE["${run_accession}"]="${tissue_type}"
done < "${META_TSV}"

echo "[INFO] Samples to process: ${#SAMPLES[@]}"

for sample in "${SAMPLES[@]}"; do
  per_sample_output="${PER_SAMPLE_DIR}/${sample}.sample_state_methylation.tsv"

  if [[ "${REUSE_VALID_RESULTS}" == "1" ]] &&
     is_valid_per_sample_result "${per_sample_output}"; then
    printf '%s\t%s\n' "${sample}" "${per_sample_output}" \
      >> "${LOG_DIR}/reused_valid_results.log"
    echo "[INFO] Reusing valid E1-E11 result: ${sample}"
    continue
  fi

  sample_tissue="${SAMPLE_TO_TISSUE[${sample}]:-}"
  if [[ -z "${sample_tissue}" ]]; then
    printf '%s\tNO_TISSUE_IN_METADATA\n' "${sample}" \
      >> "${LOG_DIR}/missing_tissue_in_metadata.log"
    echo "[WARN] No tissue metadata was found for ${sample}." >&2
    continue
  fi

  state_resolution="$(resolve_state_bed "${sample_tissue}" || true)"
  if [[ -z "${state_resolution}" ]]; then
    printf '%s\t%s\tNO_MATCHED_STATE_BED\n' "${sample}" "${sample_tissue}" \
      >> "${LOG_DIR}/missing_state_bed.log"
    echo "[WARN] No tissue-matched E1-E11 state BED was found for ${sample} (${sample_tissue})." >&2
    continue
  fi

  IFS=$'\t' read -r state_tissue state_bed <<< "${state_resolution}"

  mapfile -t cov_files < <(
    find "${RESULTS_BASE}/${sample}" \
      -maxdepth 3 \
      -type f \
      -name "*.deduplicated.bismark.cov.gz" \
      -print 2>/dev/null | LC_ALL=C sort
  )

  if (( ${#cov_files[@]} == 0 )); then
    mapfile -t cov_files < <(
      find "${RESULTS_BASE}/${sample}" \
        -maxdepth 3 \
        -type f \
        -name "*.bismark.cov.gz" \
        -print 2>/dev/null | LC_ALL=C sort
    )
  fi

  if (( ${#cov_files[@]} == 0 )); then
    printf '%s\t%s\tNO_COV_FILE\n' "${sample}" "${sample_tissue}" \
      >> "${LOG_DIR}/missing_cov.log"
    echo "[WARN] No Bismark coverage file was found for ${sample}." >&2
    continue
  fi

  if (( ${#cov_files[@]} > 1 )); then
    echo "[WARN] Multiple coverage files were found for ${sample}; using ${cov_files[0]}." >&2
  fi
  cov_file="${cov_files[0]}"

  cov_bed="${SAMPLE_TMP_DIR}/${sample}.cov.bed"
  state_summary="${SAMPLE_TMP_DIR}/${sample}.state_summary.tsv"

  zcat "${cov_file}" |
    awk -v OFS='\t' '
      BEGIN {
        FS = OFS = "\t"
      }
      NF >= 6 &&
      $2 ~ /^[0-9]+$/ &&
      $3 ~ /^[0-9]+$/ &&
      $5 ~ /^[0-9]+([.][0-9]+)?$/ &&
      $6 ~ /^[0-9]+([.][0-9]+)?$/ {
        start = $2 - 1
        if (start < 0) {
          start = 0
        }
        print $1, start, $3, $4, $5, $6
      }
    ' |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "${cov_bed}"

  bedtools intersect \
    -sorted \
    -a "${cov_bed}" \
    -b "${state_bed}" \
    -wa \
    -wb |
    awk -v OFS='\t' '
      BEGIN {
        FS = OFS = "\t"
      }
      {
        state = $10
        methylated[state] += $5
        unmethylated[state] += $6
        cpg_sites[state]++
      }
      END {
        for (state in methylated) {
          total = methylated[state] + unmethylated[state]
          if (total > 0) {
            weighted_percentage = 100 * methylated[state] / total
          } else {
            weighted_percentage = "NA"
          }
          print state, methylated[state], unmethylated[state], total, cpg_sites[state], weighted_percentage
        }
      }
    ' |
    LC_ALL=C sort -k1,1V > "${state_summary}"

  if [[ ! -s "${state_summary}" ]]; then
    printf '%s\t%s\t%s\tNO_STATE_OVERLAP\n' \
      "${sample}" "${sample_tissue}" "${state_tissue}" \
      >> "${LOG_DIR}/no_state_overlap.log"
    echo "[WARN] No CpG sites overlapped chromatin states for ${sample}." >&2
  fi

  awk \
    -v OFS='\t' \
    -v sample="${sample}" \
    -v sample_tissue="${sample_tissue}" \
    -v state_tissue="${state_tissue}" '
      BEGIN {
        FS = OFS = "\t"
      }
      {
        observed[$1] = 1
        methylated[$1] = $2
        unmethylated[$1] = $3
        total[$1] = $4
        cpg_sites[$1] = $5
        weighted_percentage[$1] = $6
      }
      END {
        for (i = 1; i <= 11; i++) {
          state = "E" i
          if (state in observed) {
            print sample, sample_tissue, state_tissue, state, \
              methylated[state], unmethylated[state], total[state], \
              cpg_sites[state], weighted_percentage[state]
          } else {
            print sample, sample_tissue, state_tissue, state, \
              0, 0, 0, 0, "NA"
          }
        }
      }
    ' "${state_summary}" > "${per_sample_output}"

  if ! is_valid_per_sample_result "${per_sample_output}"; then
    echo "[ERROR] The generated result does not contain exactly one row for each state E1-E11: ${per_sample_output}" >&2
    exit 1
  fi

  printf '%s\t%s\t%s\t%s\n' \
    "${sample}" "${sample_tissue}" "${state_tissue}" "${per_sample_output}" \
    >> "${LOG_DIR}/processed_results.log"

  rm -f -- "${cov_bed}" "${state_summary}"
  echo "[INFO] Completed: ${sample}"
done

master_output="${OUTDIR}/all_samples_sample_state_weighted_methylation.tsv"
printf '%s\n' \
  $'sample\tsample_tissue\tstate_tissue\tstate\tmethylated_count\tunmethylated_count\ttotal_count\tcpg_sites\tweighted_meth_pct' \
  > "${master_output}"

merged_sample_count=0
for sample in "${SAMPLES[@]}"; do
  per_sample_output="${PER_SAMPLE_DIR}/${sample}.sample_state_methylation.tsv"

  if is_valid_per_sample_result "${per_sample_output}"; then
    cat "${per_sample_output}" >> "${master_output}"
    merged_sample_count=$((merged_sample_count + 1))
  else
    printf '%s\tNO_VALID_E1_E11_RESULT\n' "${sample}" \
      >> "${LOG_DIR}/missing_result_when_merging.log"
  fi
done

if (( merged_sample_count == 0 )); then
  echo "[ERROR] No valid sample-level E1-E11 results were available for merging." >&2
  exit 1
fi

echo "[INFO] Merged samples: ${merged_sample_count}"
echo "[INFO] Figure 2c input table: ${master_output}"

