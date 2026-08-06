#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 11: validate Bismark outputs and the E1-E11 panel table.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
SAMPLESHEET="${SAMPLESHEET:-${BASE}/wgbs_bismark/samples.tsv}"
RESULTS_ROOT="${RESULTS_ROOT:-${BASE}/wgbs_bismark/results}"
MASTER_TABLE="${MASTER_TABLE:-${BASE}/wgbs_bismark/state_methylation_by_sample_tissue/all_samples_sample_state_weighted_methylation.tsv}"
REPORT="${REPORT:-${BASE}/wgbs_bismark/state_methylation_by_sample_tissue/methylation_validation.tsv}"

[[ -s "${SAMPLESHEET}" ]] || {
  echo "[ERROR] Sample sheet was not found: ${SAMPLESHEET}" >&2
  exit 1
}

mkdir -p "$(dirname "${REPORT}")"
printf 'check\tstatus\tdetail\n' > "${REPORT}"

sample_failures=0
while IFS=$'\t' read -r sample_id library_type _ _; do
  [[ "${sample_id}" == "sample_id" ]] && continue
  sample_dir="${RESULTS_ROOT}/${sample_id}"
  cov_count="$(
    find "${sample_dir}" -maxdepth 3 -type f -name '*.bismark.cov.gz' \
      2>/dev/null | wc -l || true
  )"
  complete_marker="${sample_dir}/PIPELINE_COMPLETE"

  if [[ "${cov_count}" -ge 1 ]]; then
    printf 'sample_output\tPASS\t%s (%s; cov=%s; marker=%s)\n' \
      "${sample_id}" "${library_type}" "${cov_count}" \
      "$([[ -e "${complete_marker}" ]] && echo yes || echo no)" \
      >> "${REPORT}"
  else
    printf 'sample_output\tFAIL\t%s (%s; cov=%s; marker=%s)\n' \
      "${sample_id}" "${library_type}" "${cov_count}" \
      "$([[ -e "${complete_marker}" ]] && echo yes || echo no)" \
      >> "${REPORT}"
    sample_failures=$((sample_failures + 1))
  fi
done < "${SAMPLESHEET}"

if [[ ! -s "${MASTER_TABLE}" ]]; then
  printf 'master_table\tFAIL\tMissing: %s\n' "${MASTER_TABLE}" >> "${REPORT}"
  echo "[ERROR] Master methylation table was not found: ${MASTER_TABLE}" >&2
  exit 1
fi

if awk -F'\t' '
  NR == 1 {
    expected = "sample\tsample_tissue\tstate_tissue\tstate\tmethylated_count\tunmethylated_count\ttotal_count\tcpg_sites\tweighted_meth_pct"
    if ($0 != expected) {
      exit 1
    }
    next
  }
  {
    state = $4
    if (state !~ /^E([1-9]|10|11)$/) {
      exit 1
    }
    key = $1 SUBSEP state
    if (++seen[key] != 1) {
      exit 1
    }
    samples[$1] = 1
    states[state] = 1
    if ($9 != "NA" && ($9 + 0 < 0 || $9 + 0 > 100)) {
      exit 1
    }
  }
  END {
    for (i = 1; i <= 11; i++) {
      state = "E" i
      if (!(state in states)) {
        exit 1
      }
    }
    for (sample in samples) {
      for (i = 1; i <= 11; i++) {
        state = "E" i
        if (!(sample SUBSEP state in seen)) {
          exit 1
        }
      }
    }
  }
' "${MASTER_TABLE}"; then
  printf 'master_table\tPASS\tE1-E11 completeness, uniqueness, and range checks passed.\n' \
    >> "${REPORT}"
else
  printf 'master_table\tFAIL\tE1-E11 completeness, uniqueness, or range check failed.\n' \
    >> "${REPORT}"
  echo "[ERROR] Master table validation failed." >&2
  exit 1
fi

if (( sample_failures > 0 )); then
  echo "[ERROR] ${sample_failures} sample-level outputs failed validation." >&2
  exit 1
fi

echo "[DONE] Validation report: ${REPORT}"
