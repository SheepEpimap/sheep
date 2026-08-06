#!/usr/bin/env bash

set -euo pipefail

# Figure 2j-2k, step 10: validate sample and merged Hi-C maps with Juicer dump.

JUICER_ROOT="${JUICER_ROOT:-/vol2/zhangshiwen/juicer}"
WORK_ROOT="${WORK_ROOT:-${JUICER_ROOT}/work}"
SAMPLE_FILE="${SAMPLE_FILE:-${JUICER_ROOT}/tissue.txt}"
GROUP_TSV="${GROUP_TSV:-07_mega_group_samples.tsv}"
GROUP_ROOT="${GROUP_ROOT:-${WORK_ROOT}/mega_merge}"
JUICER_TOOLS="${JUICER_TOOLS:-juicer_tools}"
TEST_CHROMOSOME="${TEST_CHROMOSOME:-chr6}"
TEST_RESOLUTION="${TEST_RESOLUTION:-100000}"
REPORT="${REPORT:-${WORK_ROOT}/juicer_qc/hic_output_validation.tsv}"

[[ -s "${SAMPLE_FILE}" ]] || {
  echo "[ERROR] Sample file was not found: ${SAMPLE_FILE}" >&2
  exit 1
}
[[ -s "${GROUP_TSV}" ]] || {
  echo "[ERROR] Group table was not found: ${GROUP_TSV}" >&2
  exit 1
}

run_juicer_tools() {
  if [[ "${JUICER_TOOLS}" == *.jar ]]; then
    java -Xmx16g -jar "${JUICER_TOOLS}" "$@"
  else
    "${JUICER_TOOLS}" "$@"
  fi
}

mkdir -p "$(dirname "${REPORT}")"
printf 'level\tname\thic_file\tstatus\tdump_rows\n' > "${REPORT}"
failures=0

validate_hic() {
  local level="$1"
  local name="$2"
  local hic_file="$3"
  local temporary_dump row_count

  if [[ ! -s "${hic_file}" ]]; then
    printf '%s\t%s\t%s\tFAIL_MISSING\t0\n' \
      "${level}" "${name}" "${hic_file}" >> "${REPORT}"
    failures=$((failures + 1))
    return
  fi

  temporary_dump="$(mktemp)"
  if run_juicer_tools dump observed NONE \
    "${hic_file}" \
    "${TEST_CHROMOSOME}" \
    "${TEST_CHROMOSOME}" \
    BP "${TEST_RESOLUTION}" \
    "${temporary_dump}" >/dev/null 2>&1; then
    row_count="$(awk 'NF { n++ } END { print n + 0 }' "${temporary_dump}")"
  else
    row_count=0
  fi
  rm -f -- "${temporary_dump}"

  if (( row_count > 0 )); then
    printf '%s\t%s\t%s\tPASS\t%s\n' \
      "${level}" "${name}" "${hic_file}" "${row_count}" >> "${REPORT}"
  else
    printf '%s\t%s\t%s\tFAIL_DUMP\t0\n' \
      "${level}" "${name}" "${hic_file}" >> "${REPORT}"
    failures=$((failures + 1))
  fi
}

while read -r sample_id _; do
  [[ -z "${sample_id}" ]] && continue
  validate_hic \
    "sample" \
    "${sample_id}" \
    "${WORK_ROOT}/${sample_id}/aligned/inter_30.hic"
done < "${SAMPLE_FILE}"

mapfile -t groups < <(
  cut -f1 "${GROUP_TSV}" | awk 'NF' | LC_ALL=C sort -u
)
for group in "${groups[@]}"; do
  validate_hic \
    "group" \
    "${group}" \
    "${GROUP_ROOT}/${group}/${group}.MAPQ30.hic"
done

echo "[DONE] Hi-C validation report: ${REPORT}"
if (( failures > 0 )); then
  echo "[ERROR] Hi-C output failures: ${failures}" >&2
  exit 1
fi
