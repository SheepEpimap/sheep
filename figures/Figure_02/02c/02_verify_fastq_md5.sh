#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 2: verify downloaded FASTQ files against ENA MD5 values.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
FASTQDIR="${FASTQDIR:-${BASE}/fastq}"
MD5_TSV="${MD5_TSV:-${BASE}/md5sum.txt}"
THREADS="${THREADS:-8}"

OUT_DETAIL="${FASTQDIR}/md5_check_detail.tsv"
OUT_FAILED="${FASTQDIR}/md5_failed_files.txt"

[[ -d "${FASTQDIR}" ]] || {
  echo "[ERROR] FASTQ directory was not found: ${FASTQDIR}" >&2
  exit 1
}
[[ -s "${MD5_TSV}" ]] || {
  echo "[ERROR] MD5 metadata table was not found: ${MD5_TSV}" >&2
  exit 1
}

declare -A expected
while IFS=$'\t' read -r run_accession _ md5_pair _; do
  run_accession="${run_accession%$'\r'}"
  md5_pair="${md5_pair%$'\r'}"
  [[ -z "${run_accession}" || "${run_accession}" == "run_accession" ]] && continue

  IFS=';' read -r md5_1 md5_2 <<< "${md5_pair}"
  expected["${run_accession}_1"]="$(
    printf '%s' "${md5_1:-}" | tr -d ' ' | tr '[:upper:]' '[:lower:]'
  )"
  expected["${run_accession}_2"]="$(
    printf '%s' "${md5_2:-}" | tr -d ' ' | tr '[:upper:]' '[:lower:]'
  )"
done < "${MD5_TSV}"

tmp_md5="$(mktemp)"
trap 'rm -f -- "${tmp_md5}"' EXIT

find "${FASTQDIR}" -maxdepth 1 -type f -name 'SRR*_*.fastq.gz' -print0 \
  | xargs -0 -r -n 1 -P "${THREADS}" md5sum \
  > "${tmp_md5}"

printf 'run\tmate\tfilename\texpected_md5\tactual_md5\tstatus\n' > "${OUT_DETAIL}"
: > "${OUT_FAILED}"

declare -A observed
while read -r actual path; do
  filename="$(basename "${path}")"
  if [[ "${filename}" =~ ^(SRR[0-9]+)_([12])\.fastq\.gz$ ]]; then
    run_accession="${BASH_REMATCH[1]}"
    mate="${BASH_REMATCH[2]}"
    key="${run_accession}_${mate}"
    observed["${key}"]=1
    actual="$(printf '%s' "${actual}" | tr '[:upper:]' '[:lower:]')"
    expected_md5="${expected[${key}]:-}"

    if [[ -z "${expected_md5}" ]]; then
      status="NO_EXPECTED_MD5"
    elif [[ "${actual}" != "${expected_md5}" ]]; then
      status="MISMATCH"
    else
      status="OK"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${run_accession}" "${mate}" "${filename}" \
      "${expected_md5}" "${actual}" "${status}" >> "${OUT_DETAIL}"

    [[ "${status}" == "OK" ]] || printf '%s\n' "${filename}" >> "${OUT_FAILED}"
  fi
done < "${tmp_md5}"

for key in "${!expected[@]}"; do
  if [[ -z "${observed[${key}]:-}" ]]; then
    run_accession="${key%_*}"
    mate="${key##*_}"
    filename="${run_accession}_${mate}.fastq.gz"
    printf '%s\t%s\t%s\t%s\t\tMISSING_FILE\n' \
      "${run_accession}" "${mate}" "${filename}" "${expected[${key}]}" \
      >> "${OUT_DETAIL}"
    printf '%s\n' "${filename}" >> "${OUT_FAILED}"
  fi
done

LC_ALL=C sort -u -o "${OUT_FAILED}" "${OUT_FAILED}"
failed_count="$(wc -l < "${OUT_FAILED}" | tr -d ' ')"

echo "[DONE] Detailed MD5 report: ${OUT_DETAIL}"
echo "[INFO] Failed or missing FASTQ files: ${failed_count}"
[[ "${failed_count}" -eq 0 ]]
