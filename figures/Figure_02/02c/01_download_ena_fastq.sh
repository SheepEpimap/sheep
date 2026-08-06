#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 1: download RRBS/WGBS paired-end FASTQ files from ENA.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
METADATA_TSV="${METADATA_TSV:-${BASE}/select_RRBSresults_read_run_tsv.txt}"
ASPERA_COLUMN="${ASPERA_COLUMN:-2}"
OUTDIR="${OUTDIR:-${BASE}/fastq}"
PARALLEL_DOWNLOADS="${PARALLEL_DOWNLOADS:-6}"
RATE_LIMIT="${RATE_LIMIT:-400m}"
ASPERA_PORT="${ASPERA_PORT:-33001}"
RESUME_POLICY="${RESUME_POLICY:-2}"
MAX_RETRIES="${MAX_RETRIES:-3}"

[[ -s "${METADATA_TSV}" ]] || {
  echo "[ERROR] Metadata table was not found or is empty: ${METADATA_TSV}" >&2
  exit 1
}

command -v ascp >/dev/null 2>&1 || {
  echo "[ERROR] ascp was not found in PATH." >&2
  exit 1
}

key_candidates=(
  "${ASPERA_KEY:-}"
  "${CONDA_PREFIX:-}/etc/asperaweb_id_dsa.openssh"
  "${CONDA_PREFIX:-}/etc/aspera/asperaweb_id_dsa.openssh"
  "${HOME}/.aspera/connect/etc/asperaweb_id_dsa.openssh"
  "/etc/asperaweb_id_dsa.openssh"
)

aspera_key=""
for candidate in "${key_candidates[@]}"; do
  if [[ -n "${candidate}" && -f "${candidate}" ]]; then
    aspera_key="${candidate}"
    break
  fi
done

[[ -n "${aspera_key}" ]] || {
  echo "[ERROR] The ENA Aspera private key was not found." >&2
  exit 1
}

mkdir -p "${OUTDIR}" "${OUTDIR}/logs"
paths_file="${OUTDIR}/aspera_paths.txt"

cut -f"${ASPERA_COLUMN}" "${METADATA_TSV}" \
  | tr -d '\r' \
  | sed 's/[;,]/\n/g' \
  | sed 's/^fasp\.sra\.ebi\.ac\.uk://' \
  | awk 'NF' \
  | LC_ALL=C sort -u \
  > "${paths_file}"

path_count="$(wc -l < "${paths_file}" | tr -d ' ')"
if [[ "${path_count}" -eq 0 ]]; then
  echo "[ERROR] No ENA Aspera paths were parsed from column ${ASPERA_COLUMN}." >&2
  exit 1
fi

export OUTDIR RATE_LIMIT ASPERA_PORT RESUME_POLICY MAX_RETRIES aspera_key

download_one() {
  local remote_path="$1"
  local filename log_file attempt
  filename="$(basename "${remote_path}")"
  log_file="${OUTDIR}/logs/${filename}.log"

  if [[ -s "${OUTDIR}/${filename}" &&
        ! -e "${OUTDIR}/${filename}.aspx" &&
        ! -e "${OUTDIR}/${filename}.aspera-ckpt" ]]; then
    echo "[SKIP] Complete file: ${filename}" >> "${log_file}"
    return 0
  fi

  for ((attempt = 1; attempt <= MAX_RETRIES; attempt++)); do
    ascp_args=(
      -QT
      -P "${ASPERA_PORT}"
      -k "${RESUME_POLICY}"
      --overwrite=diff
      -i "${aspera_key}"
    )
    if [[ -n "${RATE_LIMIT}" ]]; then
      ascp_args+=(-l "${RATE_LIMIT}")
    fi

    if ascp "${ascp_args[@]}" \
      "era-fasp@fasp.sra.ebi.ac.uk:${remote_path}" \
      "${OUTDIR}/" >> "${log_file}" 2>&1; then
      return 0
    fi

    echo "[WARN] Attempt ${attempt} failed: ${filename}" >> "${log_file}"
    sleep "$((attempt * 2))"
  done

  echo "[ERROR] Download failed after ${MAX_RETRIES} attempts: ${filename}" \
    >> "${log_file}"
  return 1
}

export -f download_one
xargs -r -P "${PARALLEL_DOWNLOADS}" -n 1 bash -c \
  'download_one "$1"' _ < "${paths_file}"

remaining_metadata="$(
  find "${OUTDIR}" -maxdepth 1 -type f \
    \( -name '*.aspx' -o -name '*.aspera-ckpt' \) \
    | wc -l
)"

echo "[DONE] FASTQ directory: ${OUTDIR}"
echo "[INFO] Parsed files: ${path_count}"
echo "[INFO] Remaining transfer metadata files: ${remaining_metadata}"
