#!/usr/bin/env bash
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/00_config.sh"

missing=0

check_file() {
  local label="$1" path="$2"
  if [[ -f "${path}" ]]; then
    printf "[OK]      %-28s %s\n" "${label}" "${path}"
  else
    printf "[MISSING] %-28s %s\n" "${label}" "${path}"
    missing=$((missing + 1))
  fi
}

check_dir() {
  local label="$1" path="$2"
  if [[ -d "${path}" ]]; then
    printf "[OK]      %-28s %s\n" "${label}" "${path}"
  else
    printf "[MISSING] %-28s %s\n" "${label}" "${path}"
    missing=$((missing + 1))
  fi
}

check_cmd() {
  local cmd="$1"
  if command -v "${cmd}" >/dev/null 2>&1; then
    printf "[OK]      command %-20s %s\n" "${cmd}" "$(command -v "${cmd}")"
  else
    printf "[MISSING] command %-20s\n" "${cmd}"
    missing=$((missing + 1))
  fi
}

echo "AS pipeline portability/environment check"
echo "Project root: ${AS_PROJECT_ROOT}"
echo

check_dir  "Input BAM directory" "${INPUT_BAM_DIR}"
check_dir  "RNA-seq raw reads" "${RNASEQ_RAW_READS_DIR}"
check_file "Reference FASTA" "${REFERENCE_FASTA}"
check_file "Reference GTF" "${REFERENCE_GTF}"
check_file "Phasing VCF" "${PHASING_VCF}"
check_dir  "STAR index" "${STAR_INDEX}"
check_dir  "Salmon index" "${SALMON_INDEX}"
check_dir  "WASP installation" "${WASP_PATH}"
check_file "SPP runner" "${SPP_RUNNER}"
for generated in "${WASP_HAPS_H5}" "${WASP_SNP_TAB_H5}" "${WASP_SNP_INDEX_H5}"; do
  if [[ -f "${generated}" ]]; then
    printf "[OK]      %-28s %s\n" "generated WASP HDF5" "${generated}"
  else
    printf "[PENDING] %-28s %s (created by step 09)\n" "generated WASP HDF5" "${generated}"
  fi
done

echo
for cmd in bash awk sed grep find samtools bedtools bcftools bgzip tabix; do
  check_cmd "${cmd}"
done

echo
if [[ "${missing}" -eq 0 ]]; then
  echo "All checked external paths and base commands are available."
  exit 0
fi

echo "${missing} required external paths or base commands are missing."
echo "Edit scripts/00_config.sh or export replacement paths before running the affected stage."
exit 2
