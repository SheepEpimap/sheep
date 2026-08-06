#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 5: create a paired-end sample sheet with RRBS/WGBS labels.

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
FASTQDIR="${FASTQDIR:-${BASE}/fastq}"
SAMPLESHEET="${SAMPLESHEET:-${BASE}/wgbs_bismark/samples.tsv}"

RRBS_IDS=(
  SRR15616174 SRR15616175 SRR15616176 SRR15616177
  SRR15616178 SRR15616179 SRR15616180 SRR15616181
  SRR15616182 SRR15616183 SRR15616184 SRR15616185
)

is_rrbs() {
  local sample_id="$1"
  local rrbs_id
  for rrbs_id in "${RRBS_IDS[@]}"; do
    [[ "${sample_id}" == "${rrbs_id}" ]] && return 0
  done
  return 1
}

mkdir -p "$(dirname "${SAMPLESHEET}")"
printf 'sample_id\tlibtype\tr1\tr2\n' > "${SAMPLESHEET}"

shopt -s nullglob
r1_files=("${FASTQDIR}"/*_1.fastq.gz)
if (( ${#r1_files[@]} == 0 )); then
  echo "[ERROR] No paired-end FASTQ files were found in ${FASTQDIR}." >&2
  exit 1
fi

for r1 in "${r1_files[@]}"; do
  sample_id="$(basename "${r1}" _1.fastq.gz)"
  r2="${FASTQDIR}/${sample_id}_2.fastq.gz"
  if [[ ! -s "${r2}" ]]; then
    echo "[WARN] Missing R2 file for ${sample_id}; sample omitted." >&2
    continue
  fi

  library_type="WGBS"
  is_rrbs "${sample_id}" && library_type="RRBS"
  printf '%s\t%s\t%s\t%s\n' \
    "${sample_id}" "${library_type}" "${r1}" "${r2}" \
    >> "${SAMPLESHEET}"
done

sample_count="$(( $(wc -l < "${SAMPLESHEET}") - 1 ))"
echo "[DONE] Sample sheet: ${SAMPLESHEET}"
echo "[INFO] Samples: ${sample_count}"
