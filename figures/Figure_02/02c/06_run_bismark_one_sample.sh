#!/usr/bin/env bash

set -euo pipefail

# Figure 2c, step 6: process one paired-end RRBS or WGBS sample with Bismark.

if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 <sample_id> <RRBS|WGBS> <R1.fastq.gz> <R2.fastq.gz>" >&2
  exit 1
fi

SAMPLE_ID="$1"
LIBRARY_TYPE="$2"
R1="$3"
R2="$4"

BASE="${BASE:-/vol2/zhangshiwen/rrbs}"
SHEEP_REFERENCE="${SHEEP_REFERENCE:-${BASE}/wgbs_bismark/ref/sheep}"
LAMBDA_REFERENCE="${LAMBDA_REFERENCE:-${BASE}/wgbs_bismark/ref/lambda}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-16}}"
NON_DIRECTIONAL="${NON_DIRECTIONAL:-0}"
RUN_LAMBDA="${RUN_LAMBDA:-0}"
BISMARK_PARALLEL="${BISMARK_PARALLEL:-1}"
BUFFER_SIZE="${BUFFER_SIZE:-10G}"

[[ "${LIBRARY_TYPE}" =~ ^(RRBS|WGBS)$ ]] || {
  echo "[ERROR] Library type must be RRBS or WGBS: ${LIBRARY_TYPE}" >&2
  exit 1
}
[[ -s "${R1}" && -s "${R2}" ]] || {
  echo "[ERROR] Paired FASTQ input is missing or empty." >&2
  exit 1
}

OUT="${BASE}/wgbs_bismark/results/${SAMPLE_ID}"
TMP="${BASE}/wgbs_bismark/tmp/${SAMPLE_ID}"
TRIMDIR="${OUT}/01_trim"
ALIGN_SHEEP="${OUT}/02_align_sheep"
DEDUP_SHEEP="${OUT}/03_dedup_sheep"
METH_SHEEP="${OUT}/04_meth_sheep"
ALIGN_LAMBDA="${OUT}/05_align_lambda"
METH_LAMBDA="${OUT}/06_meth_lambda"

mkdir -p \
  "${TRIMDIR}" "${ALIGN_SHEEP}" "${DEDUP_SHEEP}" "${METH_SHEEP}" \
  "${ALIGN_LAMBDA}" "${METH_LAMBDA}" "${TMP}"

trim_cores="${THREADS}"
(( trim_cores > 8 )) && trim_cores=8
extractor_cores="$((THREADS / 3))"
(( extractor_cores < 1 )) && extractor_cores=1

trim_extra=()
[[ "${LIBRARY_TYPE}" == "RRBS" ]] && trim_extra+=(--rrbs)

trim_galore \
  --paired \
  --gzip \
  --quality 20 \
  --length 20 \
  --cores "${trim_cores}" \
  --fastqc \
  --basename "${SAMPLE_ID}" \
  "${trim_extra[@]}" \
  -o "${TRIMDIR}" \
  "${R1}" "${R2}"

TRIMMED_R1="${TRIMDIR}/${SAMPLE_ID}_val_1.fq.gz"
TRIMMED_R2="${TRIMDIR}/${SAMPLE_ID}_val_2.fq.gz"
[[ -s "${TRIMMED_R1}" && -s "${TRIMMED_R2}" ]] || {
  echo "[ERROR] Trim Galore did not produce both paired outputs." >&2
  exit 1
}

bismark_extra=()
[[ "${NON_DIRECTIONAL}" == "1" ]] && bismark_extra+=(--non_directional)

bismark \
  --genome "${SHEEP_REFERENCE}" \
  --bowtie2 \
  --parallel "${BISMARK_PARALLEL}" \
  --temp_dir "${TMP}/sheep" \
  "${bismark_extra[@]}" \
  -1 "${TRIMMED_R1}" \
  -2 "${TRIMMED_R2}" \
  -o "${ALIGN_SHEEP}"

mapfile -t sheep_bams < <(
  find "${ALIGN_SHEEP}" -maxdepth 1 -type f -name '*bismark*pe.bam' \
    -print | LC_ALL=C sort
)
(( ${#sheep_bams[@]} > 0 )) || {
  echo "[ERROR] Sheep-aligned Bismark BAM was not produced." >&2
  exit 1
}
sheep_bam="${sheep_bams[0]}"

extraction_bam="${sheep_bam}"
if [[ "${LIBRARY_TYPE}" == "WGBS" ]]; then
  (
    cd "${DEDUP_SHEEP}"
    deduplicate_bismark --bam --paired "${sheep_bam}"
  )
  mapfile -t deduplicated_bams < <(
    find "${DEDUP_SHEEP}" -maxdepth 1 -type f -name '*deduplicated.bam' \
      -print | LC_ALL=C sort
  )
  (( ${#deduplicated_bams[@]} > 0 )) || {
    echo "[ERROR] Deduplicated WGBS BAM was not produced." >&2
    exit 1
  }
  extraction_bam="${deduplicated_bams[0]}"
fi

bismark_methylation_extractor \
  -p \
  --multicore "${extractor_cores}" \
  --gzip \
  --bedGraph \
  --counts \
  --comprehensive \
  --genome_folder "${SHEEP_REFERENCE}" \
  --no_overlap \
  --buffer_size "${BUFFER_SIZE}" \
  -o "${METH_SHEEP}" \
  "${extraction_bam}"

if [[ "${RUN_LAMBDA}" == "1" ]]; then
  bismark \
    --genome "${LAMBDA_REFERENCE}" \
    --bowtie2 \
    --parallel "${BISMARK_PARALLEL}" \
    --temp_dir "${TMP}/lambda" \
    "${bismark_extra[@]}" \
    -1 "${TRIMMED_R1}" \
    -2 "${TRIMMED_R2}" \
    -o "${ALIGN_LAMBDA}"

  mapfile -t lambda_bams < <(
    find "${ALIGN_LAMBDA}" -maxdepth 1 -type f -name '*bismark*pe.bam' \
      -print | LC_ALL=C sort
  )
  if (( ${#lambda_bams[@]} > 0 )); then
    bismark_methylation_extractor \
      -p \
      --multicore "${extractor_cores}" \
      --gzip \
      --bedGraph \
      --counts \
      --CX_context \
      --comprehensive \
      --genome_folder "${LAMBDA_REFERENCE}" \
      --no_overlap \
      --buffer_size "${BUFFER_SIZE}" \
      -o "${METH_LAMBDA}" \
      "${lambda_bams[0]}"
  else
    echo "[WARN] Lambda-aligned BAM was not produced; conversion check skipped." >&2
  fi
fi

touch "${OUT}/PIPELINE_COMPLETE"
echo "[DONE] ${SAMPLE_ID}"
