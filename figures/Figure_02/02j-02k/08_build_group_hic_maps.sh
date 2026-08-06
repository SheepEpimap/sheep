#!/usr/bin/env bash

set -euo pipefail

# Figure 2j-2k, step 8: merge selected samples and build tissue-group .hic maps.

JUICER_ROOT="${JUICER_ROOT:-/vol2/zhangshiwen/juicer}"
WORK_ROOT="${WORK_ROOT:-${JUICER_ROOT}/work}"
GROUP_TSV="${GROUP_TSV:-07_mega_group_samples.tsv}"
OUT_ROOT="${OUT_ROOT:-${WORK_ROOT}/mega_merge}"
CHROM_SIZES="${CHROM_SIZES:-${JUICER_ROOT}/restriction_sites/sheep.chrom.sizes}"
JUICER_TOOLS="${JUICER_TOOLS:-juicer_tools}"
MAPQ="${MAPQ:-30}"
RESOLUTIONS="${RESOLUTIONS:-2500000,1000000,500000,250000,100000,50000,25000,10000,5000}"
REBUILD="${REBUILD:-0}"

[[ -s "${GROUP_TSV}" ]] || {
  echo "[ERROR] Group table was not found: ${GROUP_TSV}" >&2
  exit 1
}
[[ -s "${CHROM_SIZES}" ]] || {
  echo "[ERROR] Chromosome-size file was not found: ${CHROM_SIZES}" >&2
  exit 1
}

run_juicer_tools() {
  if [[ "${JUICER_TOOLS}" == *.jar ]]; then
    java -Xmx48g -jar "${JUICER_TOOLS}" "$@"
  else
    "${JUICER_TOOLS}" "$@"
  fi
}

if [[ "${JUICER_TOOLS}" == *.jar ]]; then
  [[ -s "${JUICER_TOOLS}" ]] || {
    echo "[ERROR] Juicer Tools JAR was not found: ${JUICER_TOOLS}" >&2
    exit 1
  }
  command -v java >/dev/null 2>&1 || {
    echo "[ERROR] java was not found." >&2
    exit 1
  }
else
  command -v "${JUICER_TOOLS}" >/dev/null 2>&1 || {
    echo "[ERROR] Juicer Tools command was not found: ${JUICER_TOOLS}" >&2
    exit 1
  }
fi

mkdir -p "${OUT_ROOT}"
mapfile -t groups < <(
  cut -f1 "${GROUP_TSV}" | awk 'NF' | LC_ALL=C sort -u
)

for group in "${groups[@]}"; do
  group_dir="${OUT_ROOT}/${group}"
  merged="${group_dir}/merged_nodups.txt"
  output_hic="${group_dir}/${group}.MAPQ${MAPQ}.hic"

  mkdir -p "${group_dir}"
  if [[ "${REBUILD}" != "1" && -s "${output_hic}" ]]; then
    echo "[SKIP] Existing group map: ${output_hic}"
    continue
  fi

  : > "${merged}"
  sample_count=0
  while IFS=$'\t' read -r table_group sample_id; do
    [[ "${table_group}" == "${group}" ]] || continue
    source_pairs="${WORK_ROOT}/${sample_id}/aligned/merged_nodups.txt"
    if [[ ! -s "${source_pairs}" ]]; then
      echo "[ERROR] Valid-pair file was not found: ${source_pairs}" >&2
      exit 1
    fi
    awk 'NF' "${source_pairs}" >> "${merged}"
    sample_count=$((sample_count + 1))
  done < "${GROUP_TSV}"

  (( sample_count > 0 )) || {
    echo "[ERROR] No samples were assigned to group ${group}." >&2
    exit 1
  }
  [[ -s "${merged}" ]] || {
    echo "[ERROR] Merged valid-pair file is empty: ${merged}" >&2
    exit 1
  }

  run_juicer_tools pre \
    -n \
    -q "${MAPQ}" \
    -r "${RESOLUTIONS}" \
    "${merged}" \
    "${output_hic}" \
    "${CHROM_SIZES}"

  run_juicer_tools addNorm "${output_hic}"
  [[ -s "${output_hic}" ]] || {
    echo "[ERROR] Group Hi-C map was not generated: ${output_hic}" >&2
    exit 1
  }

  echo "[DONE] ${group}: ${sample_count} samples -> ${output_hic}"
done
