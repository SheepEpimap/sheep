#!/usr/bin/env bash

set -euo pipefail

# Figure 2j-2k, step 9: export the two chr6 contact-map regions shown in Figure 2.

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 <input.hic>" >&2
  exit 1
fi

INPUT_HIC="$1"
JUICER_TOOLS="${JUICER_TOOLS:-juicer_tools}"
CHR6="${CHR6:-chr6}"
NORMALIZATION="${NORMALIZATION:-KR}"
OUTDIR="${OUTDIR:-figure2_hic_regions}"

FIGURE_2J_REGION="${FIGURE_2J_REGION:-${CHR6}:65000000:75000000}"
FIGURE_2K_REGION="${FIGURE_2K_REGION:-${CHR6}:84445000:84505000}"
FIGURE_2J_RESOLUTION="${FIGURE_2J_RESOLUTION:-50000}"
FIGURE_2K_RESOLUTION="${FIGURE_2K_RESOLUTION:-5000}"

[[ -s "${INPUT_HIC}" ]] || {
  echo "[ERROR] Input .hic file was not found: ${INPUT_HIC}" >&2
  exit 1
}

run_juicer_tools() {
  if [[ "${JUICER_TOOLS}" == *.jar ]]; then
    java -Xmx16g -jar "${JUICER_TOOLS}" "$@"
  else
    "${JUICER_TOOLS}" "$@"
  fi
}

mkdir -p "${OUTDIR}"

run_juicer_tools dump observed "${NORMALIZATION}" \
  "${INPUT_HIC}" \
  "${FIGURE_2J_REGION}" \
  "${FIGURE_2J_REGION}" \
  BP "${FIGURE_2J_RESOLUTION}" \
  "${OUTDIR}/Figure_02j_chr6_65Mb_75Mb.${FIGURE_2J_RESOLUTION}bp.tsv"

run_juicer_tools dump observed "${NORMALIZATION}" \
  "${INPUT_HIC}" \
  "${FIGURE_2K_REGION}" \
  "${FIGURE_2K_REGION}" \
  BP "${FIGURE_2K_RESOLUTION}" \
  "${OUTDIR}/Figure_02k_TMPRSS11A.${FIGURE_2K_RESOLUTION}bp.tsv"

echo "[DONE] Figure 2j and 2k contact matrices: ${OUTDIR}"
