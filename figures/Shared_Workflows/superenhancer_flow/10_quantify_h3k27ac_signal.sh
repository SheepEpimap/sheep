#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd multiBigwigSummary
mapfile -t BIGWIGS < <(find "${H3K27AC_BW_DIR}" -maxdepth 1 -type f -name 'H3K27ac*ZScores.bw' | sort)
((${#BIGWIGS[@]} > 0)) || die "No H3K27ac Z-score bigWig files found in ${H3K27AC_BW_DIR}"

LABELS=()
if [[ -s "${H3K27AC_LABELS_FILE}" ]]; then
  mapfile -t LABELS < <(grep 'H3K27ac' "${H3K27AC_LABELS_FILE}" | awk 'NF{print $0}')
fi
if ((${#LABELS[@]} != ${#BIGWIGS[@]})); then
  log "Label count does not match bigWig count; using file basenames as labels."
  LABELS=()
  for bw in "${BIGWIGS[@]}"; do
    LABELS+=("$(basename "${bw}" .bw)")
  done
fi

run_summary() {
  local catalog="$1" stem="$2"
  require_file "${catalog}"
  multiBigwigSummary BED-file \
    -p "${THREADS}" \
    --BED "${catalog}" \
    -b "${BIGWIGS[@]}" \
    --labels "${LABELS[@]}" \
    -out "${RESULT_DIR}/${stem}.npz" \
    --outRawCounts "${RESULT_DIR}/${stem}.txt"

  sed "s/'//g" "${RESULT_DIR}/${stem}.txt" > "${RESULT_DIR}/${stem}_last.txt"
  awk '!/nan/' "${RESULT_DIR}/${stem}_last.txt" > "${RESULT_DIR}/${stem}_last_clean.txt"
}

run_summary "${RESULT_DIR}/all_super_enhancer_combine_10000.csv" "all_super_enhancer_count_10000"
run_summary "${RESULT_DIR}/all_super_enhancer_combine_1.csv" "all_super_enhancer_count_100"

log "H3K27ac signal matrices completed. No plotting is performed by this script."
