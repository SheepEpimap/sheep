#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"; source "${SCRIPT_DIR}/00_config.sh"
summary="${BINOMIAL_DIR}/summary_line_counts.tsv"; printf 'file\tdata_rows\n' > "${summary}"
find "${BINOMIAL_DIR}" -maxdepth 1 -type f -name '*_binomial_results.txt' -print0 | sort -z | while IFS= read -r -d '' file; do n=$(wc -l < "${file}"); ((n>0))&&n=$((n-1)); printf '%s\t%s\n' "$(basename "${file}")" "${n}"; done >> "${summary}"
echo "${summary}"
