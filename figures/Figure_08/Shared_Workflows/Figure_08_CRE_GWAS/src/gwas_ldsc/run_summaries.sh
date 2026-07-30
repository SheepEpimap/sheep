#!/usr/bin/env bash
set -euo pipefail

: "${GWAS_BASE:?Set GWAS_BASE before running summaries}"
PY_SUMMARY="${PY_SUMMARY:-python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for script in \
    summarize_plot_A1_B1_152.py \
    summarize_plot_A2_B2_152.py \
    trait_qc_152.py \
    data_driven_discovery.py \
    data_driven_diagnostic.py \
    stat_plot_A2_B2_all_tissues_enrichment.py; do
    "$PY_SUMMARY" "$SCRIPT_DIR/$script"
done
