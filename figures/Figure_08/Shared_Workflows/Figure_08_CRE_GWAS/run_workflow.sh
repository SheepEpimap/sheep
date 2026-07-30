#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"

usage() {
    echo "Usage: bash run_workflow.sh <paths.tsv> <run_id> [--dry-run|--execute]" >&2
    echo "Example: bash run_workflow.sh config/paths.server.tsv 20260729_fig8 --dry-run" >&2
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
    usage
    exit 2
fi

PATHS_CONFIG="$1"
RUN_ID="$2"
RUN_MODE="${3:---dry-run}"

if [[ "${PATHS_CONFIG}" != /* ]]; then
    PATHS_CONFIG="${PROJECT_DIR}/${PATHS_CONFIG}"
fi

if [[ ! -f "${PATHS_CONFIG}" ]]; then
    echo "Path configuration not found: ${PATHS_CONFIG}" >&2
    echo "Copy config/paths.example.tsv to config/paths.server.tsv and edit it first." >&2
    exit 1
fi

case "${RUN_MODE}" in
    --dry-run)
        SUBMIT_ARGS=(--dry-run)
        ;;
    --execute)
        SUBMIT_ARGS=()
        ;;
    *)
        usage
        exit 2
        ;;
esac

cd "${PROJECT_DIR}"

"${PYTHON_BIN}" workflow.py \
    --paths "${PATHS_CONFIG}" \
    --stage preflight \
    --run-id "${RUN_ID}"

"${PYTHON_BIN}" workflow.py \
    --paths "${PATHS_CONFIG}" \
    --stage cre_submit_slurm \
    --run-id "${RUN_ID}" \
    "${SUBMIT_ARGS[@]}"

"${PYTHON_BIN}" workflow.py \
    --paths "${PATHS_CONFIG}" \
    --stage gwas_submit_slurm \
    --run-id "${RUN_ID}" \
    "${SUBMIT_ARGS[@]}"
