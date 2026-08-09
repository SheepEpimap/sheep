#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
for y in atac_tools.yml rnaseq_tools.yml wasp_env.yml ase_env.yml; do conda env create -f "${HERE}/${y}" || conda env update -f "${HERE}/${y}" --prune; done
