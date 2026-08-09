#!/usr/bin/env bash
# Shared shell helpers. Source only; do not execute directly.

as_die() {
  echo "ERROR: $*" >&2
  exit 1
}

as_warn() {
  echo "WARNING: $*" >&2
}

as_info() {
  echo "[$(date '+%F %T')] $*" >&2
}

as_require_file() {
  local path="$1"
  local label="${2:-file}"
  [[ -f "${path}" ]] || as_die "${label} not found: ${path}"
}

as_require_dir() {
  local path="$1"
  local label="${2:-directory}"
  [[ -d "${path}" ]] || as_die "${label} not found: ${path}"
}

as_require_cmds() {
  local cmd
  for cmd in "$@"; do
    command -v "${cmd}" >/dev/null 2>&1 || as_die "Required command not found in PATH: ${cmd}"
  done
}

as_activate_conda() {
  local env_name="$1"
  local conda_sh="${CONDA_ROOT}/etc/profile.d/conda.sh"
  [[ -f "${conda_sh}" ]] || as_die "Conda initialization script not found: ${conda_sh}"
  # shellcheck disable=SC1090
  source "${conda_sh}"
  conda activate "${env_name}" || as_die "Failed to activate conda environment: ${env_name}"
}

as_sample_is_rnaseq() {
  [[ "$1" == RNASeq* ]]
}

as_sample_is_atac() {
  [[ "$1" == ATAC* ]]
}

as_sample_is_peak_assay() {
  local sample="$1"
  ! as_sample_is_rnaseq "${sample}" && [[ "${sample}" != Input* ]]
}
