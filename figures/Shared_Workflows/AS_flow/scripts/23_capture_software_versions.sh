#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

OUT_DIR="${LOG_DIR}/reproducibility"
mkdir -p "${OUT_DIR}"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
REPORT="${OUT_DIR}/software_versions_${STAMP}.txt"

capture_command() {
    local label="$1"
    shift
    printf '\n[%s]\n' "${label}" >> "${REPORT}"
    if command -v "$1" >/dev/null 2>&1; then
        "$@" >> "${REPORT}" 2>&1 || true
    else
        printf 'NOT FOUND: %s\n' "$1" >> "${REPORT}"
    fi
}

{
    printf 'Pipeline software environment report\n'
    printf 'Generated (UTC): %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    printf 'Host: %s\n' "$(hostname 2>/dev/null || true)"
    printf 'Working directory: %s\n' "$(pwd)"
    printf 'Project root: %s\n' "${AS_PROJECT_ROOT}"
    printf 'Kernel: %s\n' "$(uname -a 2>/dev/null || true)"
    printf 'Bash: %s\n' "${BASH_VERSION}"
    if command -v git >/dev/null 2>&1 && git -C "${AS_PROJECT_ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        printf 'Git commit: %s\n' "$(git -C "${AS_PROJECT_ROOT}" rev-parse HEAD)"
        printf 'Git status:\n'
        git -C "${AS_PROJECT_ROOT}" status --short || true
    else
        printf 'Git commit: unavailable\n'
    fi
} > "${REPORT}"

capture_command "Conda" conda --version
capture_command "BWA" bwa 2>&1
capture_command "SAMtools" samtools --version
capture_command "BCFtools" bcftools --version
capture_command "BEDTools" bedtools --version
capture_command "MACS2" macs2 --version
capture_command "Trim Galore" trim_galore --version
capture_command "STAR" STAR --version
capture_command "Salmon" salmon --version
capture_command "StringTie" stringtie --version
capture_command "GenMap" genmap --version
capture_command "GATK" gatk --version
capture_command "R" R --version
capture_command "Python" python --version
capture_command "Java" java -version

if [[ -r "${CONDA_ROOT}/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "${CONDA_ROOT}/etc/profile.d/conda.sh"
    for env_name in "${ATAC_ENV}" "${RNASEQ_ENV}" "${WASP_ENV}" "${ASE_ENV}"; do
        printf '\n[Conda environment: %s]\n' "${env_name}" >> "${REPORT}"
        conda env export -n "${env_name}" >> "${REPORT}" 2>&1 || \
            printf 'Unable to export environment %s\n' "${env_name}" >> "${REPORT}"
    done
else
    printf '\nConda initialization script not found at %s\n' "${CONDA_ROOT}/etc/profile.d/conda.sh" >> "${REPORT}"
fi

printf 'Software environment report: %s\n' "${REPORT}"
