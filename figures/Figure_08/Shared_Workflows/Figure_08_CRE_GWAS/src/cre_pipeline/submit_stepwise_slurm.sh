#!/usr/bin/env bash
set -euo pipefail

: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR before submission}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${WORKDIR}/logs"
mkdir -p "$LOG_DIR"

submit() {
    local step="$1" partition="$2" cpus="$3" mem="$4" dependency="$5" script="$6"
    local args=(
        --parsable
        --job-name="E1E9_${step}"
        --partition="$partition"
        --nodes=1
        --ntasks=1
        --cpus-per-task="$cpus"
        --mem="$mem"
        --chdir="$WORKDIR"
        --output="$LOG_DIR/${step}_%j.out"
        --error="$LOG_DIR/${step}_%j.err"
        --export=ALL
    )
    [[ -n "$dependency" ]] && args+=(--dependency="afterok:${dependency}")
    local jid
    jid=$(sbatch "${args[@]}" --wrap="bash '$SCRIPT_DIR/$script'")
    printf '%s\n' "${jid%%;*}"
}

j1=$(submit 01_reclassify low 32 64G "" 01_reclassify_e1_e9.sh)
j2=$(submit 02_hg38_to_hg19 low 32 64G "$j1" 02_hg38_to_hg19.sh)
j3=$(submit 03_collect_e5 smp 1 8G "$j2" 03_collect_e5_hg19.sh)
j4=$(submit 04_reuse_human_tsr smp 1 8G "$j3" 04_reuse_human_tsr.sh)
j5=$(submit 05_label_tsr low 32 64G "$j4" 05_label_cre_with_tsr.sh)
j6=$(submit 06_reuse_projection smp 1 8G "$j5" 06_reuse_sheep_projection.sh)
j7=$(submit 07_split_ldsc low 8 32G "$j6" 07_split_ldsc_beds.sh)
j8=$(submit 08_collect_gwas_beds low 8 16G "$j7" 08_collect_gwas_beds.sh)

printf 'step\tjob_id\n'
printf '01_reclassify\t%s\n02_hg38_to_hg19\t%s\n03_collect_e5\t%s\n04_reuse_human_tsr\t%s\n05_label_tsr\t%s\n06_reuse_projection\t%s\n07_split_ldsc\t%s\n08_collect_gwas_beds\t%s\n' \
    "$j1" "$j2" "$j3" "$j4" "$j5" "$j6" "$j7" "$j8"
