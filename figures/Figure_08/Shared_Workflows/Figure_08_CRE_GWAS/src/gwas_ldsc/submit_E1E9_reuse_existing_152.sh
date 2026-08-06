#!/usr/bin/env bash
set -euo pipefail

: "${GWAS_BASE:?Set GWAS_BASE in the paths configuration}"
: "${LDSC_ROOT:?Set LDSC_ROOT in the paths configuration}"
: "${LDSC_PYTHON:?Set LDSC_PYTHON in the paths configuration}"
: "${LDSC_REF:?Set LDSC_REF in the paths configuration}"
PY_SUMMARY="${PY_SUMMARY:-python}"
ARRAY_CHUNK_SIZE="${ARRAY_CHUNK_SIZE:-1000}"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-80}"

if [[ "$ARRAY_CHUNK_SIZE" -lt 1 ]]; then
    echo "[ERROR] ARRAY_CHUNK_SIZE must be >= 1" >&2
    exit 1
fi

if [[ "$ARRAY_CHUNK_SIZE" -gt 1000 ]]; then
    echo "[WARN] ARRAY_CHUNK_SIZE=$ARRAY_CHUNK_SIZE exceeds the conservative 1000-task array size; using 1000" >&2
    ARRAY_CHUNK_SIZE=1000
fi

if [[ "$ARRAY_CONCURRENCY" -lt 1 ]]; then
    echo "[ERROR] ARRAY_CONCURRENCY must be >= 1" >&2
    exit 1
fi

if [[ "$ARRAY_CONCURRENCY" -gt 80 ]]; then
    echo "[WARN] ARRAY_CONCURRENCY=$ARRAY_CONCURRENCY exceeds 80; using 80" >&2
    ARRAY_CONCURRENCY=80
fi

cd "$GWAS_BASE"
mkdir -p logs outputs summary figures

line_count() {
    local manifest="$1"
    if [[ ! -s "$manifest" ]]; then
        echo "[ERROR] missing or empty manifest: $manifest" >&2
        exit 1
    fi
    wc -l < "$manifest"
}

expected_first_column() {
    case "$(basename "$1")" in
        make_annot_jobs.tsv|control_make_annot_jobs.tsv|ldscore_jobs.tsv|control_ldscore_jobs.tsv|A1_h2_manifest_152.tsv|B1_h2_manifest_152.tsv)
            printf 'annotation_id\n' ;;
        A2_B2_cts_manifest_152.tsv)
            printf 'ldcts_name\n' ;;
        *)
            echo "[ERROR] no registered schema for manifest: $1" >&2
            return 1 ;;
    esac
}

job_count() {
    local manifest="$1"
    if [[ ! -s "$manifest" ]]; then
        echo "[ERROR] missing or empty manifest: $manifest" >&2
        exit 1
    fi

    local expected
    local observed
    local n
    expected=$(expected_first_column "$manifest")
    observed=$(head -n 1 "$manifest" | cut -f1 | tr -d '\r')
    if [[ "$observed" != "$expected" ]]; then
        echo "[ERROR] manifest header mismatch: $manifest expected=$expected observed=$observed" >&2
        echo "[ERROR] normalize manifests with prepare_E1E9_reuse_existing_recalc.sh before submission" >&2
        exit 1
    fi
    n=$(($(wc -l < "$manifest") - 1))

    if [[ "$n" -le 0 ]]; then
        echo "[ERROR] manifest has no data rows: $manifest" >&2
        exit 1
    fi

    echo "$n"
}

submit_array() {
    local manifest="$1"
    local script="$2"
    local name="$3"
    shift 3
    local n
    n=$(job_count "$manifest")
    if [[ "$n" -le 0 ]]; then
        echo "[ERROR] manifest has no jobs: $manifest" >&2
        exit 1
    fi

    local offset=0
    local chunk_count
    local last
    local chunk_idx=1
    local jid
    local ids=()

    while [[ "$offset" -lt "$n" ]]; do
        chunk_count=$((n - offset))
        if [[ "$chunk_count" -gt "$ARRAY_CHUNK_SIZE" ]]; then
            chunk_count="$ARRAY_CHUNK_SIZE"
        fi
        last=$((chunk_count - 1))
        jid=$(
            sbatch --parsable \
                --job-name="${name}_${chunk_idx}" \
                --array="0-${last}%${ARRAY_CONCURRENCY}" \
                --chdir="$GWAS_BASE" \
                --output="$GWAS_BASE/logs/${name}_%A_%a.log" \
                --error="$GWAS_BASE/logs/${name}_%A_%a.err" \
                --export=ALL,GWAS_BASE="$GWAS_BASE",OFFSET="$offset",LDSC_ROOT="$LDSC_ROOT",LDSC_PYTHON="$LDSC_PYTHON",LDSC_REF="$LDSC_REF" \
                "$@" \
                "$script"
        )
        jid="${jid%%;*}"
        ids+=("$jid")
        echo "[INFO] submitted $name chunk=$chunk_idx offset=$offset array=0-${last}%${ARRAY_CONCURRENCY} job=$jid" >&2
        offset=$((offset + chunk_count))
        chunk_idx=$((chunk_idx + 1))
    done

    local IFS=:
    echo "${ids[*]}"
}

echo "[INFO] submitting E1-E9 reuse-existing 152-trait GWAS recalculation"
echo "[INFO] GWAS_BASE=$GWAS_BASE"
echo "[INFO] ARRAY_CHUNK_SIZE=$ARRAY_CHUNK_SIZE"
echo "[INFO] ARRAY_CONCURRENCY=$ARRAY_CONCURRENCY"

make_jid=$(submit_array manifests/make_annot_jobs.tsv scripts/run_make_annots_array.sh E1E9_make_annot)
ctrl_make_jid=$(submit_array manifests/control_make_annot_jobs.tsv scripts/run_control_make_annots.sh E1E9_ctrl_annot)

ld_dep="afterok:${make_jid}:${ctrl_make_jid}"
ld_jid=$(submit_array manifests/ldscore_jobs.tsv scripts/run_ldscore_array.sh E1E9_ldscore --dependency="$ld_dep")
ctrl_ld_jid=$(submit_array manifests/control_ldscore_jobs.tsv scripts/run_control_ldscore.sh E1E9_ctrl_ldscore --dependency="$ld_dep")

h2_dep="afterok:${ld_jid}:${ctrl_ld_jid}"
a1_jid=$(submit_array manifests/A1_h2_manifest_152.tsv scripts/run_A1_h2_152.sh E1E9_A1_h2_152 --dependency="$h2_dep")
b1_jid=$(submit_array manifests/B1_h2_manifest_152.tsv scripts/run_B1_h2_152.sh E1E9_B1_h2_152 --dependency="$h2_dep")
cts_jid=$(submit_array manifests/A2_B2_cts_manifest_152.tsv scripts/run_A2_B2_cts_152.sh E1E9_A2B2_cts_152 --dependency="$h2_dep")

summary_dep="afterok:${a1_jid}:${b1_jid}:${cts_jid}"
summary_jid=$(
    sbatch --parsable \
        --job-name=E1E9_summary_152 \
        --dependency="$summary_dep" \
        --output="$GWAS_BASE/logs/summary_152_%j.log" \
        --error="$GWAS_BASE/logs/summary_152_%j.err" \
        --export=ALL,GWAS_BASE="$GWAS_BASE" \
        --wrap="cd '$GWAS_BASE' && $PY_SUMMARY scripts/summarize_plot_A1_B1_152.py && $PY_SUMMARY scripts/summarize_plot_A2_B2_152.py && $PY_SUMMARY scripts/trait_qc_152.py && $PY_SUMMARY scripts/data_driven_discovery.py && $PY_SUMMARY scripts/data_driven_diagnostic.py && $PY_SUMMARY scripts/stat_plot_A2_B2_all_tissues_enrichment.py"
)

echo "[DONE] submitted jobs:"
echo "  make_annot:      $make_jid"
echo "  control_annot:   $ctrl_make_jid"
echo "  ldscore:         $ld_jid"
echo "  control_ldscore: $ctrl_ld_jid"
echo "  A1 h2 152:       $a1_jid"
echo "  B1 h2 152:       $b1_jid"
echo "  A2/B2 CTS 152:   $cts_jid"
echo "  summary:         $summary_jid"
