#!/bin/bash
#SBATCH -p smp,low
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -o logs/make_annot_%A_%a.log
#SBATCH -e logs/make_annot_%A_%a.err

set -euo pipefail

: "${GWAS_BASE:?Set GWAS_BASE in the paths configuration}"
: "${LDSC_REF:?Set LDSC_REF in the paths configuration}"
BASE="$GWAS_BASE"
mkdir -p "$BASE/logs"
cd "$BASE"

MANIFEST="${MANIFEST:-$BASE/manifests/make_annot_jobs.tsv}"
TASK=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
OFFSET=${OFFSET:-0}
GLOBAL_TASK=$((TASK + OFFSET))

echo "[INFO] array_task=$TASK"
echo "[INFO] offset=$OFFSET"
echo "[INFO] global_task=$GLOBAL_TASK"

LINE=$(awk -v n="$((GLOBAL_TASK+2))" 'BEGIN{FS="\t"} NR==n{print}' "$MANIFEST")
[[ -n "$LINE" ]] || { echo "[ERROR] no manifest row for global_task=$GLOBAL_TASK in $MANIFEST" >&2; exit 1; }

IFS=$'\t' read -r ANNOT_ID GROUP BED CHR OUT LD_PREFIX <<< "$LINE"
[[ -n "$ANNOT_ID" && -n "$BED" && -n "$CHR" && -n "$OUT" ]] || { echo "[ERROR] incomplete manifest row: $LINE" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"

echo "[INFO] annotation_id=$ANNOT_ID"
echo "[INFO] group=$GROUP"
echo "[INFO] chr=$CHR"
echo "[INFO] bed=$BED"
echo "[INFO] out=$OUT"
echo "[INFO] start=$(date '+%F %T')"

BIM="$LDSC_REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}.bim"

if [[ ! -s "$BIM" ]]; then
    echo "[ERROR] missing BIM file: $BIM" >&2
    exit 1
fi

EXPECTED_LINES=$(( $(wc -l < "$BIM") + 1 ))

if [[ -s "$OUT" ]]; then
    ACTUAL_LINES=$( (zcat "$OUT" 2>/dev/null || true) | wc -l )

    if [[ "$ACTUAL_LINES" -eq "$EXPECTED_LINES" ]]; then
        echo "[SKIP] existing valid annot: $OUT lines=$ACTUAL_LINES"
        exit 0
    else
        echo "[WARN] existing invalid annot: $OUT lines=$ACTUAL_LINES expected=$EXPECTED_LINES; regenerating" >&2
        rm -f "$OUT"
    fi
fi

"$BASE/scripts/make_thin_annot_python_compatible.sh" "$BED" "$CHR" "$OUT"

echo "[DONE] $(date '+%F %T')"
