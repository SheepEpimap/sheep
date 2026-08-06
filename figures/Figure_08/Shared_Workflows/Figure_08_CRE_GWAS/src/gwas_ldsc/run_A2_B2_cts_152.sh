#!/bin/bash
#SBATCH -p smp,low
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -o logs/A2_B2_cts_152_%A_%a.log
#SBATCH -e logs/A2_B2_cts_152_%A_%a.err

set -euo pipefail
: "${GWAS_BASE:?Set GWAS_BASE in the paths configuration}"
: "${LDSC_ROOT:?Set LDSC_ROOT in the paths configuration}"
: "${LDSC_PYTHON:?Set LDSC_PYTHON in the paths configuration}"
: "${LDSC_REF:?Set LDSC_REF in the paths configuration}"
BASE="$GWAS_BASE"
mkdir -p "$BASE/logs"
cd "$BASE"
LDSC="$LDSC_ROOT"
PY="$LDSC_PYTHON"
REF="$LDSC_REF"
MANIFEST="${MANIFEST:-$BASE/manifests/A2_B2_cts_manifest_152.tsv}"

TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
OFFSET=${OFFSET:-0}
GLOBAL_TASK=$((TASK_ID + OFFSET))
echo "[INFO] array_task=$TASK_ID"
echo "[INFO] offset=$OFFSET"
echo "[INFO] global_task=$GLOBAL_TASK"
MAP=$(awk -v n="$((GLOBAL_TASK+2))" 'BEGIN{FS="\t"; OFS="|"} NR==n{print $1,$2,$3,$4,$5}' "$MANIFEST" | tr -d '\r')
[[ -n "$MAP" ]] || { echo "[ERROR] no manifest row for global_task=$GLOBAL_TASK in $MANIFEST" >&2; exit 1; }
IFS='|' read -r LDCTS_NAME TRAIT_ID SUMSTATS LDCTS_FILE OUTDIR <<< "$MAP"
[[ -n "$LDCTS_NAME" && -n "$TRAIT_ID" && -n "$SUMSTATS" && -n "$LDCTS_FILE" && -n "$OUTDIR" ]] || { echo "[ERROR] incomplete manifest row: $MAP" >&2; exit 1; }

mkdir -p "$OUTDIR"
OUT="$OUTDIR/${TRAIT_ID}_${LDCTS_NAME}_cts"
rm -f "$OUT".*

echo "ldcts_name=$LDCTS_NAME"
echo "trait_id=$TRAIT_ID"
echo "ldcts_file=$LDCTS_FILE"
echo "start=$(date '+%F %T')"

"$PY" "$LDSC/ldsc.py" \
    --h2-cts "$SUMSTATS" \
    --ref-ld-chr "$REF/baselineLD/baselineLD." \
    --ref-ld-chr-cts "$LDCTS_FILE" \
    --w-ld-chr "$REF/weights/weights.hm3_noMHC." \
    --out "$OUT"

echo "done=$(date '+%F %T')"
