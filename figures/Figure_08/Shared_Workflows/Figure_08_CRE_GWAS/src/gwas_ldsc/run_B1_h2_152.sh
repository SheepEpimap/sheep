#!/bin/bash
#SBATCH -p smp,low
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -o logs/B1_h2_152_%A_%a.log
#SBATCH -e logs/B1_h2_152_%A_%a.err

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
MANIFEST="${MANIFEST:-$BASE/manifests/B1_h2_manifest_152.tsv}"
TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
OFFSET=${OFFSET:-0}
GLOBAL_TASK=$((TASK_ID + OFFSET))
echo "[INFO] array_task=$TASK_ID"
echo "[INFO] offset=$OFFSET"
echo "[INFO] global_task=$GLOBAL_TASK"
MAP=$(awk -v n="$((GLOBAL_TASK+2))" 'BEGIN{FS="\t"; OFS="|"} NR==n{print $1,$2,$3,$4,$5,$6}' "$MANIFEST" | tr -d '\r')
[[ -n "$MAP" ]] || { echo "[ERROR] no manifest row for global_task=$GLOBAL_TASK in $MANIFEST" >&2; exit 1; }
IFS='|' read -r ANNOT_ID CRE_CLASS TRAIT_ID SUMSTATS ANNOT_PREFIX OUTDIR <<< "$MAP"
[[ -n "$ANNOT_ID" && -n "$TRAIT_ID" && -n "$SUMSTATS" && -n "$ANNOT_PREFIX" && -n "$OUTDIR" ]] || { echo "[ERROR] incomplete manifest row: $MAP" >&2; exit 1; }
mkdir -p "$OUTDIR"
OUT="$OUTDIR/${TRAIT_ID}_${ANNOT_ID}_baseline_custom"
rm -f "$OUT".*
echo "annotation_id=$ANNOT_ID"; echo "cre_class=$CRE_CLASS"; echo "trait_id=$TRAIT_ID"; echo "start=$(date '+%F %T')"
"$PY" "$LDSC/ldsc.py" --h2 "$SUMSTATS" --ref-ld-chr "$REF/baselineLD/baselineLD.,$ANNOT_PREFIX" --w-ld-chr "$REF/weights/weights.hm3_noMHC." --frqfile-chr "$REF/1000G_Phase3_frq/1000G.EUR.QC." --overlap-annot --out "$OUT"
echo "done=$(date '+%F %T')"
