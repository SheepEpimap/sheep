#!/bin/bash
#SBATCH -p smp,low
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -o logs/ldscore_control_%A_%a.log
#SBATCH -e logs/ldscore_control_%A_%a.err

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
MANIFEST="${MANIFEST:-$BASE/manifests/control_ldscore_jobs.tsv}"

TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
OFFSET=${OFFSET:-0}
GLOBAL_TASK=$((TASK_ID + OFFSET))

echo "[INFO] array_task=$TASK_ID"
echo "[INFO] offset=$OFFSET"
echo "[INFO] global_task=$GLOBAL_TASK"

MAP=$(awk -v n="$((GLOBAL_TASK+2))" 'BEGIN{FS="\t"; OFS="|"} NR==n{print $1,$2,$3,$4}' "$MANIFEST")
[[ -n "$MAP" ]] || { echo "[ERROR] no manifest row for global_task=$GLOBAL_TASK in $MANIFEST" >&2; exit 1; }
IFS='|' read -r ANNOT_ID CHROM ANNOT OUT_PREFIX <<< "$MAP"
[[ -n "$ANNOT_ID" && -n "$CHROM" && -n "$ANNOT" && -n "$OUT_PREFIX" ]] || { echo "[ERROR] incomplete manifest row: $MAP" >&2; exit 1; }

mkdir -p "$(dirname "$OUT_PREFIX")"

ln -sf "$ANNOT" "$(dirname "$OUT_PREFIX")/${ANNOT_ID}.${CHROM}.annot.gz"

rm -f "$OUT_PREFIX".log \
      "$OUT_PREFIX".l2.ldscore.gz \
      "$OUT_PREFIX".l2.M \
      "$OUT_PREFIX".l2.M_5_50

echo "annotation_id=$ANNOT_ID"
echo "chrom=$CHROM"
echo "annot=$ANNOT"
echo "out_prefix=$OUT_PREFIX"
echo "start=$(date '+%F %T')"

"$PY" "$LDSC/ldsc.py" \
  --l2 \
  --bfile "$REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHROM" \
  --ld-wind-cm 1 \
  --annot "$ANNOT" \
  --thin-annot \
  --print-snps "$REF/baselineLD_snps/baselineLD.$CHROM.snp" \
  --out "$OUT_PREFIX"

echo "done=$(date '+%F %T')"
