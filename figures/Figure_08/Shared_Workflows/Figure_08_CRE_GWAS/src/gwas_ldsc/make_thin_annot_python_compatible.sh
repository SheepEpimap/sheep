#!/bin/bash
set -euo pipefail

BED=${1:?Usage: make_thin_annot_python_compatible.sh BED CHR OUT}
CHR=${2:?Usage: make_thin_annot_python_compatible.sh BED CHR OUT}
OUT=${3:?Usage: make_thin_annot_python_compatible.sh BED CHR OUT}

: "${LDSC_REF:?Set LDSC_REF in the paths configuration}"
BIM="$LDSC_REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}.bim"

awk -v c="$CHR" '$1==c {print "chr"$1"\t"$4"\t"($4+1)}' "$BIM" \
  | bedtools intersect -a - -b "$BED" -c \
  | awk 'BEGIN{print "ANNOT"} {print ($4>0)?1:0}' \
  | gzip -c > "$OUT"

n_bim=$(awk -v c="$CHR" '$1==c {n++} END{print n+0}' "$BIM")
n_annot=$(zcat "$OUT" | tail -n +2 | wc -l)

if [[ "$n_bim" -ne "$n_annot" ]]; then
    echo "[ERROR] row mismatch: BIM=$n_bim annot=$n_annot chr=$CHR out=$OUT" >&2
    exit 1
fi
