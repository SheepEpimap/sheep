#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# =========================================================
#
# input:
#   /data/home/sczd644/run/zsw_chrombpnet/track/hits/*_hits_tf.bed
#
# output:
# =========================================================

HITS_DIR="/data/home/sczd644/run/zsw_chrombpnet/track/hits/"
OUTDIR="/data/home/sczd644/run/zsw_chrombpnet/phylop"
mkdir -p "$OUTDIR"

OUT_UNION_GZ="${OUTDIR}/union_all_tissues_hits_tf.tsv.gz"
OUT_COUNTS="${OUTDIR}/union_motif_counts.tsv"

TMP_UNION="$(mktemp "${OUTDIR}/.tmp.union.XXXXXX.tsv")"
trap 'rm -f "$TMP_UNION"' EXIT

shopt -s nullglob
files=( "${HITS_DIR}"/*_hits_tf.bed )
if [[ ${#files[@]} -eq 0 ]]; then
  echo "[ERROR] No *_hits_tf.bed found in $HITS_DIR" >&2
  exit 1
fi

: > "$TMP_UNION"

for f in "${files[@]}"; do
  base="$(basename "$f")"
  tissue="${base%_hits_tf.bed}"

  awk -v tissue="$tissue" '
    BEGIN{FS="[ \t]+"; OFS="\t"}
    NF>=6{
      chr=$1; s=$2; e=$3; m=$4; score=$5; strand=$6;

      if(chr=="" || s+0<0 || e+0<=s+0) next;

      print chr, s, e, m, score, strand, tissue
    }
  ' "$f" >> "$TMP_UNION"
done

gzip -c "$TMP_UNION" > "$OUT_UNION_GZ"

echo "[DONE] Union saved: $OUT_UNION_GZ" >&2
echo "[INFO] Total hits (all tissues): $(zcat "$OUT_UNION_GZ" | wc -l)" >&2

zcat "$OUT_UNION_GZ" \
  | awk 'BEGIN{FS=OFS="\t"} {cnt[$4]++} END{for(m in cnt) print m, cnt[m]}' \
  | sort -k2,2nr > "$OUT_COUNTS"

echo "[DONE] Motif counts saved: $OUT_COUNTS" >&2
