#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

LIST=${1:?Usage: $0 <tissue_list.txt>}

OUTBASE="/data/home/sczd644/run/zsw_chrombpnet/network/TSR"
OUTDIR="${OUTBASE}/ALL_fullChain"
mkdir -p "$OUTDIR"

OUT_EDGE="${OUTDIR}/ALL.edge.fullChain.txt"   #  :id type tissue
OUT_NODE="${OUTDIR}/ALL.node.fullChain.txt"   #  :node1 node2 weight tissue
OUT_TAB="${OUTDIR}/ALL.all_tables.with_header.fullChain.txt"

: > "$OUT_EDGE"
: > "$OUT_NODE"
: > "$OUT_TAB"

first_edge=1
first_node=1
first_tab=1

while IFS= read -r tissue; do
  [[ -z "$tissue" || "$tissue" =~ ^# ]] && continue
  tissue="${tissue//$'\r'/}"

  TDIR="${OUTBASE}/${tissue}/tables"
  EDGE="${TDIR}/${tissue}.edge.fullChain.txt"
  NODE="${TDIR}/${tissue}.node.fullChain.txt"
  TAB="${TDIR}/${tissue}.all_tables.with_header.fullChain.txt"

  if [[ ! -s "$EDGE" || ! -s "$NODE" || ! -s "$TAB" ]]; then
    echo "[WARN] skip ${tissue}: missing final tables in ${TDIR}" >&2
    continue
  fi

  if [[ $first_edge -eq 1 ]]; then
    cat "$EDGE" >> "$OUT_EDGE"
    first_edge=0
  else
    tail -n +2 "$EDGE" >> "$OUT_EDGE"
  fi

  if [[ $first_node -eq 1 ]]; then
    cat "$NODE" >> "$OUT_NODE"
    first_node=0
  else
    tail -n +2 "$NODE" >> "$OUT_NODE"
  fi

  if [[ $first_tab -eq 1 ]]; then
    cat "$TAB" >> "$OUT_TAB"
    first_tab=0
  else
    tail -n +2 "$TAB" >> "$OUT_TAB"
  fi

done < "$LIST"

echo "[DONE] merged:"
echo "  $OUT_EDGE"
echo "  $OUT_NODE"
echo "  $OUT_TAB"


bash 2.sh tissue.txt
