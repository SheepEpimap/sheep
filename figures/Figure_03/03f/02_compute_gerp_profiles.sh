#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# =========================
# =========================
HITS_DIR="/data/home/sczd644/run/zsw_chrombpnet/track/hits"
BIGWIG="/vol2/zhuangziyang/4_compair/conservation_score/gerp_conservation_scores.ovis_aries.ARS-UI_Ramb_v3.0.bw"
OUTDIR="/data/home/sczd644/run/zsw_chrombpnet/phylop_gerp_all"
WINDOW="${1:-25}"                 #   ±25
KEEP_TMP="${KEEP_TMP:-0}"         # 0= file;1=

mkdir -p "$OUTDIR/results"

# =========================
# =========================
shopt -s nullglob
BEDS=( "$HITS_DIR"/*_hits_tf.bed )
if [[ ${#BEDS[@]} -eq 0 ]]; then
  echo "[ERROR] No *_hits_tf.bed found in: $HITS_DIR" >&2
  exit 1
fi

# =========================
#    output:motif \t n_hits
# =========================
STATS_TSV="$OUTDIR/motif_counts.tsv"
echo "[INFO] build motif counts: $STATS_TSV"
awk 'BEGIN{FS=OFS="\t"} {cnt[$4]++} END{for(m in cnt) print m, cnt[m]}' "${BEDS[@]}" \
  | sort -k1,1 > "$STATS_TSV"

MOTIFS_TXT="$OUTDIR/motifs.txt"
cut -f1 "$STATS_TSV" > "$MOTIFS_TXT"
echo "[INFO] total motifs: $(wc -l < "$MOTIFS_TXT")"

# =========================
# =========================
while IFS=$'\t' read -r MOTIF NHITS; do
  [[ -z "${MOTIF:-}" ]] && continue

  SAFE="$(echo "$MOTIF" | sed 's/[^A-Za-z0-9._-]/_/g')"
  MDIR="$OUTDIR/results/$SAFE"
  mkdir -p "$MDIR"

  echo "[INFO] motif=$MOTIF  n_hits=$NHITS"

  HITS_BED="$MDIR/${SAFE}.hits.bed"
  WIN_BED="$MDIR/${SAFE}.win.w${WINDOW}.bed"
  OUT_TSV="$MDIR/${SAFE}.gerp_by_offset.w${WINDOW}.tsv"
  META_TSV="$MDIR/${SAFE}.meta.tsv"

awk -v m="$MOTIF" 'BEGIN{FS=OFS="\t"} $4==m{print $1,$2,$3,$6}' "${BEDS[@]}" > "$HITS_BED"


  if [[ ! -s "$HITS_BED" ]]; then
    echo "[WARN] empty hits: $MOTIF"
    continue
  fi

  awk -v W="$WINDOW" 'BEGIN{FS=OFS="\t"; id=0}
  {
    chr=$1; st=$2; en=$3; strand=$4;
    mid=int((st+en)/2);
    id++;
    for(off=-W; off<=W; off++){
      if(strand=="-"){ s=mid-off } else { s=mid+off }
      if(s<0) s=0;
      e=s+1;
      name="h"id"_"off;
      print chr, s, e, name;
    }
  }' "$HITS_BED" > "$WIN_BED"


  bigWigAverageOverBed "$BIGWIG" "$WIN_BED" stdout | \
  awk -v W="$WINDOW" -v N="$NHITS" 'BEGIN{FS=OFS="\t"}
    {
      # name: h<id>_<off>
      split($1, a, "_");
      off=a[2];

      covered=$3+0;
      mean=$6;

      if(covered>0 && mean!="nan" && mean!="NaN"){
        sum[off]+=mean;
        n[off]+=1;
      }
    }
    END{
      print "off","mean_gerp","n_cov","n_total","cov_rate";
      for(off=-W; off<=W; off++){
        if(n[off]>0){
          m=sum[off]/n[off];
          rate=n[off]/N;
          print off, m, n[off], N, rate;
        } else {
          print off, "NA", 0, N, 0;
        }
      }
    }' > "$OUT_TSV"

  {
    echo -e "key\tvalue"
    echo -e "motif\t$MOTIF"
    echo -e "safe\t$SAFE"
    echo -e "window\t$WINDOW"
    echo -e "n_hits\t$NHITS"
  } > "$META_TSV"

  if [[ "$KEEP_TMP" == "0" ]]; then
    rm -f "$WIN_BED"
    # rm -f "$HITS_BED"
  fi

done < "$STATS_TSV"

echo "[DONE] all motif GERP results in: $OUTDIR/results"
