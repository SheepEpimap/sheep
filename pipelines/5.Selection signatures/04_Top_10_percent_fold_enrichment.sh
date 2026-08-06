####0. Generate the genome-size file####
####/vol2/mengzhu/genome/reference####
awk 'BEGIN{OFS="\t"}
{
  chr=$1
  if(chr !~ /^chr/) chr="chr"chr
  print chr,0,$2,"ws"
}' sheep1.size \
| sort -k1,1 -k2,2n \
| gzip -c > sheep_workspace_allchr.bed.gz

####1. Extract the top 10% of selection-signal regions####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature####
in="moderneurope_modernasia.10000_10000.windowed.weir.fst.filter.t0.01.annotation"
out="fst.europe_asia_10000.decile.top10pct.tsv"

n=$(( $(wc -l < "$in") - 1 ))      # Number of data rows (excluding the header)
k=$(( (n + 9) / 10 ))              # Top 10%: ceil(n/10)
(( k < 1 )) && k=1                 # Keep at least one row (for very small files)

{
  head -n 1 "$in"
  tail -n +2 "$in" | sort -t $'\t' -k5,5gr | head -n "$k"
} > "$out"

####2. Generate fst.ancient_CEA.top10pct.bed.gz (process selection signatures)####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature####
awk 'BEGIN{FS=OFS="\t"}
NR==1{next}
{
  chr=$1;
  if(chr !~ /^chr/) chr="chr"chr;   # 13 -> chr13 (leave unchanged if already chr13)
  start=$2-1;                       # 1-based -> BED 0-based
  end=$3;
  if(start<0) start=0;
  print chr,start,end,"selection"   # Use one track name in column 4 so intervals are not treated as separate tracks
}' ancientasia_4000y.10000_10000.windowed.weir.fst.filter.t0.01.annotation \
| sort -k1,1 -k2,2n \
| gzip -c > GAT_fst.ancientasia_10000.decile.top10pct_1.bed.gz

####3. Generate GAT_TSR_*_E5.bed.gz files (process TSR enhancers)####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5####
outdir="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature"
mkdir -p "$outdir"

shopt -s nullglob  # If no files match, do not pass the wildcard literally

for inbed in TSR_*_E5_id.bed; do
  state="${inbed%_id.bed}"   # For example: TSR_thymus_E5
  out="$outdir/GAT_${state}.bed.gz"

  awk -v OFS='\t' -v s="$state" '{print $1,$2,$3,s}' "$inbed" \
  | sort -k1,1 -k2,2n \
  | gzip -c > "$out"

  echo "[OK] $inbed -> $out"
done

####4. Run####
####Use GAT to calculate fold enrichment between selection signals and tissue-specific enhancers
vim 02_GAT_europe_asia.sh
#!/bin/bash
set -euo pipefail

SEG="GAT_fst.europe_asia_10000.decile.top10pct.bed.gz"
WS="/vol2/mengzhu/genome/reference/sheep_workspace_allchr.bed.gz"

OUTDIR="GAT_results_ancient_10000_europe_asia"
LOGDIR="GAT_logs_ancient_10000_europe_asia"
mkdir -p "$OUTDIR" "$LOGDIR"

shopt -s nullglob

for ANN in GAT_*.bed.gz; do
  # Skip the segment file itself (it also matches GAT_*.bed.gz)
  [[ "$ANN" == "$SEG" ]] && continue

  base="$(basename "$ANN" .bed.gz)"   # e.g. GAT_TSR_thymus_E5
  tag="${base#GAT_}"                  # e.g. TSR_thymus_E5

  echo "[RUN] annotation=$ANN -> ${OUTDIR}/${tag}.tsv"

  gat-run.py \
    --segment-file="$SEG" \
    --annotation-file="$ANN" \
    --workspace-file="$WS" \
    --counter=nucleotide-overlap \
    --num-samples=10000 \
    --qvalue-method=BH \
    --num-threads=10 \
    --nbuckets=100001 \
    --log="${LOGDIR}/${tag}.log" \
    > "${OUTDIR}/${tag}.europe_asia.tsv"
done
#sbatch -c 10  -p low --mem 20G 02_GAT_europe_asia.sh

###Domestication analysis###
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "selection"}' Demestication_domestic_sheep.bed > Demestication_domestic_sheep.selection.bed
gzip -c Demestication_domestic_sheep.selection.bed >Demestication_domestic_sheep.selection.bed.gz

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "selection"}' Demestication_mouflon.bed > Demestication_mouflon_sheep.selection.bed
gzip -c Demestication_mouflon_sheep.selection.bed >Demestication_mouflon_sheep.selection.bed.gz

vim 02_GAT_Demestication_mouflon.sh
#!/bin/bash
set -euo pipefail

SEG="Demestication_mouflon_sheep.selection.bed.gz"
WS="/vol2/mengzhu/genome/reference/sheep_workspace_allchr.bed.gz"

OUTDIR="GAT_results_Demestication_mouflon"
LOGDIR="GAT_logs_Demestication_mouflon"
mkdir -p "$OUTDIR" "$LOGDIR"

shopt -s nullglob

for ANN in GAT_*.bed.gz; do
  # Skip the segment file itself (it also matches GAT_*.bed.gz)
  [[ "$ANN" == "$SEG" ]] && continue

  base="$(basename "$ANN" .bed.gz)"   # e.g. GAT_TSR_thymus_E5
  tag="${base#GAT_}"                  # e.g. TSR_thymus_E5

  echo "[RUN] annotation=$ANN -> ${OUTDIR}/${tag}.tsv"

  gat-run.py \
    --segment-file="$SEG" \
    --annotation-file="$ANN" \
    --workspace-file="$WS" \
    --counter=nucleotide-overlap \
    --num-samples=10000 \
    --qvalue-method=BH \
    --num-threads=10 \
    --nbuckets=100001 \
    --log="${LOGDIR}/${tag}.log" \
    > "${OUTDIR}/${tag}.Demestication_mouflon.tsv"
done
#sbatch -c 10  -p low --mem 20G 02_GAT_Demestication_mouflon.sh

#Prepare data and generate the plotting input file
awk 'FNR>1 || NR==1' TSR_*_E5.ancient_CEA_1.tsv > All_TSR_E5.ancient_CEA_1.tsv
awk 'FNR>1 || NR==1' TSR_*_E5.ancient_EUR_1.tsv > All_TSR_E5.ancient_EUR_1.tsv
awk 'FNR>1 || NR==1' TSR_*_E5.europe_asia_1.tsv > All_TSR_E5.europe_asia_1.tsv
awk 'FNR>1 || NR==1' TSR_*.Demestication_domestic.tsv > All_TSR_E5.Demestication_domestic.tsv
awk 'BEGIN{
  OFS="\t";
  print "Demestication_domestic_annotation","Demestication_domestic_fold","Demestication_domestic_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.Demestication_domestic.tsv > All_TSR_E5.Demestication_domestic_1.tsv

awk 'FNR>1 || NR==1' TSR_*.Demestication_mouflon.tsv > All_TSR_E5.Demestication_mouflon.tsv
awk 'BEGIN{
  OFS="\t";
  print "Demestication_mouflon_annotation","Demestication_mouflon_fold","Demestication_mouflon_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.Demestication_mouflon.tsv > All_TSR_E5.Demestication_mouflon_1.tsv

awk 'BEGIN{
  OFS="\t";
  print "ancient_EUR_1000_annotation","ancient_EUR_1000_fold","ancient_EUR_1000_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.ancient_EUR_1.tsv > All_TSR_E5.ancient_EUR_2.tsv

awk 'BEGIN{
  OFS="\t";
  print "europe_asia_1000_annotation","europe_asia_1000_fold","europe_asia_1000_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.europe_asia_1.tsv > All_TSR_E5.europe_asia_2.tsv

awk 'BEGIN{
  OFS="\t";
  print "ancient_CEA_1000_annotation","ancient_CEA_1000_fold","ancient_CEA_1000_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.ancient_CEA_1.tsv > All_TSR_E5.ancient_CEA_2.tsv

paste \
  All_TSR_E5.ancient_CEA_2.tsv \
  <(cut -f2-3 All_TSR_E5.ancient_EUR_2.tsv) \
  <(cut -f2-3 All_TSR_E5.europe_asia_2.tsv) \
  <(cut -f2-3 All_TSR_E5.Demestication_domestic_1.tsv) \
> All_TSR_E5_Selection_signatures_1.tsv

