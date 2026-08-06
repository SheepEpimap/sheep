#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# ============================================================
#
# Output directory:
#   /data/home/sczd644/run/zsw_chrombpnet/network/<tissue>/
#     tmp/ logs/
# ============================================================

tissue=${1:?Usage: $0 <tissue>}

OUTBASE="/data/home/sczd644/run/zsw_chrombpnet/network/TSR"
OUTDIR="${OUTBASE}/${tissue}"
TMPDIR="${OUTDIR}/tmp"
EDGEDIR="${OUTDIR}/edges"
TABDIR="${OUTDIR}/tables"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$TMPDIR" "$EDGEDIR" "$TABDIR" "$LOGDIR"

LOG="${LOGDIR}/${tissue}.build_network.log"
exec > >(tee -i "$LOG") 2>&1

echo "[INFO] tissue=${tissue}"
date

# ---------- input ----------
PAIRS="/vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_output_E5_confident.tsv"
E5_BED="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/TSR_${tissue}_E5_id.bed"
E1_BED="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E1/TSR_${tissue}_E1_id.bed"
HITS="/data/home/sczd644/run/zsw_chrombpnet/finemo/${tissue}_finemo/hits_with_tf_names.tsv"
TSS2K="/vol2/mengzhu/genome/part_change_esemb100/TSS_esemble100_colin.bed_up2k.bed"

for cmd in bedtools gawk awk sort cut sed tail head; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "[ERROR] missing command: $cmd" >&2; exit 1; }
done
for f in "$PAIRS" "$E5_BED" "$E1_BED" "$HITS" "$TSS2K"; do
  [[ -s "$f" ]] || { echo "[ERROR] Missing/empty: $f" >&2; exit 1; }
done

DEDUP_HITS=${DEDUP_HITS:-1}

# ============================================================
# ============================================================
E5_3="${TMPDIR}/${tissue}.E5.3col.bed"
E1_3="${TMPDIR}/${tissue}.E1.3col.bed"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$E5_BED" > "$E5_3"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$E1_BED" > "$E1_3"

# ============================================================
#   chr st en enhancer_id corr pval gene dist fdr
# ============================================================
PAIRS_BED="${TMPDIR}/E5_pairs.all.bed"
awk '
  BEGIN{FS="[ \t]+"; OFS="\t"}
  {
    enh=$1
    split(enh,a,":"); chr=a[1]
    split(a[2],b,"-"); st=b[1]; en=b[2]
    corr=$2; pval=$3; gene=$4; dist=$5; fdr=$6
    print chr,st,en,enh,corr,pval,gene,dist,fdr
  }
' "$PAIRS" > "$PAIRS_BED"

# ============================================================
# ============================================================
PAIRS_TISSUE="${TMPDIR}/${tissue}.pairs.in_tissueE5.bed9"
bedtools intersect -wa -a "$PAIRS_BED" -b "$E5_3" > "$PAIRS_TISSUE"

# ============================================================
# ============================================================
ENH_CAND_BED="${TMPDIR}/${tissue}.enhancers.cand.nonLOC.bed4"
awk -v OFS="\t" '$7 !~ /^LOC/ {print $1,$2,$3,$4}' "$PAIRS_TISSUE" | awk '!seen[$0]++' > "$ENH_CAND_BED"

# ============================================================
# Step 5) Fi-NeMo hits -> bed6: chr start end motif hit_coef_global strand
# ============================================================
HITS_BED="${TMPDIR}/${tissue}.hits.bed6"
awk -F'\t' -v OFS='\t' 'NR>1 {print $1,$2,$3,$6,$8,$13}' "$HITS" > "$HITS_BED"

HITS_USE="$HITS_BED"
if [[ "$DEDUP_HITS" -eq 1 ]]; then
  HITS_DEDUP="${TMPDIR}/${tissue}.hits.dedup.max.bed6"
  gawk -F'\t' -v OFS='\t' '
    {
      key=$1 FS $2 FS $3 FS $4 FS $6
      w=$5+0
      if(!(key in mx) || w>mx[key]){
        mx[key]=w
        line[key]=$1 OFS $2 OFS $3 OFS $4 OFS w OFS $6
      }
    }
    END{for(k in line) print line[k]}
  ' "$HITS_BED" | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n -k4,4 -k6,6 > "$HITS_DEDUP"
  HITS_USE="$HITS_DEDUP"
fi

# ============================================================
# ============================================================
HIT_X_ENH="${TMPDIR}/${tissue}.hit_x_enh.tsv"
bedtools intersect -wa -wb -a "$HITS_USE" -b "$ENH_CAND_BED" > "$HIT_X_ENH"

# ============================================================
# ============================================================
MOTIF_ENH_EDGE_ALL="${EDGEDIR}/${tissue}.motif_enhancer.edges.tsv"
MOTIF_ENH_QC_ALL="${EDGEDIR}/${tissue}.motif_enhancer.qc.tsv"
ENH_WITH_HIT="${TMPDIR}/${tissue}.enhancers.with_motifHit.txt"

gawk -v OFS="\t" -v tissue="$tissue" '
  {
    motif=$4
    w=$5+0
    enh=$10
    key=motif SUBSEP enh
    sum[key]+=w
    if(!(key in mx) || w>mx[key]) mx[key]=w
    n[key]+=1
  }
  END{
    for(k in sum){
      split(k,a,SUBSEP)
      motif=a[1]; enh=a[2]
      print motif, enh, tissue, sum[k], mx[k], n[k]
    }
  }
' "$HIT_X_ENH" \
| LC_ALL=C sort -t $'\t' -k1,1 -k2,2 \
> "${TMPDIR}/${tissue}.motif_enh.agg.tsv"

{
  echo -e "source\ttarget\ttissue\tweight_sum\tweight_max\tn_hits"
  cat "${TMPDIR}/${tissue}.motif_enh.agg.tsv"
} > "$MOTIF_ENH_QC_ALL"

{
  echo -e "source\ttarget\ttissue\tedge_type\tweight_sum"
  awk -F'\t' -v OFS='\t' 'NR>=1 {print $1,$2,$3,"motif-enhancer",$4}' "${TMPDIR}/${tissue}.motif_enh.agg.tsv"
} > "$MOTIF_ENH_EDGE_ALL"

cut -f2 "$MOTIF_ENH_EDGE_ALL" | tail -n +2 | LC_ALL=C sort -u > "$ENH_WITH_HIT"

# ============================================================
# ============================================================
ENH_GENE_EDGE_ALL="${EDGEDIR}/${tissue}.enhancer_gene.edges.tsv"
{
  echo -e "source\ttarget\ttissue\tedge_type\tweight_corr\tpvalue\tdistance\tfdr"
  awk -v OFS="\t" -v tissue="$tissue" '
    NR==FNR {keep[$1]=1; next}
    ($4 in keep) && ($7 !~ /^LOC/) {print $4,$7,tissue,"enhancer-gene",$5,$6,$8,$9}
  ' "$ENH_WITH_HIT" "$PAIRS_TISSUE" | awk '!seen[$0]++'
} > "$ENH_GENE_EDGE_ALL"

GENES_EG="${TMPDIR}/${tissue}.genes.from_EG.txt"
cut -f2 "$ENH_GENE_EDGE_ALL" | tail -n +2 | LC_ALL=C sort -u > "$GENES_EG"

# ============================================================
# ============================================================
PROMO_CAND="${TMPDIR}/${tissue}.promoter_from_TSS2k.bed6"
# TSS2K: chr start end gene1 gene2 strand biotype
awk -v OFS="\t" '
  NR==FNR {keep[$1]=1; next}
  (($4 in keep) || ($5 in keep)) {
    gene=$4; if(!(gene in keep)) gene=$5
    print $1,$2,$3,gene,$6,$7
  }
' "$GENES_EG" "$TSS2K" | awk '!seen[$0]++' > "$PROMO_CAND"

# ============================================================
# ============================================================

PROMO_INT="${TMPDIR}/${tissue}.promoter.E1_within2kb.bed7"

bedtools intersect -wa -wb \
  -a "$E1_3" -b "$PROMO_CAND" \
| awk -v OFS="\t" '
    {
      # A: E1_3 => $1 $2 $3
      chr=$1; st=$2; en=$3
      promoter=chr ":" st "-" en   # promoter   E1  ( )

      gene=$7; strand=$8; biotype=$9

      print chr,st,en,promoter,gene,strand,biotype
    }
  ' | awk '!seen[$0]++' > "$PROMO_INT"

PROMO_GENE_EDGE_ALL="${EDGEDIR}/${tissue}.promoter_gene.edges.tsv"
{
  echo -e "source\ttarget\ttissue\tedge_type\tweight"
  awk -F'\t' -v OFS="\t" -v tissue="$tissue" '
    {print $4,$5,tissue,"promoter-gene",1}
  ' "$PROMO_INT" | awk '!seen[$0]++'
} > "$PROMO_GENE_EDGE_ALL"

# ============================================================
# ============================================================
PROMO_GENES="${TMPDIR}/${tissue}.genes.with_promoter.txt"
cut -f2 "$PROMO_GENE_EDGE_ALL" | tail -n +2 | LC_ALL=C sort -u > "$PROMO_GENES"

ENH_GENE_EDGE_PF1="${EDGEDIR}/${tissue}.enhancer_gene.promoterFiltered.tsv"
awk -F'\t' -v OFS='\t' '
  NR==FNR {keep[$1]=1; next}
  FNR==1 {print; next}
  ($2 in keep) {print}
' "$PROMO_GENES" "$ENH_GENE_EDGE_ALL" > "$ENH_GENE_EDGE_PF1"

ENH_PF1="${TMPDIR}/${tissue}.enhancers.after_promoterFilter.txt"
cut -f1 "$ENH_GENE_EDGE_PF1" | tail -n +2 | LC_ALL=C sort -u > "$ENH_PF1"

MOTIF_ENH_EDGE_PF1="${EDGEDIR}/${tissue}.motif_enhancer.promoterFiltered.tsv"
awk -F'\t' -v OFS='\t' '
  NR==FNR {keep[$1]=1; next}
  FNR==1 {print; next}
  ($2 in keep) {print}
' "$ENH_PF1" "$MOTIF_ENH_EDGE_ALL" > "$MOTIF_ENH_EDGE_PF1"

PROMO_GENE_EDGE_PF1="${EDGEDIR}/${tissue}.promoter_gene.promoterFiltered.tsv"
awk -F'\t' -v OFS='\t' '
  NR==FNR {keep[$1]=1; next}
  FNR==1 {print; next}
  ($2 in keep) {print}
' "$PROMO_GENES" "$PROMO_GENE_EDGE_ALL" > "$PROMO_GENE_EDGE_PF1"

ALL_EDGES_PF1="${EDGEDIR}/${tissue}.all_edges.promoterFiltered.tsv"
{
  echo -e "source\ttarget\ttissue\tedge_type\tweight\tpvalue_or_na\tdistance_or_na\tfdr_or_na"
  awk -F'\t' -v OFS='\t' 'NR>1 {print $1,$2,$3,$4,$5,"NA","NA","NA"}' "$MOTIF_ENH_EDGE_PF1"
  awk -F'\t' -v OFS='\t' 'NR>1 {print $1,$2,$3,$4,$5,$6,$7,$8}' "$ENH_GENE_EDGE_PF1"
  awk -F'\t' -v OFS='\t' 'NR>1 {print $1,$2,$3,$4,$5,"NA","NA","NA"}' "$PROMO_GENE_EDGE_PF1"
} > "$ALL_EDGES_PF1"

# ============================================================
# ============================================================
TABLE_FULL="${TABDIR}/${tissue}.all_tables.with_header.fullChain.txt"

gawk -F'\t' -v OFS='\t' -v tissue="$tissue" '
  ARGIND==1 && FNR==1 {next}
  ARGIND==1 {
    motif=$1; enh=$2; me=$5+0
    key=enh SUBSEP motif
    me_w[key]=me
    me_list[enh]=(enh in me_list ? me_list[enh] ";" motif : motif)
    next
  }
  ARGIND==2 && FNR==1 {next}
  ARGIND==2 {
    enh=$1; gene=$2; eg=$5+0
    key=enh SUBSEP gene
    eg_w[key]=eg
    eg_list[gene]=(gene in eg_list ? eg_list[gene] ";" enh : enh)
    next
  }
  ARGIND==3 && FNR==1 {
    print "Tissue","Motif","Enhancer","Promoter","Gene","Motif_Enhancer","Enhancer_Gene","Promoter_Gene"
    next
  }
  ARGIND==3 {
    promoter=$1; gene=$2; pg=$5+0
    if(!(gene in eg_list)) next
    n=split(eg_list[gene], enhs, ";")
    for(i=1;i<=n;i++){
      enh=enhs[i]
      if(!(enh in me_list)) continue
      eg=eg_w[enh SUBSEP gene]+0
      m=split(me_list[enh], motifs, ";")
      for(j=1;j<=m;j++){
        motif=motifs[j]
        me=me_w[enh SUBSEP motif]+0
        print tissue,motif,enh,promoter,gene,me,eg,pg
      }
    }
  }
' "$MOTIF_ENH_EDGE_PF1" "$ENH_GENE_EDGE_PF1" "$PROMO_GENE_EDGE_PF1" > "$TABLE_FULL"

# ============================================================
# ============================================================
EDGE_FULL="${TABDIR}/${tissue}.edge.fullChain.txt"   #  :id type tissue
NODE_FULL="${TABDIR}/${tissue}.node.fullChain.txt"   #  :node1 node2 weight tissue

{
  echo -e "id\ttype\ttissue"
  awk -F'\t' -v OFS='\t' 'NR>1 {print $2,"Motif",$1; print $3,"Enhancer",$1; print $4,"Promoter",$1; print $5,"Gene",$1}' "$TABLE_FULL" \
  | awk '!seen[$0]++'
} > "$EDGE_FULL"

{
  echo -e "node1\tnode2\tweight\ttissue"
  awk -F'\t' -v OFS='\t' '
    NR==1 {next}
    {
      tissue=$1
      motif=$2; enh=$3; prom=$4; gene=$5
      me=$6+0; eg=$7+0; pg=$8+0
      k1=motif SUBSEP enh SUBSEP me SUBSEP tissue
      k2=enh   SUBSEP gene SUBSEP eg SUBSEP tissue
      k3=prom  SUBSEP gene SUBSEP pg SUBSEP tissue
      if(!(k1 in seen)){seen[k1]=1; print motif,enh,me,tissue}
      if(!(k2 in seen)){seen[k2]=1; print enh,gene,eg,tissue}
      if(!(k3 in seen)){seen[k3]=1; print prom,gene,pg,tissue}
    }
  ' "$TABLE_FULL"
} > "$NODE_FULL"

echo "[DONE] ${tissue}"
echo "  edges(all):          $EDGEDIR"
echo "  fullChain table:     $TABLE_FULL"
echo "  fullChain edge(node):$EDGE_FULL"
echo "  fullChain node(link):$NODE_FULL"
date
