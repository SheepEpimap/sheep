#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail
export LC_ALL=C

############################################
# inputfile
############################################
GENE_COUNT_FILE="H3K27ac_output_E5_gene_count.txt"
EXPR_FILE="/vol2/mengzhu/snakemake_sheep/expressiondir/rna/all_tisssues_expression_tpm.csv"
TAU_FILE="/vol2/mengzhu/snakemake_sheep/expressiondir/all_tisssues_expression_tpm.median.tau.csv"
CONS_FILE="/vol2/mengzhu/genome/conservation_1/sheep_to_human_conservation_last.txt"

############################################
############################################
if [ $# -ne 4 ]; then
    echo " : $0 < > < > < > < >"
    echo " : $0 80 20 30 5"
    exit 1
fi

############################################
############################################
HIGH_THRESHOLD="$1"
MEDIUM_LOW="$2"
MEDIUM_HIGH="$3"
LOW_THRESHOLD="$4"

echo " grouping :"
echo " : >= $HIGH_THRESHOLD"
echo " : $MEDIUM_LOW ~ $MEDIUM_HIGH"
echo " : <= $LOW_THRESHOLD"

############################################
############################################
for f in "$GENE_COUNT_FILE" "$EXPR_FILE" "$TAU_FILE" "$CONS_FILE"; do
    if [ ! -f "$f" ]; then
        echo "[ERROR] file : $f" >&2
        exit 1
    fi
done

############################################
############################################
TMPDIR="tmp_fix_join_${HIGH_THRESHOLD}_${MEDIUM_LOW}_${MEDIUM_HIGH}_${LOW_THRESHOLD}"
mkdir -p "$TMPDIR"

############################################
############################################
awk '
BEGIN{FS="[[:space:]]+"; OFS="\t"}
{
    sub(/\r$/, "")
}
NF{
    $1=$1
    print
}' "$GENE_COUNT_FILE" > "$TMPDIR/gene_count.norm.tsv"

############################################
############################################
awk '
BEGIN{FS="[[:space:]]+"; OFS="\t"}
{
    sub(/\r$/, "")
}
NR==1 && $1=="ID" {next}
NF{
    $1=$1
    print
}' "$EXPR_FILE" > "$TMPDIR/expression.noheader.tsv"

############################################
############################################
awk '
BEGIN{FS="[[:space:]]+"; OFS="\t"}
{
    sub(/\r$/, "")
}
NF{
    $1=$1
    print
}' "$TAU_FILE" > "$TMPDIR/tau.norm.tsv"

############################################
#    1 gene
#    2 human_id
#    3 humantochicken
#    4 conservation
#    5 q
############################################
awk '
BEGIN{FS="[[:space:]]+"; OFS="\t"}
{
    sub(/\r$/, "")
}
NF{
    $1=$1
    gene=$1
    cons=$4+0
    if (!(gene in best_cons) || cons > best_cons[gene]) {
        best_cons[gene]=cons
        best_line[gene]=$0
    }
}
END{
    for (g in best_line) {
        print best_line[g]
    }
}' "$CONS_FILE" \
| sort -t $'\t' -k1,1 > "$TMPDIR/conservation.dedup.tsv"

############################################
############################################
sort -t $'\t' -k1,1 "$TMPDIR/gene_count.norm.tsv" > "$TMPDIR/gene_count.sorted.tsv"
sort -t $'\t' -k1,1 "$TMPDIR/expression.noheader.tsv" > "$TMPDIR/expression.sorted.tsv"
sort -t $'\t' -k1,1 "$TMPDIR/tau.norm.tsv" > "$TMPDIR/tau.sorted.tsv"

############################################
############################################
awk -F '\t' -v OFS='\t' -v threshold="$HIGH_THRESHOLD" '
$2 >= threshold {print}
' "$TMPDIR/gene_count.norm.tsv" > "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}.txt"

join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}.txt") \
    "$TMPDIR/expression.sorted.tsv" \
    > "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_expression.txt"

cut -f 8-93 "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_expression.txt" \
| awk '
BEGIN{FS=OFS="\t"}
{
    sum=0
    for(i=1;i<=NF;i++) sum += $i
    print sum, sum/86
}' > high_1.txt

cut -f 1-7 "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_expression.txt" > high_3.txt

paste high_3.txt high_1.txt \
| awk 'BEGIN{FS=OFS="\t"} {$NF=$NF; $(NF+1)="Top"; print}' \
> "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_expression_average.txt"

############################################
############################################
awk -F '\t' -v OFS='\t' -v low="$MEDIUM_LOW" -v high="$MEDIUM_HIGH" '
$2 >= low && $2 <= high {print}
' "$TMPDIR/gene_count.norm.tsv" > "H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}.txt"

join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 "H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}.txt") \
    "$TMPDIR/expression.sorted.tsv" \
    > "H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}_expression.txt"

cut -f 8-93 "H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}_expression.txt" \
| awk '
BEGIN{FS=OFS="\t"}
{
    sum=0
    for(i=1;i<=NF;i++) sum += $i
    print sum, sum/86
}' > medium_1.txt

cut -f 1-7 "H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}_expression.txt" > medium_3.txt

paste medium_3.txt medium_1.txt \
| awk 'BEGIN{FS=OFS="\t"} {$NF=$NF; $(NF+1)="Middle"; print}' \
> "H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}_expression_average.txt"

############################################
############################################
awk -F '\t' -v OFS='\t' -v threshold="$LOW_THRESHOLD" '
$2 <= threshold {print}
' "$TMPDIR/gene_count.norm.tsv" > "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}.txt"

join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}.txt") \
    "$TMPDIR/expression.sorted.tsv" \
    > "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_expression.txt"

cut -f 8-93 "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_expression.txt" \
| awk '
BEGIN{FS=OFS="\t"}
{
    sum=0
    for(i=1;i<=NF;i++) sum += $i
    print sum, sum/86
}' > low_1.txt

cut -f 1-7 "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_expression.txt" > low_3.txt

paste low_3.txt low_1.txt \
| awk 'BEGIN{FS=OFS="\t"} {$NF=$NF; $(NF+1)="Bottom"; print}' \
> "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_expression_average.txt"

############################################
############################################
join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_expression_average.txt") \
    "$TMPDIR/tau.sorted.tsv" \
    > high_final.txt

join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 "H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}_expression_average.txt") \
    "$TMPDIR/tau.sorted.tsv" \
    > medium_final.txt

join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_expression_average.txt") \
    "$TMPDIR/tau.sorted.tsv" \
    > low_final.txt

cat high_final.txt medium_final.txt low_final.txt > H3K27ac_output_E5_gene_count_top005_expression_all_1.txt

############################################
############################################
join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 high_final.txt) \
    "$TMPDIR/conservation.dedup.tsv" \
    > "H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_conservation.txt"

join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 medium_final.txt) \
    "$TMPDIR/conservation.dedup.tsv" \
    > "H3K27ac_output_E5_gene_count_median_${MEDIUM_LOW}_${MEDIUM_HIGH}_conservation.txt"

join -t $'\t' -1 1 -2 1 \
    <(sort -t $'\t' -k1,1 low_final.txt) \
    "$TMPDIR/conservation.dedup.tsv" \
    > "H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_conservation.txt"

cat \
"H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_conservation.txt" \
"H3K27ac_output_E5_gene_count_median_${MEDIUM_LOW}_${MEDIUM_HIGH}_conservation.txt" \
"H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_conservation.txt" \
> H3K27ac_output_E5_gene_count_top005_expression_conservation_1.txt

############################################
############################################
awk 'BEGIN{FS=OFS="\t"} { $9=sprintf("%.2f", $9); $11=sprintf("%.4f", $11); print $0 }' \
H3K27ac_output_E5_gene_count_top005_expression_all_1.txt \
> H3K27ac_output_E5_gene_count_top005_expression_all_2.txt

awk 'BEGIN{FS=OFS="\t"} { $9=sprintf("%.2f", $9); $11=sprintf("%.4f", $11); $14=sprintf("%.2f", $14); print $0 }' \
H3K27ac_output_E5_gene_count_top005_expression_conservation_1.txt \
> H3K27ac_output_E5_gene_count_top005_expression_conservation_2.txt

############################################
############################################
{
    echo "=====  check ====="
    echo "Top  grouping: $(wc -l < H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}.txt)"
    echo "Top  : $(wc -l < H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_expression.txt)"
    echo "Top tau : $(wc -l < high_final.txt)"
    echo "Top conservation : $(wc -l < H3K27ac_output_E5_gene_count_top005_${HIGH_THRESHOLD}_conservation.txt)"
    echo
    echo "Middle  grouping: $(wc -l < H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}.txt)"
    echo "Middle  : $(wc -l < H3K27ac_output_E5_gene_count_medium_${MEDIUM_LOW}_${MEDIUM_HIGH}_expression.txt)"
    echo "Middle tau : $(wc -l < medium_final.txt)"
    echo "Middle conservation : $(wc -l < H3K27ac_output_E5_gene_count_median_${MEDIUM_LOW}_${MEDIUM_HIGH}_conservation.txt)"
    echo
    echo "Bottom  grouping: $(wc -l < H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}.txt)"
    echo "Bottom  : $(wc -l < H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_expression.txt)"
    echo "Bottom tau : $(wc -l < low_final.txt)"
    echo "Bottom conservation : $(wc -l < H3K27ac_output_E5_gene_count_down005_${LOW_THRESHOLD}_conservation.txt)"
} > match_summary.txt

############################################
############################################
rm -f high_1.txt high_3.txt medium_1.txt medium_3.txt low_1.txt low_3.txt
rm -f high_final.txt medium_final.txt low_final.txt
rm -rf "$TMPDIR"

echo " !"
echo "results : match_summary.txt"
