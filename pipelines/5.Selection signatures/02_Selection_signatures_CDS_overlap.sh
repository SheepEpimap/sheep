#!/bin/bash
set -euo pipefail

# =========================
# 1. Input file
# =========================
fst_file="Demestication_domestic/Demestication_domestic_sheep.bed"
cds_file="/vol2/mengzhu/genome/part_change_esemb100/CDS_esemble100_colin.bed"

# Output directory
outdir="fst_cds_stat"
mkdir -p "${outdir}"

# =========================
# 2. Convert the FST window file to BED
#    Original file:
#    CHROM BIN_START BIN_END ...
#    BED requires:
#    chrom start end
#
#    Assumptions:
#    BIN_START/BIN_END are commonly used interval boundaries,
#    convert to BED with start = BIN_START - 1 and end = BIN_END
# =========================
#awk 'BEGIN{FS=OFS="\t"}
#NR==1 {next}
#{
#    start = $2 - 1
#    if (start < 0) start = 0
    # Retain the original information for downstream tracking
#    print $1, start, $3, $1 ":" $2 "-" $3, $2, $3, $4, $5, $6, $7, $8
#}' "${fst_file}" \
#| sort -k1,1 -k2,2n -k3,3n \
#> "${outdir}/Demestication_domestic_windows.sorted.bed"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {next}
NF>=3 && $1!="" && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ {
    start = $2 - 1
    if (start < 0) start = 0
    print $1, start, $3
}' "${fst_file}" \
| sed 's/\r$//' \
| sort -k1,1 -k2,2n -k3,3n \
> "${outdir}/Demestication_domestic_windows.sorted.bed"

# Field descriptions:
# col1  chrom
# col2  bed_start
# col3  bed_end
# col4  window_id
# col5  original_BIN_START
# col6  original_BIN_END
# col7  N_VARIANTS
# col8  WEIGHTED_FST
# col9  MEAN_FST
# col10 type
# col11 ann

# =========================
# 3. Preprocess the CDS file
#    Keep only the first three columns: chr/start/end
#    Then sort and merge to prevent duplicated CDS segments from affecting the statistics
# =========================
cut -f1-3 "${cds_file}" \
| awk 'BEGIN{FS=OFS="\t"} $2 < $3 {print $1, $2, $3}' \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> "${outdir}/cds.merged.bed"

# =========================
# 4. Count windows that overlap CDS regions
#    -u: Output a window from A once if it has any overlap with B
# =========================
bedtools intersect \
    -a "${outdir}/Demestication_domestic_windows.sorted.bed" \
    -b "${outdir}/cds.merged.bed" \
    -sorted \
    -u \
> "${outdir}/Demestication_domestic_windows.in_CDS.bed"

# =========================
# 5. Count windows with no CDS overlap
#    -v: Output only windows from A that have no overlap at all
# =========================
bedtools intersect \
    -a "${outdir}/Demestication_domestic_windows.sorted.bed" \
    -b "${outdir}/cds.merged.bed" \
    -sorted \
    -v \
> "${outdir}/Demestication_domestic_windows.in_nonCDS.bed"

# =========================
# 6. Label each window as CDS or nonCDS
#    -c: Count the number of CDS overlaps for each window
# =========================
bedtools intersect \
    -a "${outdir}/Demestication_domestic_windows.sorted.bed" \
    -b "${outdir}/cds.merged.bed" \
    -sorted \
    -c \
| awk 'BEGIN{FS=OFS="\t"; print "chrom","bed_start","bed_end","window_id","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST","type","ann","cds_overlap_count","class"}
{
    cls = ($NF > 0 ? "CDS" : "nonCDS")
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,cls
}' \
> "${outdir}/Demestication_domestic_windows.CDS_classification.tsv"

# =========================
# 7. Calculate totals and proportions
# =========================
total_windows=$(wc -l < "${outdir}/Demestication_domestic_windows.sorted.bed")
cds_windows=$(wc -l < "${outdir}/Demestication_domestic_windows.in_CDS.bed")
noncds_windows=$(wc -l < "${outdir}/Demestication_domestic_windows.in_nonCDS.bed")

cds_ratio=$(awk -v a="${cds_windows}" -v b="${total_windows}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')
noncds_ratio=$(awk -v a="${noncds_windows}" -v b="${total_windows}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')

# Percentage format
cds_pct=$(awk -v x="${cds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')
noncds_pct=$(awk -v x="${noncds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')

# Write the summary table
{
    echo -e "category\tcount\tratio\tpercent"
    echo -e "CDS\t${cds_windows}\t${cds_ratio}\t${cds_pct}"
    echo -e "nonCDS\t${noncds_windows}\t${noncds_ratio}\t${noncds_pct}"
    echo -e "total\t${total_windows}\t1.000000\t100.00%"
} > "${outdir}/Demestication_domestic_summary_CDS_vs_nonCDS.tsv"

# Also print to the screen
echo "===== Summary ====="
cat "${outdir}/Demestication_domestic_summary_CDS_vs_nonCDS.tsv"

echo
echo "Output files:"
echo "${outdir}/Demestication_domestic_windows.sorted.bed"
echo "${outdir}/cds.merged.bed"
echo "${outdir}/Demestication_domestic_windows.in_CDS.bed"
echo "${outdir}/Demestication_domestic_windows.in_nonCDS.bed"
echo "${outdir}/Demestication_domestic_windows.CDS_classification.tsv"
echo "${outdir}/Demestication_domestic_summary_CDS_vs_nonCDS.tsv"

cut -f1-3 "Demestication_domestic_windows.sorted.bed" \
| sed 's/\r$//' \
| awk 'BEGIN{FS=OFS="\t"} NF>=3 && $2 < $3 {print $1,$2,$3}' \
| sort -k1,1 -k2,2n -k3,3n \
> "Demestication_domestic_windows_3col.sorted.bed"

bedtools coverage \
    -a "Demestication_domestic_windows_3col.sorted.bed" \
    -b "cds.merged.bed" \
> "Demestication_domestic.coverage_by_CDS.tsv"

read total_bp cds_bp < <(
    awk 'BEGIN{total=0; cds=0}
    {
        # coverage output:
        # col1 chr
        # col2 start
        # col3 end
        # col4 number of overlapping B features
        # col5 number of covered bp in A
        # col6 length of the A interval
        # col7 coverage fraction
        cds += $5
        total += $6
    }
    END{
        print total, cds
    }' "Demestication_domestic.coverage_by_CDS.tsv"
)

noncds_bp=$(( total_bp - cds_bp ))

cds_ratio=$(awk -v a="${cds_bp}" -v b="${total_bp}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')
noncds_ratio=$(awk -v a="${noncds_bp}" -v b="${total_bp}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')

cds_pct=$(awk -v x="${cds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')
noncds_pct=$(awk -v x="${noncds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')

# =========================
# 7. Write results
# =========================
{
    echo -e "category\tbp\tratio\tpercent"
    echo -e "CDS\t${cds_bp}\t${cds_ratio}\t${cds_pct}"
    echo -e "nonCDS\t${noncds_bp}\t${noncds_ratio}\t${noncds_pct}"
    echo -e "total_selected_region\t${total_bp}\t1.000000\t100.00%"
} > "Demestication_domestic_selected_region_CDS_percent.tsv"


for f in *_selected_region_CDS_percent.tsv; do
    sample=${f%_selected_region_CDS_percent.tsv}
    cds_pct=$(awk '$1=="CDS" {print $4}' "$f")
    noncds_pct=$(awk '$1=="nonCDS" {print $4}' "$f")
    echo -e "${sample}\t${cds_pct}\t${noncds_pct}"
done > all_selected_region_CDS_percent_summary.tsv

