#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] line $LINENO, exit code $?" >&2' ERR
shopt -s nullglob

# ========================================================================
# Annotate hg19 CRE files with TSR label via hg38 key in col4.
# Output: col1-3 hg19, col4 hg38_key, col5 CRE_class, col6 TSR_label, col7 TSR_bp
# Tie-break priority: all_common > broad > intermediate > ts_*
# ========================================================================

NCPU="${SLURM_CPUS_PER_TASK:-8}"

: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR in the paths configuration}"
CRE_HG19_DIR="${WORKDIR}/E5_hg19_all"
TSR_DIR="${WORKDIR}/human_TSR_E5_hg38"
LABELED="${TSR_DIR}/labeled.bed"

OUT_DIR="${WORKDIR}/CRE_TSR_E5_hg19"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

[[ -s "${LABELED}" ]] || { echo "ERROR: ${LABELED} not found or empty" >&2; exit 1; }
export LABELED TMP_DIR OUT_DIR

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

annotate_file() {
    local f="$1"
    local base; base=$(basename "$f" .bed)   # e.g. Adipose_adipose_E5.sfCRE.hg19
    local cre_class
    case "${base}" in
        *.sfCRE.hg19) cre_class="sfCRE" ;;
        *.sdCRE.hg19) cre_class="sdCRE" ;;
        *.soCRE.hg19) cre_class="soCRE" ;;
        *.ssCRE.hg19) cre_class="ssCRE" ;;
        *)            cre_class="unknown" ;;
    esac
    local out="${OUT_DIR}/${base}.TSR.bed"

    # 1) decode col4 (hg38_key) → hg38 bed with CRE ID = col4
    local hg38="${TMP_DIR}/${base}.hg38.bed"
    awk 'BEGIN{FS=OFS="\t"}
         NF>=4{ split($4,a,":"); if(length(a)<2) next
                split(a[2],b,"-"); if(length(b)<2) next
                print a[1],b[1],b[2],$4 }' "$f" \
      | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u > "${hg38}"

    if [[ ! -s "${hg38}" ]]; then
        : > "${out}"
        echo -e "${base}\t0\t0\tempty_input"
        return 0
    fi

    # 2) intersect with TSR labeled → (cre_id, tsr_label, overlap_bp) rows
    #    with -wao: every -a row gets at least one -b row (label '.' if no overlap)
    local isx="${TMP_DIR}/${base}.isx.tsv"
    bedtools intersect -a "${hg38}" -b "${LABELED}" -wao \
      | awk 'BEGIN{FS=OFS="\t"} {print $4,$8,$9}' > "${isx}"

    # 3) dominant label per CRE (bp-weighted, tie-break by priority)
    local dom="${TMP_DIR}/${base}.dom.tsv"
    awk 'BEGIN{FS=OFS="\t"}
         function prio(l){
             if(l=="all_common")   return 4
             if(l=="broad")        return 3
             if(l=="intermediate") return 2
             if(l ~ /^ts_/)        return 1
             return 0
         }
         {
             cid=$1; lab=$2; bp=$3+0
             seen[cid]=1
             if(lab=="." || lab=="") next
             key=cid SUBSEP lab
             acc[key]+=bp
         }
         END{
             for(k in acc){
                 split(k,a,SUBSEP); c=a[1]; l=a[2]; b=acc[k]; p=prio(l)
                 if(!(c in bb) || b>bb[c] || (b==bb[c] && p>bp_p[c])){
                     bb[c]=b; bl[c]=l; bp_p[c]=p
                 }
             }
             for(c in seen){
                 if(c in bl) print c,bl[c],bb[c]
                 else        print c,"none",0
             }
         }' "${isx}" > "${dom}"

    # 4) join back to hg19 by col4 (CRE ID = hg38 key)
    LC_ALL=C sort -k4,4 "$f"     > "${TMP_DIR}/${base}.hg19.s"
    LC_ALL=C sort -k1,1 "${dom}" > "${TMP_DIR}/${base}.dom.s"
    join -t $'\t' -1 4 -2 1 -a 1 -e "none" \
         -o '1.1 1.2 1.3 1.4 2.2 2.3' \
         "${TMP_DIR}/${base}.hg19.s" "${TMP_DIR}/${base}.dom.s" \
      | awk -v C="${cre_class}" 'BEGIN{FS=OFS="\t"}
                                  {if($6=="none")$6=0; print $1,$2,$3,$4,C,$5,$6}' \
      | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "${out}"

    rm -f "${hg38}" "${isx}" "${dom}" \
          "${TMP_DIR}/${base}.hg19.s" "${TMP_DIR}/${base}.dom.s"

    local n_in n_out n_none
    n_in=$(wc -l < "$f")
    n_out=$(wc -l < "${out}")
    n_none=$(awk -F'\t' '$6=="none"' "${out}" | wc -l)
    echo -e "${base}\t${cre_class}\t${n_in}\t${n_out}\t${n_none}"
}
export -f annotate_file

log "===== Annotate CRE files with TSR ====="
TASKS="${TMP_DIR}/tasks.txt"; : > "${TASKS}"
for f in "${CRE_HG19_DIR}"/*.sfCRE.hg19.bed \
         "${CRE_HG19_DIR}"/*.sdCRE.hg19.bed \
         "${CRE_HG19_DIR}"/*.soCRE.hg19.bed \
         "${CRE_HG19_DIR}"/*.ssCRE.hg19.bed; do
    [[ -e "$f" ]] && echo "$f" >> "${TASKS}"
done
log "  total tasks: $(wc -l < "${TASKS}")"

SUMMARY="${OUT_DIR}/annotation_summary.tsv"
{
echo -e "base\tcre_class\tn_input\tn_output\tn_none"
xargs -P "${NCPU}" -I {} bash -c 'annotate_file "$1"' _ {} < "${TASKS}"
} > "${SUMMARY}"

log "Summary: ${SUMMARY}"
log "===== ALL DONE ====="
