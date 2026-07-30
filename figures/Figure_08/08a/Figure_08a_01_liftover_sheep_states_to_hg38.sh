#!/usr/bin/env bash
# Figure 8a, step 01: project sheep chromatin-state intervals to hg38.
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# =========================================================
# =========================================================
INDIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/All_chromatin_state"
OUTROOT="/vol2/zhangshiwen/sheep_cor/liftover_to_hg38"

TMPDIR="${OUTROOT}/tmp"
STEP1_DIR="${OUTROOT}/01_v2_sourceNC"
STEP2_DIR="${OUTROOT}/02_v3_raw"
STEP3_DIR="${OUTROOT}/03_v3_for_hg38"
STEP4_DIR="${OUTROOT}/04_hg38_raw"
STEP5_DIR="${OUTROOT}/05_hg38_final"
UNMAP1_DIR="${OUTROOT}/unmapped_v2_to_v3"
UNMAP2_DIR="${OUTROOT}/unmapped_v3_to_hg38"

mkdir -p "${TMPDIR}" "${STEP1_DIR}" "${STEP2_DIR}" "${STEP3_DIR}" \
         "${STEP4_DIR}" "${STEP5_DIR}" "${UNMAP1_DIR}" "${UNMAP2_DIR}"

SUMMARY="${OUTROOT}/liftover_summary.tsv"

# =========================================================
# =========================================================
CHAIN_V2_TO_V3="/vol2/mengzhu/genome/GCF_016772045.1ToGCF_016772045.2.over.chain.gz"
CHAIN_V3_TO_HG38="/vol2/mengzhu/genome/GCF_016772045.2ToHg38_chr.over.chain.gz"

MAP_V2_NC_TO_CHR="/data/home/sczd644/run/zsw_chrombpnet/phylop/v2NCtochr.txt"   # NC_056054.1 -> chr1
MAP_V3_NUM_TO_NC="/data/home/sczd644/run/zsw_chrombpnet/phylop/v3NCtochr.txt"   # 1 -> NC_056054.1

LIFTOVER_MINMATCH="0.8"

# =========================================================
# =========================================================
for cmd in liftOver awk sort zcat wc head basename; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] command not found: ${cmd}" >&2
        exit 1
    }
done

# =========================================================
# =========================================================
log() {
    echo "[INFO] $*" >&2
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

done_file() {
    [[ -s "$1" ]]
}

exists_file() {
    [[ -e "$1" ]]
}

line_count() {
    if [[ -e "$1" ]]; then
        wc -l < "$1"
    else
        echo 0
    fi
}

calc_ratio() {
    local num="$1"
    local den="$2"
    awk -v n="$num" -v d="$den" 'BEGIN{
        if (d == 0) {
            print "NA"
        } else {
            printf "%.6f\n", n/d
        }
    }'
}

detect_style() {
    local name="$1"
    if [[ "$name" =~ ^chr ]]; then
        echo "chr"
    elif [[ "$name" =~ ^NC_ ]]; then
        echo "nc"
    elif [[ "$name" =~ ^[0-9]+$ ]]; then
        echo "num"
    else
        echo "other"
    fi
}

rename_first_col_by_map() {
    local mapfile="$1"
    local infile="$2"
    local outfile="$3"

    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR {
        map[$1]=$2
        next
    }
    {
        if ($1 in map) {
            $1 = map[$1]
            print
        }
    }' "$mapfile" "$infile" > "$outfile"
}

adjust_nc_suffix_first_col() {
    local infile="$1"
    local outfile="$2"
    local suffix="$3"

    awk -v suf="$suffix" 'BEGIN{FS=OFS="\t"}
    {
        sub(/\.[0-9]+$/, suf, $1)
        print
    }' "$infile" > "$outfile"
}

sort_bed_file() {
    local infile="$1"
    local outfile="$2"
    LC_ALL=C sort -k1,1V -k2,2n -k3,3n "$infile" > "$outfile"
}

# =========================================================
# UCSC chain:
# =========================================================
read -r V2V3_TARGET_NAME V2V3_SOURCE_NAME < <(
    zcat "${CHAIN_V2_TO_V3}" | awk '/^chain/ {print $3, $8; exit}'
)

read -r V3HG_TARGET_NAME V3HG_SOURCE_NAME < <(
    zcat "${CHAIN_V3_TO_HG38}" | awk '/^chain/ {print $3, $8; exit}'
)

log "CHAIN_V2_TO_V3 target : ${V2V3_TARGET_NAME}"
log "CHAIN_V2_TO_V3 source : ${V2V3_SOURCE_NAME}"
log "CHAIN_V3_TO_HG38 target : ${V3HG_TARGET_NAME}"
log "CHAIN_V3_TO_HG38 source : ${V3HG_SOURCE_NAME}"

STYLE_V2V3_TARGET=$(detect_style "${V2V3_TARGET_NAME}")
STYLE_V2V3_SOURCE=$(detect_style "${V2V3_SOURCE_NAME}")
STYLE_V3HG_SOURCE=$(detect_style "${V3HG_SOURCE_NAME}")
STYLE_V3HG_TARGET=$(detect_style "${V3HG_TARGET_NAME}")

log "style chain1 source = ${STYLE_V2V3_SOURCE}"
log "style chain1 target = ${STYLE_V2V3_TARGET}"
log "style chain2 source = ${STYLE_V3HG_SOURCE}"
log "style chain2 target = ${STYLE_V3HG_TARGET}"

V2V3_SOURCE_SUFFIX=$(echo "${V2V3_SOURCE_NAME}" | sed -E 's/^.*(\.[0-9]+)$/\1/')
V2V3_TARGET_SUFFIX=$(echo "${V2V3_TARGET_NAME}" | sed -E 's/^.*(\.[0-9]+)$/\1/')

# =========================================================
# =========================================================

CHR_TO_V2SRC_NC="${TMPDIR}/v2_chr_to_chain1_sourceNC.tsv"
if done_file "${CHR_TO_V2SRC_NC}"; then
    log "skip map: ${CHR_TO_V2SRC_NC}"
else
    log "building map: ${CHR_TO_V2SRC_NC}"
    awk -v suf="${V2V3_SOURCE_SUFFIX}" 'BEGIN{OFS="\t"}
    {
        nc=$1
        sub(/\.[0-9]+$/, suf, nc)
        print $2, nc
    }' "${MAP_V2_NC_TO_CHR}" > "${CHR_TO_V2SRC_NC}"
fi

V3_NC_TO_CHR="${TMPDIR}/v3_nc_to_chr.tsv"
if done_file "${V3_NC_TO_CHR}"; then
    log "skip map: ${V3_NC_TO_CHR}"
else
    log "building map: ${V3_NC_TO_CHR}"
    awk 'BEGIN{OFS="\t"} {print $2, "chr"$1}' "${MAP_V3_NUM_TO_NC}" > "${V3_NC_TO_CHR}"
fi

V3_NUM_TO_CHR="${TMPDIR}/v3_num_to_chr.tsv"
if done_file "${V3_NUM_TO_CHR}"; then
    log "skip map: ${V3_NUM_TO_CHR}"
else
    log "building map: ${V3_NUM_TO_CHR}"
    awk 'BEGIN{OFS="\t"} {print $1, "chr"$1}' "${MAP_V3_NUM_TO_NC}" > "${V3_NUM_TO_CHR}"
fi

# =========================================================
# =========================================================
echo -e "file\torig_lines\tv2_sourceNC_lines\tv3_raw_lines\tv3_for_hg38_lines\thg38_raw_lines\thg38_final_lines\tunmapped_v2_to_v3\tunmapped_v3_to_hg38\tv2_to_v3_rate\tv3_to_hg38_rate\tfinal_to_orig_rate" > "${SUMMARY}"

# =========================================================
# =========================================================
shopt -s nullglob
files=( "${INDIR}"/*_E*.bed )

if [[ ${#files[@]} -eq 0 ]]; then
    die "no files matched: ${INDIR}/*_E*.bed"
fi

for bed in "${files[@]}"; do
    base=$(basename "${bed}" .bed)
    log "processing: ${base}"

    [[ -s "${bed}" ]] || {
        log "skip empty file: ${bed}"
        continue
    }

    NCOLS=$(awk 'NR==1{print NF; exit}' "${bed}")
    [[ -n "${NCOLS}" ]] || die "cannot detect column count: ${bed}"

    in_style=$(detect_style "$(awk 'NR==1{print $1; exit}' "${bed}")")
    log "${base}: input style = ${in_style}, columns = ${NCOLS}"

    STEP1_OUT="${STEP1_DIR}/${base}.v2sourceNC.bed"
    STEP2_OUT="${STEP2_DIR}/${base}.v3raw.bed"
    STEP3_OUT="${STEP3_DIR}/${base}.v3forHg38.bed"
    STEP4_OUT="${STEP4_DIR}/${base}.hg38.raw.bed"
    STEP5_OUT="${STEP5_DIR}/${base}.hg38.bed"

    UNMAP1="${UNMAP1_DIR}/${base}.v2_to_v3.unmapped.bed"
    UNMAP2="${UNMAP2_DIR}/${base}.v3_to_hg38.unmapped.bed"

    # -----------------------------------------------------
    # Step 1: input -> chain1 source
    # -----------------------------------------------------
    if done_file "${STEP1_OUT}"; then
        log "${base}: skip step1"
    else
        log "${base}: step1 input -> chain1 source"

        step1_key="${in_style},${STYLE_V2V3_SOURCE}"
        case "${step1_key}" in
            chr,nc)
                rename_first_col_by_map "${CHR_TO_V2SRC_NC}" "${bed}" "${STEP1_OUT}"
                ;;
            nc,nc)
                adjust_nc_suffix_first_col "${bed}" "${STEP1_OUT}" "${V2V3_SOURCE_SUFFIX}"
                ;;
            chr,chr)
                cp "${bed}" "${STEP1_OUT}"
                ;;
            *)
                die "${base}: invalid step1 style conversion: input=${in_style}, chain1_source=${STYLE_V2V3_SOURCE}"
                ;;
        esac
    fi

    # -----------------------------------------------------
    # Step 2: v2 -> v3
    # -----------------------------------------------------
    if done_file "${STEP2_OUT}" && exists_file "${UNMAP1}"; then
        log "${base}: skip step2"
    else
        log "${base}: step2 liftOver v2 -> v3"
        liftOver -minMatch="${LIFTOVER_MINMATCH}" -bedPlus="${NCOLS}" \
            "${STEP1_OUT}" \
            "${CHAIN_V2_TO_V3}" \
            "${STEP2_OUT}" \
            "${UNMAP1}"
    fi

    # -----------------------------------------------------
    # Step 3: chain1 target -> chain2 source
    # -----------------------------------------------------
    if done_file "${STEP3_OUT}"; then
        log "${base}: skip step3"
    else
        log "${base}: step3 bridge chain1 target -> chain2 source"

        step3_key="${STYLE_V2V3_TARGET},${STYLE_V3HG_SOURCE}"
        case "${step3_key}" in
            nc,chr)
                rename_first_col_by_map "${V3_NC_TO_CHR}" "${STEP2_OUT}" "${STEP3_OUT}"
                ;;
            num,chr)
                rename_first_col_by_map "${V3_NUM_TO_CHR}" "${STEP2_OUT}" "${STEP3_OUT}"
                ;;
            chr,chr)
                cp "${STEP2_OUT}" "${STEP3_OUT}"
                ;;
            nc,nc)
                cp "${STEP2_OUT}" "${STEP3_OUT}"
                ;;
            *)
                die "${base}: invalid bridge style: chain1_target=${STYLE_V2V3_TARGET}, chain2_source=${STYLE_V3HG_SOURCE}"
                ;;
        esac
    fi

    # -----------------------------------------------------
    # Step 4: v3 -> hg38
    # -----------------------------------------------------
    if done_file "${STEP4_OUT}" && exists_file "${UNMAP2}"; then
        log "${base}: skip step4"
    else
        log "${base}: step4 liftOver v3 -> hg38"
        liftOver -minMatch="${LIFTOVER_MINMATCH}" -bedPlus="${NCOLS}" \
            "${STEP3_OUT}" \
            "${CHAIN_V3_TO_HG38}" \
            "${STEP4_OUT}" \
            "${UNMAP2}"
    fi

    # -----------------------------------------------------
    # -----------------------------------------------------
    if done_file "${STEP5_OUT}"; then
        log "${base}: skip step5"
    else
        log "${base}: step5 finalize hg38 names"

        case "${STYLE_V3HG_TARGET}" in
            chr)
                sort_bed_file "${STEP4_OUT}" "${STEP5_OUT}"
                ;;
            num)
                awk 'BEGIN{FS=OFS="\t"} {$1="chr"$1; print}' "${STEP4_OUT}" \
                | LC_ALL=C sort -k1,1V -k2,2n -k3,3n \
                > "${STEP5_OUT}"
                ;;
            *)
                die "${base}: invalid human target style: ${STYLE_V3HG_TARGET}"
                ;;
        esac
    fi

    orig_n=$(line_count "${bed}")
    step1_n=$(line_count "${STEP1_OUT}")
    step2_n=$(line_count "${STEP2_OUT}")
    step3_n=$(line_count "${STEP3_OUT}")
    step4_n=$(line_count "${STEP4_OUT}")
    step5_n=$(line_count "${STEP5_OUT}")
    unmap1_n=$(line_count "${UNMAP1}")
    unmap2_n=$(line_count "${UNMAP2}")

    rate_v2_to_v3=$(calc_ratio "${step2_n}" "${orig_n}")
    rate_v3_to_hg38=$(calc_ratio "${step4_n}" "${step3_n}")
    rate_final_to_orig=$(calc_ratio "${step5_n}" "${orig_n}")

    echo -e "${base}\t${orig_n}\t${step1_n}\t${step2_n}\t${step3_n}\t${step4_n}\t${step5_n}\t${unmap1_n}\t${unmap2_n}\t${rate_v2_to_v3}\t${rate_v3_to_hg38}\t${rate_final_to_orig}" >> "${SUMMARY}"
done

log "all done"
log "final hg38 files: ${STEP5_DIR}"
log "summary: ${SUMMARY}"
