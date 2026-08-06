#!/usr/bin/env bash
# Figure 8a, step 02: convert projected intervals from hg38 to hg19 for LDSC.
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# =========================================================
# =========================================================
INDIR="/vol2/zhangshiwen/sheep_cor/liftover_to_hg38/05_hg38_final"

OUTROOT="/vol2/zhangshiwen/sheep_cor/liftover_hg38_to_hg19"

STEP1_DIR="${OUTROOT}/01_hg38_for_chain"
STEP2_DIR="${OUTROOT}/02_hg19_raw"
STEP3_DIR="${OUTROOT}/03_hg19_final"
UNMAP_DIR="${OUTROOT}/unmapped_hg38_to_hg19"

mkdir -p "${STEP1_DIR}" "${STEP2_DIR}" "${STEP3_DIR}" "${UNMAP_DIR}"

SUMMARY="${OUTROOT}/liftover_hg38_to_hg19.summary.tsv"

# =========================================================
# 1. chain file
# =========================================================
CHAIN_HG38_TO_HG19="/vol2/mengzhu/genome/hg38ToHg19.over.chain.gz"

LIFTOVER_MINMATCH="0.8"

# =========================================================
# =========================================================
for cmd in liftOver awk sort zcat wc head basename sed; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] command not found: ${cmd}" >&2
        exit 1
    }
done

[[ -s "${CHAIN_HG38_TO_HG19}" ]] || {
    echo "[ERROR] chain file not found or empty: ${CHAIN_HG38_TO_HG19}" >&2
    exit 1
}

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

sort_bed_file() {
    local infile="$1"
    local outfile="$2"
    LC_ALL=C sort -k1,1V -k2,2n -k3,3n "$infile" > "$outfile"
}

normalize_to_chain_source() {
    local infile="$1"
    local outfile="$2"
    local in_style="$3"
    local chain_source_style="$4"

    case "${in_style},${chain_source_style}" in
        chr,chr)
            cp "${infile}" "${outfile}"
            ;;
        num,num)
            cp "${infile}" "${outfile}"
            ;;
        chr,num)
            awk 'BEGIN{FS=OFS="\t"}{
                sub(/^chr/, "", $1)
                print
            }' "${infile}" > "${outfile}"
            ;;
        num,chr)
            awk 'BEGIN{FS=OFS="\t"}{
                $1 = "chr"$1
                print
            }' "${infile}" > "${outfile}"
            ;;
        *)
            die "invalid input/source style conversion: input=${in_style}, chain_source=${chain_source_style}"
            ;;
    esac
}

finalize_hg19_names() {
    local infile="$1"
    local outfile="$2"
    local chain_target_style="$3"

    case "${chain_target_style}" in
        chr)
            sort_bed_file "${infile}" "${outfile}"
            ;;
        num)
            awk 'BEGIN{FS=OFS="\t"}{
                $1 = "chr"$1
                print
            }' "${infile}" \
            | LC_ALL=C sort -k1,1V -k2,2n -k3,3n \
            > "${outfile}"
            ;;
        *)
            die "invalid human target style: ${chain_target_style}"
            ;;
    esac
}

# =========================================================
# UCSC chain:
# =========================================================
read -r HG19_TARGET_NAME HG38_SOURCE_NAME < <(
    zcat "${CHAIN_HG38_TO_HG19}" | awk '/^chain/ {print $3, $8; exit}'
)

STYLE_HG38_SOURCE=$(detect_style "${HG38_SOURCE_NAME}")
STYLE_HG19_TARGET=$(detect_style "${HG19_TARGET_NAME}")

log "CHAIN_HG38_TO_HG19 source : ${HG38_SOURCE_NAME}"
log "CHAIN_HG38_TO_HG19 target : ${HG19_TARGET_NAME}"
log "style chain source = ${STYLE_HG38_SOURCE}"
log "style chain target = ${STYLE_HG19_TARGET}"

# =========================================================
# =========================================================
echo -e "file\thg38_input_lines\thg38_for_chain_lines\thg19_raw_lines\thg19_final_lines\tunmapped_hg38_to_hg19\thg38_to_hg19_rate\tfinal_to_input_rate" > "${SUMMARY}"

# =========================================================
# =========================================================
shopt -s nullglob
files=( "${INDIR}"/*.bed )

if [[ ${#files[@]} -eq 0 ]]; then
    die "no files matched: ${INDIR}/*.bed"
fi

for bed in "${files[@]}"; do
    base=$(basename "${bed}" .bed)
    base_core="${base%.hg38}"

    log "processing: ${base}"

    [[ -s "${bed}" ]] || {
        log "skip empty file: ${bed}"
        continue
    }

    NCOLS=$(awk 'NR==1{print NF; exit}' "${bed}")
    [[ -n "${NCOLS}" ]] || die "cannot detect column count: ${bed}"

    in_style=$(detect_style "$(awk 'NR==1{print $1; exit}' "${bed}")")
    log "${base}: input style = ${in_style}, columns = ${NCOLS}"

    STEP1_OUT="${STEP1_DIR}/${base}.for_chain.bed"
    STEP2_OUT="${STEP2_DIR}/${base_core}.hg19.raw.bed"
    STEP3_OUT="${STEP3_DIR}/${base_core}.hg19.bed"
    UNMAP="${UNMAP_DIR}/${base_core}.hg38_to_hg19.unmapped.bed"

    # -----------------------------------------------------
    # -----------------------------------------------------
    if done_file "${STEP1_OUT}"; then
        log "${base}: skip step1"
    else
        log "${base}: step1 normalize input to chain source"
        normalize_to_chain_source "${bed}" "${STEP1_OUT}" "${in_style}" "${STYLE_HG38_SOURCE}"
    fi

    # -----------------------------------------------------
    # Step 2: hg38 -> hg19
    # -----------------------------------------------------
    if done_file "${STEP2_OUT}" && exists_file "${UNMAP}"; then
        log "${base}: skip step2"
    else
        log "${base}: step2 liftOver hg38 -> hg19"
        liftOver -minMatch="${LIFTOVER_MINMATCH}" -bedPlus="${NCOLS}" \
            "${STEP1_OUT}" \
            "${CHAIN_HG38_TO_HG19}" \
            "${STEP2_OUT}" \
            "${UNMAP}"
    fi

    # -----------------------------------------------------
    # -----------------------------------------------------
    if done_file "${STEP3_OUT}"; then
        log "${base}: skip step3"
    else
        log "${base}: step3 finalize hg19 names"
        finalize_hg19_names "${STEP2_OUT}" "${STEP3_OUT}" "${STYLE_HG19_TARGET}"
    fi

    in_n=$(line_count "${bed}")
    step1_n=$(line_count "${STEP1_OUT}")
    step2_n=$(line_count "${STEP2_OUT}")
    step3_n=$(line_count "${STEP3_OUT}")
    unmap_n=$(line_count "${UNMAP}")

    rate_hg38_to_hg19=$(calc_ratio "${step2_n}" "${step1_n}")
    rate_final_to_input=$(calc_ratio "${step3_n}" "${in_n}")

    echo -e "${base_core}\t${in_n}\t${step1_n}\t${step2_n}\t${step3_n}\t${unmap_n}\t${rate_hg38_to_hg19}\t${rate_final_to_input}" >> "${SUMMARY}"
done

log "all done"
log "final hg19 files: ${STEP3_DIR}"
log "summary: ${SUMMARY}"
