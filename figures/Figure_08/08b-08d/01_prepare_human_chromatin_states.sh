#!/usr/bin/env bash
# Human hg38 E1-E9 chromatin states to sheep ARS-UI_Ramb_v3.0
# (GCF_016772045.2), followed by per-tissue CRE classification in sheep V3
# coordinates. This is the workflow used for Figure 8b-d.
set -euo pipefail
trap 'echo "[ERROR] line ${LINENO}, exit code $?" >&2' ERR

NCPU="${SLURM_CPUS_PER_TASK:-12}"
LIFTOVER_BIN="${LIFTOVER_BIN:-liftOver}"
HUMAN_TO_SHEEP_CHAIN="${HUMAN_TO_SHEEP_CHAIN:?Set HUMAN_TO_SHEEP_CHAIN to hg38ToGCF_016772045.2.over.chain.gz}"
SHEEP_V3_CHROM_MAP="${SHEEP_V3_CHROM_MAP:?Set SHEEP_V3_CHROM_MAP to the V3 numeric-to-NC chromosome map}"
HUMAN_HG38_DIR="${HUMAN_HG38_DIR:?Set HUMAN_HG38_DIR to the human hg38 E1-E9 BED directory}"
SHEEP_V3_DIR="${SHEEP_V3_DIR:?Set SHEEP_V3_DIR to the sheep V3 E1-E9 BED directory}"
WORKDIR="${WORKDIR:?Set WORKDIR to an output directory}"

MINMATCH_HG38_V3="${MINMATCH_HG38_V3:-0.1}"
OVERLAP_FRAC="${OVERLAP_FRAC:-0.5}"
LIFTED_DIR="${WORKDIR}/lifted"
CLASS_DIR="${WORKDIR}/classified"
SHEEP_REF_DIR="${WORKDIR}/sheep_v3_ref"
TMP_DIR="${WORKDIR}/tmp"

STATES=(E1 E2 E3 E4 E5 E6 E7 E8 E9)
HUMAN_TISSUES=(Adipose Colon Cortex Heart Liver Lung Muscle Ovary Sintest Spleen Stomach Testis)
SHEEP_TISSUES=(adipose colon cerebral-cortex heart liver lung muscle ovary jejunum spleen abomasum testis)

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

for cmd in "${LIFTOVER_BIN}" awk bedtools sort zcat sed xargs; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] command not found: ${cmd}" >&2
        exit 1
    }
done
[[ -f "${HUMAN_TO_SHEEP_CHAIN}" ]] || {
    echo "[ERROR] chain not found: ${HUMAN_TO_SHEEP_CHAIN}" >&2
    exit 1
}
[[ -f "${SHEEP_V3_CHROM_MAP}" ]] || {
    echo "[ERROR] chromosome map not found: ${SHEEP_V3_CHROM_MAP}" >&2
    exit 1
}
mkdir -p "${LIFTED_DIR}" "${CLASS_DIR}" "${SHEEP_REF_DIR}" "${TMP_DIR}"

# The query assembly in the chain is sheep V3. Detect its chromosome naming
# convention so all lifted intervals can be normalized to numeric chromosomes.
V3_QNAME=$(set +o pipefail; zcat "${HUMAN_TO_SHEEP_CHAIN}" | awk '/^chain/{print $8; exit}')
case "${V3_QNAME}" in
    NC_*) V3_CHAIN_STYLE=nc ;;
    chr*) V3_CHAIN_STYLE=chr ;;
    [0-9]*) V3_CHAIN_STYLE=num ;;
    *) V3_CHAIN_STYLE=other ;;
esac
log "Chain query assembly: ${V3_QNAME} (style=${V3_CHAIN_STYLE})"

NC2NUM="${TMP_DIR}/v3nc2num.tsv"
awk 'BEGIN{OFS="\t"}{print $2, $1}' "${SHEEP_V3_CHROM_MAP}" > "${NC2NUM}"

# Prepare the per-tissue sheep V3 reference states in numeric chromosome style.
for sheep_tissue in "${SHEEP_TISSUES[@]}"; do
    for state in "${STATES[@]}"; do
        src="${SHEEP_V3_DIR}/${sheep_tissue}_${state}.v3forHg38.bed"
        dst="${SHEEP_REF_DIR}/${sheep_tissue}_${state}.v3.bed"
        [[ -s "${dst}" ]] && continue
        if [[ ! -s "${src}" ]]; then
            : > "${dst}"
            continue
        fi
        sed 's/^chr//' "${src}" | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "${dst}"
    done
done

process_lift() {
    local human_tissue="$1" state="$2"
    local out_v3="${LIFTED_DIR}/${human_tissue}_${state}.v3.bed"
    local out_unmapped="${LIFTED_DIR}/${human_tissue}_${state}.unmapped"
    [[ -s "${out_v3}" ]] && return 0

    local hg38_file=""
    local candidate
    for candidate in \
        "${HUMAN_HG38_DIR}/${human_tissue}_${state}.bed" \
        "${HUMAN_HG38_DIR}/${human_tissue}_${state}.hg38.bed"; do
        [[ -s "${candidate}" ]] && { hg38_file="${candidate}"; break; }
    done
    if [[ -z "${hg38_file}" ]]; then
        : > "${out_v3}"
        : > "${out_unmapped}"
        echo "[WARN] ${human_tissue}_${state}: no hg38 input found" >&2
        return 0
    fi

    local tmp="${TMP_DIR}/${human_tissue}_${state}"
    awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' "${hg38_file}" > "${tmp}.encoded.bed"
    local columns
    columns=$(awk 'NR==1{print NF;exit}' "${tmp}.encoded.bed")
    "${LIFTOVER_BIN}" -minMatch="${MINMATCH_HG38_V3}" -bedPlus="${columns}" \
        "${tmp}.encoded.bed" "${HUMAN_TO_SHEEP_CHAIN}" \
        "${tmp}.v3raw.bed" "${out_unmapped}" || true

    if [[ -s "${tmp}.v3raw.bed" ]]; then
        case "${V3_CHAIN_STYLE}" in
            nc)
                awk 'BEGIN{FS=OFS="\t"} NR==FNR{m[$1]=$2;next} ($1 in m){$1=m[$1];print}' \
                    "${NC2NUM}" "${tmp}.v3raw.bed" \
                    | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "${out_v3}"
                ;;
            chr)
                sed 's/^chr//' "${tmp}.v3raw.bed" \
                    | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "${out_v3}"
                ;;
            *)
                LC_ALL=C sort -k1,1V -k2,2n -k3,3n "${tmp}.v3raw.bed" > "${out_v3}"
                ;;
        esac
    else
        : > "${out_v3}"
    fi
    rm -f "${tmp}.encoded.bed" "${tmp}.v3raw.bed"
    echo "[DONE] ${human_tissue}_${state}: $(wc -l < "${out_v3}") lifted intervals"
}
export -f process_lift
export LIFTOVER_BIN HUMAN_TO_SHEEP_CHAIN MINMATCH_HG38_V3 HUMAN_HG38_DIR
export LIFTED_DIR TMP_DIR V3_CHAIN_STYLE NC2NUM

LIFT_TASKS="${TMP_DIR}/lift_tasks.txt"
: > "${LIFT_TASKS}"
for human_tissue in "${HUMAN_TISSUES[@]}"; do
    for state in "${STATES[@]}"; do
        echo "${human_tissue} ${state}" >> "${LIFT_TASKS}"
    done
done
xargs -P "${NCPU}" -L 1 bash -c 'process_lift "$@"' _ < "${LIFT_TASKS}"
rm -f "${LIFT_TASKS}"

is_promoter() {
    case "$1" in E1|E2|E3|E4) return 0 ;; *) return 1 ;; esac
}
to_midpoint() { awk -v OFS="\t" '{mid=int(($2+$3)/2); print $1,mid,mid+1,$4}' "$1"; }
recover_by_name() { awk 'NR==FNR{ids[$4];next} $4 in ids' "$1" "$2" | sort -u; }

process_classify() {
    local human_tissue="$1" sheep_tissue="$2" state="$3"
    local human_v3="${LIFTED_DIR}/${human_tissue}_${state}.v3.bed"
    [[ -s "${human_v3}" ]] || return 0
    local prefix="${CLASS_DIR}/${human_tissue}_${sheep_tissue}_${state}"
    local sheep_same sheep_other
    sheep_same=$(mktemp)
    sheep_other=$(mktemp)

    local ref="${SHEEP_REF_DIR}/${sheep_tissue}_${state}.v3.bed"
    if [[ -s "${ref}" ]]; then
        bedtools sort -i "${ref}" | bedtools merge -i - > "${sheep_same}"
    else
        : > "${sheep_same}"
    fi

    local other_files=() other_state other_file
    for other_state in E1 E2 E3 E4 E5 E6 E7 E8 E9; do
        [[ "${other_state}" == "${state}" ]] && continue
        other_file="${SHEEP_REF_DIR}/${sheep_tissue}_${other_state}.v3.bed"
        [[ -s "${other_file}" ]] && other_files+=("${other_file}")
    done
    if [[ ${#other_files[@]} -gt 0 ]]; then
        cat "${other_files[@]}" | bedtools sort -i - | bedtools merge -i - > "${sheep_other}"
    else
        : > "${sheep_other}"
    fi

    local sf="${prefix}.sfCRE.bed" sd="${prefix}.sdCRE.bed" so="${prefix}.soCRE.bed"
    local non_sf
    non_sf=$(mktemp)
    if is_promoter "${state}"; then
        local mid_all mid_sf mid_non_sf mid_sd
        mid_all=$(mktemp); mid_sf=$(mktemp); mid_non_sf=$(mktemp); mid_sd=$(mktemp)
        to_midpoint "${human_v3}" > "${mid_all}"
        bedtools intersect -a "${mid_all}" -b "${sheep_same}" -f 1.0 -wa | sort -u > "${mid_sf}"
        recover_by_name "${mid_sf}" "${human_v3}" > "${sf}"
        bedtools intersect -a "${human_v3}" -b "${sf}" -v > "${non_sf}"
        to_midpoint "${non_sf}" > "${mid_non_sf}"
        bedtools intersect -a "${mid_non_sf}" -b "${sheep_other}" -f 1.0 -wa | sort -u > "${mid_sd}"
        recover_by_name "${mid_sd}" "${non_sf}" > "${sd}"
        bedtools intersect -a "${non_sf}" -b "${sd}" -v > "${so}"
        rm -f "${mid_all}" "${mid_sf}" "${mid_non_sf}" "${mid_sd}"
    else
        bedtools intersect -a "${human_v3}" -b "${sheep_same}" -f "${OVERLAP_FRAC}" -wa | sort -u > "${sf}"
        bedtools intersect -a "${human_v3}" -b "${sheep_same}" -f "${OVERLAP_FRAC}" -v > "${non_sf}"
        bedtools intersect -a "${non_sf}" -b "${sheep_other}" -f "${OVERLAP_FRAC}" -wa | sort -u > "${sd}"
        bedtools intersect -a "${non_sf}" -b "${sheep_other}" -f "${OVERLAP_FRAC}" -v > "${so}"
    fi
    rm -f "${non_sf}" "${sheep_same}" "${sheep_other}"
    echo "[DONE] ${human_tissue}/${sheep_tissue} ${state}: sf=$(wc -l < "${sf}") sd=$(wc -l < "${sd}") so=$(wc -l < "${so}")"
}
export -f is_promoter to_midpoint recover_by_name process_classify
export LIFTED_DIR CLASS_DIR SHEEP_REF_DIR OVERLAP_FRAC

CLASS_TASKS="${TMP_DIR}/class_tasks.txt"
: > "${CLASS_TASKS}"
for index in "${!HUMAN_TISSUES[@]}"; do
    for state in "${STATES[@]}"; do
        echo "${HUMAN_TISSUES[$index]} ${SHEEP_TISSUES[$index]} ${state}" >> "${CLASS_TASKS}"
    done
done
xargs -P "${NCPU}" -L 1 bash -c 'process_classify "$@"' _ < "${CLASS_TASKS}"
rm -f "${CLASS_TASKS}"

SUMMARY="${WORKDIR}/per_tissue_classification.tsv"
{
    printf 'human_tissue\tsheep_tissue\tstate\tsfCRE\tsdCRE\tsoCRE\tssCRE\ttotal\n'
    for index in "${!HUMAN_TISSUES[@]}"; do
        human_tissue="${HUMAN_TISSUES[$index]}"
        sheep_tissue="${SHEEP_TISSUES[$index]}"
        for state in "${STATES[@]}"; do
            prefix="${CLASS_DIR}/${human_tissue}_${sheep_tissue}_${state}"
            sf=0; sd=0; so=0; ss=0
            [[ -s "${prefix}.sfCRE.bed" ]] && sf=$(wc -l < "${prefix}.sfCRE.bed")
            [[ -s "${prefix}.sdCRE.bed" ]] && sd=$(wc -l < "${prefix}.sdCRE.bed")
            [[ -s "${prefix}.soCRE.bed" ]] && so=$(wc -l < "${prefix}.soCRE.bed")
            unmapped="${LIFTED_DIR}/${human_tissue}_${state}.unmapped"
            [[ -s "${unmapped}" ]] && ss=$(awk '!/^#/{n++} END{print n+0}' "${unmapped}")
            total=$((sf + sd + so + ss))
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${human_tissue}" "${sheep_tissue}" "${state}" \
                "${sf}" "${sd}" "${so}" "${ss}" "${total}"
        done
    done
} > "${SUMMARY}"
log "Completed. Summary: ${SUMMARY}"
