#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] line $LINENO, exit code $?" >&2' ERR

# Reuse the final G2_human_to_sheep liftover/reference files and rebuild only
# the E1-E9 classification affected by treating E10/E11 as blank.
#
# Preserved rule:
#   E1-E4: midpoint overlap
#   E5-E9: fractional overlap, -f 0.5
#
# Rebuilt outputs:
#   classified/*.sfCRE.bed / *.sdCRE.bed / *.soCRE.bed
#   per_tissue_classification.tsv
#
# Reused inputs from SRC_WORKDIR:
#   lifted/{human_tissue}_{E1-E9}.v3.bed
#   lifted/{human_tissue}_{E1-E9}.unmapped
#   sheep_v3_ref/{sheep_tissue}_{E1-E9}.v3.bed

NCPU="${SLURM_CPUS_PER_TASK:-12}"

: "${SRC_WORKDIR:?Set SRC_WORKDIR in the paths configuration}"
: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR in the paths configuration}"

USE_SYMLINKS="${USE_SYMLINKS:-1}"
FORCE_RELINK="${FORCE_RELINK:-0}"
OVERLAP_FRAC="${OVERLAP_FRAC:-0.5}"

SRC_LIFTED_DIR="${SRC_WORKDIR}/lifted"
SRC_SHEEP_REF_DIR="${SRC_WORKDIR}/sheep_v3_ref"
SRC_SUMMARY="${SRC_WORKDIR}/per_tissue_classification.tsv"

LIFTED_DIR="${WORKDIR}/lifted"
CLASS_DIR="${WORKDIR}/classified"
SHEEP_REF_DIR="${WORKDIR}/sheep_v3_ref"
TMP_DIR="${WORKDIR}/tmp_reclassify"

STATES=(E1 E2 E3 E4 E5 E6 E7 E8 E9)
HUMAN_TISSUES=(Adipose Colon Cortex Heart Liver Lung Muscle Ovary Sintest Spleen Stomach Testis)
SHEEP_TISSUES=(adipose colon cerebral-cortex heart liver lung muscle ovary jejunum spleen abomasum testis)

log() { echo "[$(date +%F' '%T)] $*" >&2; }

mkdir -p "${LIFTED_DIR}" "${CLASS_DIR}" "${SHEEP_REF_DIR}" "${TMP_DIR}"

MANIFEST="${WORKDIR}/reuse_manifest.tsv"
{
    echo -e "kind\tsource\ttarget\taction"
} > "${MANIFEST}"

reuse_file() {
    local kind="$1" src="$2" dst="$3"

    [[ -e "${src}" ]] || { echo "ERROR: missing source ${kind}: ${src}" >&2; exit 1; }

    if [[ -e "${dst}" || -L "${dst}" ]]; then
        if [[ "${FORCE_RELINK}" == "1" ]]; then
            rm -f "${dst}"
        else
            echo -e "${kind}\t${src}\t${dst}\tkept_existing" >> "${MANIFEST}"
            return 0
        fi
    fi

    if [[ "${USE_SYMLINKS}" == "1" ]]; then
        ln -s "${src}" "${dst}"
        echo -e "${kind}\t${src}\t${dst}\tsymlinked" >> "${MANIFEST}"
    else
        cp -p "${src}" "${dst}"
        echo -e "${kind}\t${src}\t${dst}\tcopied" >> "${MANIFEST}"
    fi
}

log "===== Step 1: reuse existing E1-E9 inputs ====="
log "source: ${SRC_WORKDIR}"
log "target: ${WORKDIR}"

for ht in "${HUMAN_TISSUES[@]}"; do
    for state in "${STATES[@]}"; do
        reuse_file "lifted_v3" \
            "${SRC_LIFTED_DIR}/${ht}_${state}.v3.bed" \
            "${LIFTED_DIR}/${ht}_${state}.v3.bed"
        reuse_file "lifted_unmapped" \
            "${SRC_LIFTED_DIR}/${ht}_${state}.unmapped" \
            "${LIFTED_DIR}/${ht}_${state}.unmapped"
    done
done

for st in "${SHEEP_TISSUES[@]}"; do
    for state in "${STATES[@]}"; do
        reuse_file "sheep_ref_v3" \
            "${SRC_SHEEP_REF_DIR}/${st}_${state}.v3.bed" \
            "${SHEEP_REF_DIR}/${st}_${state}.v3.bed"
    done
done

log "reuse manifest: ${MANIFEST}"

export LIFTED_DIR CLASS_DIR SHEEP_REF_DIR TMP_DIR OVERLAP_FRAC

is_promoter() {
    case "$1" in
        E1|E2|E3|E4) return 0 ;;
        *) return 1 ;;
    esac
}
export -f is_promoter

to_midpoint() {
    awk -v OFS="\t" '{mid=int(($2+$3)/2); print $1,mid,mid+1,$4}' "$1"
}
export -f to_midpoint

recover_by_name() {
    awk 'NR==FNR{ids[$4];next} $4 in ids' "$1" "$2" | sort -u
}
export -f recover_by_name

process_classify() {
    local human_t="$1" sheep_t="$2" state="$3"
    local human_v3="${LIFTED_DIR}/${human_t}_${state}.v3.bed"
    local prefix="${CLASS_DIR}/${human_t}_${sheep_t}_${state}"
    local sf="${prefix}.sfCRE.bed"
    local sd="${prefix}.sdCRE.bed"
    local so="${prefix}.soCRE.bed"

    if [[ ! -s "${human_v3}" ]]; then
        : > "${sf}"
        : > "${sd}"
        : > "${so}"
        echo "[EMPTY] ${human_t}->${sheep_t} ${state}"
        return 0
    fi

    local sheep_same; sheep_same=$(mktemp)
    local ref="${SHEEP_REF_DIR}/${sheep_t}_${state}.v3.bed"
    if [[ -s "${ref}" ]]; then
        bedtools sort -i "${ref}" | bedtools merge -i - > "${sheep_same}"
    else
        : > "${sheep_same}"
    fi

    # E10/E11 are intentionally excluded here.
    local sheep_other; sheep_other=$(mktemp)
    local ofiles=()
    local s
    for s in E1 E2 E3 E4 E5 E6 E7 E8 E9; do
        [[ "${s}" == "${state}" ]] && continue
        local f="${SHEEP_REF_DIR}/${sheep_t}_${s}.v3.bed"
        [[ -s "${f}" ]] && ofiles+=("${f}")
    done
    if [[ ${#ofiles[@]} -gt 0 ]]; then
        cat "${ofiles[@]}" | bedtools sort -i - | bedtools merge -i - > "${sheep_other}"
    else
        : > "${sheep_other}"
    fi

    if is_promoter "${state}"; then
        local mid_all; mid_all=$(mktemp)
        to_midpoint "${human_v3}" > "${mid_all}"

        local mid_sf; mid_sf=$(mktemp)
        bedtools intersect -a "${mid_all}" -b "${sheep_same}" \
            -f 1.0 -wa | sort -u > "${mid_sf}"
        recover_by_name "${mid_sf}" "${human_v3}" > "${sf}"

        local non_sf; non_sf=$(mktemp)
        bedtools intersect -a "${human_v3}" -b "${sf}" -v > "${non_sf}"

        local mid_ns; mid_ns=$(mktemp)
        to_midpoint "${non_sf}" > "${mid_ns}"

        local mid_sd; mid_sd=$(mktemp)
        bedtools intersect -a "${mid_ns}" -b "${sheep_other}" \
            -f 1.0 -wa | sort -u > "${mid_sd}"
        recover_by_name "${mid_sd}" "${non_sf}" > "${sd}"

        bedtools intersect -a "${non_sf}" -b "${sd}" -v > "${so}"
        rm -f "${mid_all}" "${mid_sf}" "${non_sf}" "${mid_ns}" "${mid_sd}"
    else
        bedtools intersect -a "${human_v3}" -b "${sheep_same}" \
            -f "${OVERLAP_FRAC}" -wa | sort -u > "${sf}"

        local non_sf; non_sf=$(mktemp)
        bedtools intersect -a "${human_v3}" -b "${sheep_same}" \
            -f "${OVERLAP_FRAC}" -v > "${non_sf}"

        bedtools intersect -a "${non_sf}" -b "${sheep_other}" \
            -f "${OVERLAP_FRAC}" -wa | sort -u > "${sd}"

        bedtools intersect -a "${non_sf}" -b "${sheep_other}" \
            -f "${OVERLAP_FRAC}" -v > "${so}"
        rm -f "${non_sf}"
    fi

    rm -f "${sheep_same}" "${sheep_other}"
    echo "[DONE] ${human_t}->${sheep_t} ${state}: sf=$(wc -l < "${sf}") sd=$(wc -l < "${sd}") so=$(wc -l < "${so}")"
}
export -f process_classify

log "===== Step 2: rebuild E1-E9 classification ====="
TASKLIST="${TMP_DIR}/class_tasks.txt"
: > "${TASKLIST}"
for i in "${!HUMAN_TISSUES[@]}"; do
    for state in "${STATES[@]}"; do
        echo "${HUMAN_TISSUES[$i]} ${SHEEP_TISSUES[$i]} ${state}"
    done
done > "${TASKLIST}"

log "classification tasks: $(wc -l < "${TASKLIST}")"
xargs -P "${NCPU}" -L 1 bash -c 'process_classify "$@"' _ < "${TASKLIST}"

log "===== Step 3: summary ====="
SUMMARY="${WORKDIR}/per_tissue_classification.tsv"
{
    echo -e "human_tissue\tsheep_tissue\tstate\tsfCRE\tsdCRE\tsoCRE\tssCRE\ttotal"
    for i in "${!HUMAN_TISSUES[@]}"; do
        ht="${HUMAN_TISSUES[$i]}"
        st="${SHEEP_TISSUES[$i]}"
        for state in "${STATES[@]}"; do
            pf="${CLASS_DIR}/${ht}_${st}_${state}"
            sf=0; sd=0; so=0; ss=0
            [[ -s "${pf}.sfCRE.bed" ]] && sf=$(wc -l < "${pf}.sfCRE.bed")
            [[ -s "${pf}.sdCRE.bed" ]] && sd=$(wc -l < "${pf}.sdCRE.bed")
            [[ -s "${pf}.soCRE.bed" ]] && so=$(wc -l < "${pf}.soCRE.bed")
            unmap="${LIFTED_DIR}/${ht}_${state}.unmapped"
            [[ -s "${unmap}" ]] && ss=$(grep -cv "^#" "${unmap}" 2>/dev/null || echo 0)
            total=$((sf + sd + so + ss))
            echo -e "${ht}\t${st}\t${state}\t${sf}\t${sd}\t${so}\t${ss}\t${total}"
        done
    done
} > "${SUMMARY}"
log "summary: ${SUMMARY}"

if [[ -s "${SRC_SUMMARY}" ]]; then
    DELTA="${WORKDIR}/per_tissue_classification_delta_vs_source.tsv"
    awk -F'\t' 'BEGIN{OFS="\t"}
        NR==FNR {
            if (FNR > 1 && $3 != "E10") {
                k=$1 SUBSEP $2 SUBSEP $3
                osf[k]=$4+0; osd[k]=$5+0; oso[k]=$6+0; oss[k]=$7+0; ot[k]=$8+0
            }
            next
        }
        FNR==1 {
            print "human_tissue","sheep_tissue","state",
                  "old_sfCRE","new_sfCRE","delta_sfCRE",
                  "old_sdCRE","new_sdCRE","delta_sdCRE",
                  "old_soCRE","new_soCRE","delta_soCRE",
                  "old_ssCRE","new_ssCRE","delta_ssCRE",
                  "old_total","new_total","delta_total"
            next
        }
        {
            k=$1 SUBSEP $2 SUBSEP $3
            print $1,$2,$3,
                  osf[k],$4+0,($4+0)-osf[k],
                  osd[k],$5+0,($5+0)-osd[k],
                  oso[k],$6+0,($6+0)-oso[k],
                  oss[k],$7+0,($7+0)-oss[k],
                  ot[k],$8+0,($8+0)-ot[k]
        }' "${SRC_SUMMARY}" "${SUMMARY}" > "${DELTA}"
    log "delta vs source: ${DELTA}"
fi

log "===== ALL DONE ====="
