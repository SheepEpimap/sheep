#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] line $LINENO, exit code $?" >&2' ERR
shopt -s nullglob

# =========================================================
# Build LDSC annotation beds from two sources:
#
# Source 1: human side CRE (CRE_TSR_E5_hg19/)
#   cols: chr start end hg38_key CRE_class TSR_label TSR_bp
#   slice by (CRE_class, TSR_label)
#
# Source 2: sheep projection mapped back to human hg19
#           (sf_sheepE5_TSRlabeled/hg19/)
#   cols: chr start end hg38_key human_TSR_label
#   CRE_class is implicitly sfCRE (by construction),
#   slice by TSR_label only
#
# Output: one 3-col bed per slice, LDSC-ready
# =========================================================

: "${WORKDIR:?Set E1E9_WORKDIR/WORKDIR in the paths configuration}"

# ── Source 1: human CRE ─────────────────────────────────
HUMAN_DIR="${WORKDIR}/CRE_TSR_E5_hg19"
OUT_HUMAN="${WORKDIR}/LDSC_human_CRE_hg19"
mkdir -p "${OUT_HUMAN}"

# ── Source 2: sheep-projected sfCRE ─────────────────────
PROJ_DIR="${WORKDIR}/sf_sheepE5_TSRlabeled/hg19"
OUT_PROJ="${WORKDIR}/LDSC_sheep_projection_hg19"
mkdir -p "${OUT_PROJ}"

CRE_CLASSES=(sfCRE sdCRE soCRE ssCRE)
TSR_BROAD=(all_common broad intermediate tissue_specific)
TISSUES=(Adipose Colon Cortex Heart Liver Lung Muscle Ovary Sintest Spleen Stomach Testis)

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# =========================================================
# Source 1: human CRE — slice by (CRE_class, TSR_label)
# =========================================================
log "===== Source 1: human CRE slices ====="

# Build one pooled BED per CRE class and broad TSR label across 12 tissue pairs.
for cre in "${CRE_CLASSES[@]}"; do
    files=( "${HUMAN_DIR}"/*.${cre}.hg19.TSR.bed )
    (( ${#files[@]} )) || { log "  [WARN] no ${cre} files"; continue; }

    # all_common / broad / intermediate
    for L in all_common broad intermediate; do
        awk -v L="$L" 'BEGIN{FS=OFS="\t"} $6==L {print $1,$2,$3}' "${files[@]}" \
          | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u \
          | bedtools merge -i - > "${OUT_HUMAN}/${cre}_${L}.bed"
    done

    # tissue_specific pooled (any ts_*)
    awk 'BEGIN{FS=OFS="\t"} $6 ~ /^ts_/ {print $1,$2,$3}' "${files[@]}" \
      | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u \
      | bedtools merge -i - > "${OUT_HUMAN}/${cre}_tissue_specific.bed"

    # per-tissue ts (CTS LDSC)
    for T in "${TISSUES[@]}"; do
        # Collect rows labeled ts_${T} across all tissue pairs.
        # A CRE labeled ts_Adipose retains that label in every tissue-pair file.
        awk -v L="ts_${T}" 'BEGIN{FS=OFS="\t"} $6==L {print $1,$2,$3}' "${files[@]}" \
          | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u \
          | bedtools merge -i - > "${OUT_HUMAN}/${cre}_ts_${T}.bed"
    done
done

# =========================================================
# Source 2: sheep projection — slice by TSR_label only
#   (CRE_class is sfCRE by construction)
# =========================================================
log "===== Source 2: sheep projection slices ====="

files=( "${PROJ_DIR}"/*.hg19.bed )
if (( ${#files[@]} )); then
    for L in all_common broad intermediate; do
        awk -v L="$L" 'BEGIN{FS=OFS="\t"} $5==L {print $1,$2,$3}' "${files[@]}" \
          | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u \
          | bedtools merge -i - > "${OUT_PROJ}/sf_proj_${L}.bed"
    done

    awk 'BEGIN{FS=OFS="\t"} $5 ~ /^ts_/ {print $1,$2,$3}' "${files[@]}" \
      | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u \
      | bedtools merge -i - > "${OUT_PROJ}/sf_proj_tissue_specific.bed"

    for T in "${TISSUES[@]}"; do
        awk -v L="ts_${T}" 'BEGIN{FS=OFS="\t"} $5==L {print $1,$2,$3}' "${files[@]}" \
          | LC_ALL=C sort -k1,1V -k2,2n -k3,3n -u \
          | bedtools merge -i - > "${OUT_PROJ}/sf_proj_ts_${T}.bed"
    done
else
    log "  [WARN] no projection files"
fi

# =========================================================
# Summary
# =========================================================
log "===== Summary ====="
SUM="${WORKDIR}/LDSC_annotation_summary.tsv"
{
echo -e "source\tslice\tintervals\ttotal_bp"
for f in "${OUT_HUMAN}"/*.bed "${OUT_PROJ}"/*.bed; do
    [[ -s "$f" ]] || { n=0; bp=0; }
    n=$(wc -l < "$f")
    bp=$(awk '{s+=$3-$2}END{print s+0}' "$f")
    dir=$(basename "$(dirname "$f")")
    name=$(basename "$f" .bed)
    if [[ "$dir" == "LDSC_human_CRE_hg19" ]]; then
        src="human_CRE"
    else
        src="sheep_projection"
    fi
    echo -e "${src}\t${name}\t${n}\t${bp}"
done
} > "${SUM}"

log "Summary: ${SUM}"
log "Human CRE slices: ${OUT_HUMAN}"
log "Sheep projection slices: ${OUT_PROJ}"
log "===== ALL DONE ====="
