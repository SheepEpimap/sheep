#!/usr/bin/env bash
set -euo pipefail

: "${SRC:?Set E1E9_WORKDIR/SRC in the paths configuration}"
: "${GWAS_BASE:?Set GWAS_BASE in the paths configuration}"
DST="${DST:-${GWAS_BASE}/beds}"

mkdir -p "${DST}"

# A1: 12 BED files (four human CRE classes by three TSR classes).
cd "${SRC}/LDSC_human_CRE_hg19/"
for c in sfCRE sdCRE soCRE ssCRE; do
    cat "${c}_all_common.bed" "${c}_broad.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > "${DST}/${c}_common.bed"
    cp "${c}_intermediate.bed" "${DST}/"
    cp "${c}_tissue_specific.bed" "${DST}/"
done

# A2: 43 BED files after excluding four Testis files and sfCRE_Muscle.
for c in sfCRE sdCRE soCRE ssCRE; do
    for t in Adipose Colon Cortex Heart Liver Lung Muscle Ovary Sintest Spleen Stomach; do
        if [ "$c" = "sfCRE" ] && [ "$t" = "Muscle" ]; then continue; fi
        cp "${c}_ts_${t}.bed" "${DST}/"
    done
done

# B1: 3 BED
cd "${SRC}/LDSC_sheep_projection_hg19/"
cat sf_proj_all_common.bed sf_proj_broad.bed | sort -k1,1 -k2,2n | bedtools merge -i - > "${DST}/sf_proj_common.bed"
cp sf_proj_intermediate.bed "${DST}/"
cp sf_proj_tissue_specific.bed "${DST}/"

# B2: 10 BED files after excluding Testis and Muscle.
for t in Adipose Colon Cortex Heart Liver Lung Ovary Sintest Spleen Stomach; do
    cp "sf_proj_ts_${t}.bed" "${DST}/"
done

# Validate that exactly 68 BED files were collected.
ls "${DST}"/*.bed | wc -l
