#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

############################################
############################################
RESULTS_BASE="/vol2/zhangshiwen/rrbs/wgbs_bismark/results"
SAMPLE_LIST="/vol2/zhangshiwen/rrbs/wgbs_bismark/samples.tsv"
META_TSV="/vol2/zhangshiwen/rrbs/md5sum.txt"
STATE_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/All_chromatin_state"
OUTDIR="/vol2/zhangshiwen/rrbs/wgbs_bismark/state_methylation_by_sample_tissue"

mkdir -p "${OUTDIR}"/tmp
mkdir -p "${OUTDIR}"/tmp/tissue_state_beds
mkdir -p "${OUTDIR}"/per_sample
mkdir -p "${OUTDIR}"/plot
mkdir -p "${OUTDIR}"/log

############################################
############################################
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found"; exit 1; }
command -v zcat >/dev/null 2>&1 || { echo "[ERROR] zcat not found"; exit 1; }

############################################
############################################
: > "${OUTDIR}/log/skipped_existing.log"
: > "${OUTDIR}/log/processed_new.log"
: > "${OUTDIR}/log/missing_cov.log"
: > "${OUTDIR}/log/missing_tissue_in_meta.log"
: > "${OUTDIR}/log/missing_state_bed.log"
: > "${OUTDIR}/log/empty_output.log"
: > "${OUTDIR}/log/missing_result_when_merge.log"

############################################
############################################
echo "[INFO] Building tissue-level combined state BEDs ..."

rm -f "${OUTDIR}"/tmp/tissue_state_beds/*

for f in "${STATE_DIR}"/*_E*.bed; do
  [ -e "$f" ] || continue

  base=$(basename "$f")
  tissue="${base%_E*.bed}"   #   muscle
  state="${base##*_}"        #   E1.bed
  state="${state%.bed}"      #   E1

  awk -v OFS='\t' -v st="$state" '
    BEGIN{FS=OFS="\t"}
    $0 !~ /^#/ && NF >= 3 {
      s = $2
      e = $3
      if (s > e) {
        tmp = s; s = e; e = tmp
      }
      print $1, s, e, st
    }
  ' "$f" >> "${OUTDIR}/tmp/tissue_state_beds/${tissue}.unsorted.bed"
done

for uns in "${OUTDIR}"/tmp/tissue_state_beds/*.unsorted.bed; do
  [ -e "$uns" ] || continue
  sorted="${uns%.unsorted.bed}.states.sorted.bed"
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$uns" > "$sorted"
  rm -f "$uns"
done

echo "[INFO] Tissue state BEDs prepared."

############################################
############################################
mapfile -t SAMPLES < <(
  awk 'BEGIN{FS=OFS="\t"} NR>1 && $1!="" {print $1}' "$SAMPLE_LIST" \
  | sed 's/\r$//' \
  | sort -u
)

echo "[INFO] Total samples in list: ${#SAMPLES[@]}"

############################################
############################################
declare -A SAMPLE2TISSUE

while IFS=$'\t' read -r run_accession fastq_aspera fastq_md5 tissue_type; do
  [[ "$run_accession" == "run_accession" ]] && continue
  [[ -z "${run_accession}" ]] && continue
  SAMPLE2TISSUE["$run_accession"]="$tissue_type"
done < "$META_TSV"

############################################
# sample  sample_tissue  state_tissue  state
# methylated_count  unmethylated_count  total_count  cpg_sites
# weighted_meth_pct
############################################
for sample in "${SAMPLES[@]}"; do
  echo "[INFO] Checking ${sample} ..."

  OUT_SAMPLE="${OUTDIR}/per_sample/${sample}.sample_state_methylation.tsv"

  if [[ -s "$OUT_SAMPLE" ]]; then
    echo "[INFO] Existing result found, skip: ${sample}"
    echo -e "${sample}\t${OUT_SAMPLE}" >> "${OUTDIR}/log/skipped_existing.log"
    continue
  fi

  tissue="${SAMPLE2TISSUE[$sample]:-}"
  if [[ -z "$tissue" ]]; then
    echo "[WARN] No tissue found in md5sum.txt for ${sample}"
    echo -e "${sample}\tNO_TISSUE_IN_META" >> "${OUTDIR}/log/missing_tissue_in_meta.log"
    continue
  fi

  used_tissue="$tissue"
  STATE_BED="${OUTDIR}/tmp/tissue_state_beds/${used_tissue}.states.sorted.bed"

  if [[ ! -s "$STATE_BED" ]]; then
    alt1="${tissue//_/-}"
    alt2="${tissue// /-}"

    if [[ -s "${OUTDIR}/tmp/tissue_state_beds/${alt1}.states.sorted.bed" ]]; then
      used_tissue="$alt1"
      STATE_BED="${OUTDIR}/tmp/tissue_state_beds/${used_tissue}.states.sorted.bed"
    elif [[ -s "${OUTDIR}/tmp/tissue_state_beds/${alt2}.states.sorted.bed" ]]; then
      used_tissue="$alt2"
      STATE_BED="${OUTDIR}/tmp/tissue_state_beds/${used_tissue}.states.sorted.bed"
    fi
  fi

  if [[ ! -s "$STATE_BED" ]]; then
    echo "[WARN] No chromatin state BED for sample=${sample}, tissue=${tissue}"
    echo -e "${sample}\t${tissue}\tNO_MATCHED_STATE_BED" >> "${OUTDIR}/log/missing_state_bed.log"
    continue
  fi

  cov=$(find "${RESULTS_BASE}/${sample}" -maxdepth 3 -type f -name "*.deduplicated.bismark.cov.gz" | head -n 1 || true)
  if [[ -z "$cov" ]]; then
    cov=$(find "${RESULTS_BASE}/${sample}" -maxdepth 3 -type f -name "*.bismark.cov.gz" | head -n 1 || true)
  fi

  if [[ -z "$cov" ]]; then
    echo "[WARN] No cov.gz found for ${sample}"
    echo -e "${sample}\t${tissue}\tNO_COV" >> "${OUTDIR}/log/missing_cov.log"
    continue
  fi

  echo "[INFO] Processing new sample: ${sample}"
  echo "[INFO] Tissue in meta: ${tissue}"
  echo "[INFO] State tissue used: ${used_tissue}"
  echo "[INFO] Cov file: ${cov}"

  COVBED="${OUTDIR}/tmp/${sample}.cov.bed"
  SUMMARY_TMP="${OUTDIR}/tmp/${sample}.state_summary.tmp"

  ############################################
  ############################################
  zcat "$cov" \
    | awk 'BEGIN{FS=OFS="\t"} NF>=6 {s=$2-1; if(s<0)s=0; print $1, s, $3, $4, $5, $6}' \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$COVBED"

  ############################################
  ############################################
  bedtools intersect \
    -sorted \
    -a "$COVBED" \
    -b "$STATE_BED" \
    -wa -wb \
  | awk 'BEGIN{FS=OFS="\t"}
      {
        # A(cov):   1 chr, 2 start, 3 end, 4 meth_pct, 5 meth_count, 6 unmeth_count
        # B(state): 7 chr, 8 start, 9 end, 10 state
        state = $10
        meth[state]   += $5
        unmeth[state] += $6
        cpg[state]    += 1
      }
      END{
        for (s in meth) {
          total = meth[s] + unmeth[s]
          if (total > 0) {
            pct = 100 * meth[s] / total
          } else {
            pct = "NA"
          }
          print s, meth[s], unmeth[s], total, cpg[s], pct
        }
      }' \
  | LC_ALL=C sort -k1,1V > "$SUMMARY_TMP"

  ############################################
  ############################################
  awk -v OFS='\t' -v smp="$sample" -v tis="$tissue" -v stis="$used_tissue" '
    NR==FNR {
      states[$4] = 1
      next
    }
    {
      seen[$1]   = 1
      meth[$1]   = $2
      unmeth[$1] = $3
      total[$1]  = $4
      cpg[$1]    = $5
      pct[$1]    = $6
    }
    END {
      for (i=1; i<=20; i++) {
        s = "E" i
        if (s in states) {
          if (s in seen) {
            print smp, tis, stis, s, meth[s], unmeth[s], total[s], cpg[s], pct[s]
          } else {
            print smp, tis, stis, s, 0, 0, 0, 0, "NA"
          }
        }
      }
    }' "$STATE_BED" "$SUMMARY_TMP" > "$OUT_SAMPLE"

  if [[ -s "$OUT_SAMPLE" ]]; then
    echo "[INFO] Done: ${sample}"
    echo -e "${sample}\t${tissue}\t${used_tissue}\t${OUT_SAMPLE}" >> "${OUTDIR}/log/processed_new.log"
  else
    echo "[WARN] Empty output for ${sample}"
    echo -e "${sample}\t${tissue}\tEMPTY_OUTPUT" >> "${OUTDIR}/log/empty_output.log"
    rm -f "$OUT_SAMPLE"
  fi

  rm -f "$COVBED" "$SUMMARY_TMP"
done

############################################
############################################
MASTER="${OUTDIR}/all_samples_sample_state_weighted_methylation.tsv"

echo -e "sample\tsample_tissue\tstate_tissue\tstate\tmethylated_count\tunmethylated_count\ttotal_count\tcpg_sites\tweighted_meth_pct" > "$MASTER"

merged_n=0
for sample in "${SAMPLES[@]}"; do
  OUT_SAMPLE="${OUTDIR}/per_sample/${sample}.sample_state_methylation.tsv"
  if [[ -s "$OUT_SAMPLE" ]]; then
    cat "$OUT_SAMPLE" >> "$MASTER"
    merged_n=$((merged_n + 1))
  else
    echo -e "${sample}\tNO_PER_SAMPLE_RESULT" >> "${OUTDIR}/log/missing_result_when_merge.log"
  fi
done

echo "[INFO] Merge done."
echo "[INFO] Samples merged: ${merged_n}"
echo "[INFO] Master table: ${MASTER}"
echo "[INFO] Finished."
