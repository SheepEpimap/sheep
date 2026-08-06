#!/usr/bin/env bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
set -euo pipefail

# =========================
# =========================
IN_BAM_DIR="/vol2/zhangshiwen/atac_public/fastq/vail_chrombpnt/results/05_bam"
IN_PEAK_DIR="/vol2/zhangshiwen/atac_public/fastq/vail_chrombpnt/results/08_peaks"

OUT_BASE="/data/home/sczd644/run/zsw_chrombpnet/public/ATAC_bams"
BLACKLIST="/vol2/zhangshiwen/blacklist/output/sheep_blacklist.bed"

THREADS=12
GSIZE="2628104905"        # macs2   genome size(  peak)
EFF_GSIZE="2615608860"    # bamCoverage RPGC   effectiveGenomeSize

MARK="ATAC"
MARK_TYPE="narrowPeak"

CHR_LIST=(chr{1..26} chrX)
,
mkdir -p "${OUT_BASE}"

# =========================
# =========================
cd "${IN_BAM_DIR}"

tissue_list=$(
  ls ATAC_*_Rep*.last.bam 2>/dev/null \
  | sed -E 's/^ATAC_(.+)_Rep[0-9]+\.last\.bam$/\1/' \
  | sort -u
)

if [[ -z "${tissue_list}" ]]; then
  echo "[ERROR] No ATAC_*_Rep*.last.bam in ${IN_BAM_DIR}" >&2
  exit 1
fi

# =========================
# =========================
for tissue in ${tissue_list}; do
  echo "===================="
  echo "[TISSUE] ${tissue}"
  echo "===================="

  reps=$(
    ls "ATAC_${tissue}_Rep"*.last.bam 2>/dev/null \
    | sed -E "s/^ATAC_${tissue}_Rep([0-9]+)\.last\.bam$/\1/" \
    | sort -n
  )

  for r1 in ${reps}; do
    if (( r1 % 2 == 0 )); then
      continue
    fi
    r2=$((r1 + 1))

    bam1="${IN_BAM_DIR}/ATAC_${tissue}_Rep${r1}.last.bam"
    bam2="${IN_BAM_DIR}/ATAC_${tissue}_Rep${r2}.last.bam"

    if [[ ! -f "${bam1}" ]]; then
      continue
    fi
    if [[ ! -f "${bam2}" ]]; then
      echo "[WARN] Missing pair: ATAC_${tissue}_Rep${r2}.last.bam  -> skip Rep${r1}+Rep${r2}" >&2
      continue
    fi

    p1=$(printf "%02d" "${r1}")
    p2=$(printf "%02d" "${r2}")

    pair_tag="Rep${p1}-${p2}"
    pair_dir="${OUT_BASE}/${tissue}/${pair_tag}"
    data_dir="${pair_dir}/data"
    idr_dir="${pair_dir}/idr_output"
    log_dir="${pair_dir}/logs"

    mkdir -p "${data_dir}" "${idr_dir}" "${log_dir}"

    # -------------------------
    # -------------------------
    merged_unsorted="${pair_dir}/merged_unsorted.bam"
    merged_sorted="${data_dir}/merged_sorted.bam"
    merged_bam="${data_dir}/merged.bam"

    if [[ -s "${merged_bam}" && -s "${merged_bam}.bai" ]]; then
      echo "[SKIP] merged.bam exists: ${merged_bam}"
    else
      echo "[DO] merge/sort/index/filter-chr: ${tissue} ${pair_tag}"

      samtools merge -@ "${THREADS}" -f "${merged_unsorted}" "${bam1}" "${bam2}"
      samtools sort  -@ "${THREADS}" -o "${merged_sorted}" "${merged_unsorted}"
      samtools index -@ "${THREADS}" "${merged_sorted}"

      samtools view -b "${merged_sorted}" "${CHR_LIST[@]}" > "${merged_bam}"
      samtools index -@ "${THREADS}" "${merged_bam}"

      rm -f "${merged_unsorted}" "${merged_sorted}" "${merged_sorted}.bai"
    fi

    # -------------------------
    # -------------------------
    peak1="${IN_PEAK_DIR}/ATAC_${tissue}_Rep${r1}_peaks.narrowPeak"
    peak2="${IN_PEAK_DIR}/ATAC_${tissue}_Rep${r2}_peaks.narrowPeak"

    idr_out="${idr_dir}/${MARK}_${tissue}_${pair_tag}.idr.txt"
    idr_log="${idr_dir}/${MARK}_${tissue}_${pair_tag}.log"

    if [[ -s "${idr_out}" ]]; then
      echo "[SKIP] IDR exists: ${idr_out}"
    else
      if [[ ! -s "${peak1}" || ! -s "${peak2}" ]]; then
        echo "[WARN] Missing peak file(s), skip IDR/overlap: ${peak1} or ${peak2}" >&2
        continue
      fi
      echo "[DO] IDR: ${peak1} vs ${peak2}"
      idr --samples "${peak1}" "${peak2}" \
          --input-file-type "${MARK_TYPE}" \
          --output-file "${idr_out}" \
          --plot \
          --log-output-file "${idr_log}" \
          > "${log_dir}/idr.stdout.log" 2>&1
    fi

    # -------------------------
    # -------------------------
    overlap_bed="${idr_dir}/${MARK}_${tissue}_${pair_tag}.overlap_peak.bed"
    peaks_no_bl="${data_dir}/peaks_no_blacklist.bed"
    peaks_obs="${data_dir}/peaks_no_blacklist.observed.bed"

    if [[ -s "${peaks_no_bl}" ]]; then
      echo "[SKIP] peaks_no_blacklist exists: ${peaks_no_bl}"
    else
      echo "[DO] overlap (IDR & rep1 & rep2) + remove blacklist"

      tmp1="${idr_dir}/tmp.${pair_tag}.idr_vs_rep1.txt"
      tmp2="${idr_dir}/tmp.${pair_tag}.idr_vs_rep2.txt"

      bedtools intersect -wo -a "${idr_out}" -b "${peak1}" \
      | awk 'BEGIN{FS=OFS="\t"}{
              bcols=10;            # narrowPeak   10
              a=NF-(bcols+1);      # A  (NF = A + B(10) + overlap(1))
              s1=$3-$2;            # A  (A   start/end  2/3 )
              s2=$(a+3)-$(a+2);    # B  (B   start/end   A  2/3 )
              ov=$NF;              # overlap bp
              if ( (ov/s1>=0.5) || (ov/s2>=0.5) ){
                for(i=1;i<=a;i++){printf "%s%s", $i, (i==a?ORS:OFS)}
              }
            }' \
      | sort -u > "${tmp1}"

      bedtools intersect -wo -a "${tmp1}" -b "${peak2}" \
      | awk 'BEGIN{FS=OFS="\t"}{
              bcols=10;
              a=NF-(bcols+1);
              s1=$3-$2;
              s2=$(a+3)-$(a+2);
              ov=$NF;
              if ( (ov/s1>=0.5) || (ov/s2>=0.5) ){
                for(i=1;i<=a;i++){printf "%s%s", $i, (i==a?ORS:OFS)}
              }
            }' \
      | sort -u \
      | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-6]|X)$/ {print $0}' \
      > "${tmp2}"

      cp "${tmp2}" "${overlap_bed}"

      bedtools intersect -v -a "${overlap_bed}" -b "${BLACKLIST}" \
      | cut -f1-10 > "${peaks_no_bl}"

      rm -f "${tmp1}" "${tmp2}"
    fi

    # -------------------------
    # D) merged.bw + observed bed
    # -------------------------
    merged_bw="${data_dir}/merged.bw"
    if [[ -s "${merged_bw}" ]]; then
      echo "[SKIP] merged.bw exists: ${merged_bw}"
    else
      echo "[DO] bamCoverage -> ${merged_bw}"
      bamCoverage --bam "${merged_bam}" \
        -o "${merged_bw}" \
        --binSize 10 \
        --normalizeUsing RPGC --effectiveGenomeSize "${EFF_GSIZE}" \
        --extendReads -p "${THREADS}" \
        > "${log_dir}/bamCoverage.log" 2>&1
    fi

    if [[ -s "${peaks_obs}" ]]; then
      echo "[SKIP] observed bed exists: ${peaks_obs}"
    else
      echo "[DO] make observed bed (add peakID)"
      awk -v OFS="\t" '{print $1,$2,$3,$1"_"$2"_"$3}' "${peaks_no_bl}" > "${peaks_obs}"
    fi

    echo "[OK] ${tissue} ${pair_tag} done -> ${pair_dir}"
  done
done
