#!/usr/bin/env bash

set -euo pipefail

# Figure 2j-2k, step 5: collect basic Juicer QC and contact-distance decay.

JUICER_ROOT="${JUICER_ROOT:-/vol2/zhangshiwen/juicer}"
WORK_ROOT="${WORK_ROOT:-${JUICER_ROOT}/work}"
SAMPLE_FILE="${SAMPLE_FILE:-${JUICER_ROOT}/tissue.txt}"
OUTDIR="${OUTDIR:-${WORK_ROOT}/juicer_qc}"

SUMMARY="${OUTDIR}/juicer_basic_qc.tsv"
DECAY="${OUTDIR}/juicer_distance_decay.tsv"

[[ -s "${SAMPLE_FILE}" ]] || {
  echo "[ERROR] Sample file was not found: ${SAMPLE_FILE}" >&2
  exit 1
}
mkdir -p "${OUTDIR}"

printf '%s\n' \
  $'sample\tgroup\traw_read_pairs\tvalid_pairs_merged_nodups\tdups\topt_dups\ttotal_dups\tvalid_rate_vs_raw\tdup_rate_vs_valid_plus_dups\tcis_pairs\ttrans_pairs\tcis_fraction\ttrans_fraction\tcis_trans_ratio\tcis_short_lt20kb\tcis_long_ge20kb\tlong_short_ratio\tlong_cis_fraction\thas_inter30_hic\thas_inter30_loops\thas_inter30_domains' \
  > "${SUMMARY}"
printf 'sample\tgroup\tbin\tn_pairs\tmean_distance\n' > "${DECAY}"

get_group() {
  local sample_id="$1"
  case "${sample_id}" in
    HIC_AbdominalSubcutaneousAdipose_*) echo "Adipose" ;;
    HIC_blood_*) echo "Blood" ;;
    HIC_EarSkinFibroblasts_*) echo "EarSkinFibroblasts" ;;
    HIC_liver_*) echo "Liver" ;;
    HIC_Common_tissues_*) echo "Common_tissues" ;;
    *) echo "Other" ;;
  esac
}

count_lines() {
  local input_file="$1"
  if [[ ! -s "${input_file}" ]]; then
    echo 0
  elif [[ "${input_file}" == *.gz ]]; then
    zcat "${input_file}" | wc -l
  else
    wc -l < "${input_file}"
  fi
}

count_raw_pairs() {
  local fastq_dir="$1"
  local r1
  r1="$(
    find "${fastq_dir}" -maxdepth 1 -type f \
      \( -name '*_R1*.fastq' -o -name '*_R1*.fastq.gz' \
         -o -name '*_1.fastq' -o -name '*_1.fastq.gz' \) \
      -print 2>/dev/null | LC_ALL=C sort | head -n 1 || true
  )"
  if [[ -z "${r1}" || ! -s "${r1}" ]]; then
    echo 0
  elif [[ "${r1}" == *.gz ]]; then
    zcat "${r1}" | awk 'END { printf "%d\n", NR / 4 }'
  else
    awk 'END { printf "%d\n", NR / 4 }' "${r1}"
  fi
}

while read -r sample_id _; do
  [[ -z "${sample_id}" ]] && continue

  group="$(get_group "${sample_id}")"
  sample_dir="${WORK_ROOT}/${sample_id}"
  aligned="${sample_dir}/aligned"
  merged="${aligned}/merged_nodups.txt"
  dups="${aligned}/dups.txt"
  opt_dups="${aligned}/opt_dups.txt"

  if [[ ! -d "${sample_dir}" ]]; then
    echo "[WARN] Sample directory was not found: ${sample_dir}" >&2
    continue
  fi

  raw_pairs="$(count_raw_pairs "${sample_dir}/fastq")"
  valid_pairs="$(count_lines "${merged}")"
  n_dups="$(count_lines "${dups}")"
  n_opt_dups="$(count_lines "${opt_dups}")"
  total_dups="$((n_dups + n_opt_dups))"

  has_hic=0
  has_loops=0
  has_domains=0
  [[ -s "${aligned}/inter_30.hic" ]] && has_hic=1
  [[ -s "${aligned}/inter_30_loops/merged_loops.bedpe" ]] && has_loops=1
  [[ -s "${aligned}/inter_30_contact_domains/10000_blocks.bedpe" ]] \
    && has_domains=1

  if [[ ! -s "${merged}" ]]; then
    printf '%s\t%s\t%s\t0\t%s\t%s\t%s\tNA\tNA\t0\t0\tNA\tNA\tNA\t0\t0\tNA\tNA\t%s\t%s\t%s\n' \
      "${sample_id}" "${group}" "${raw_pairs}" "${n_dups}" \
      "${n_opt_dups}" "${total_dups}" "${has_hic}" "${has_loops}" \
      "${has_domains}" >> "${SUMMARY}"
    continue
  fi

  awk \
    -v sample="${sample_id}" \
    -v group="${group}" \
    -v raw="${raw_pairs}" \
    -v valid="${valid_pairs}" \
    -v dups="${n_dups}" \
    -v opt_dups="${n_opt_dups}" \
    -v total_dups="${total_dups}" \
    -v has_hic="${has_hic}" \
    -v has_loops="${has_loops}" \
    -v has_domains="${has_domains}" '
      BEGIN { OFS = "\t" }
      {
        chr1 = $2
        pos1 = $3
        chr2 = $6
        pos2 = $7
        if (chr1 == chr2) {
          cis++
          distance = pos2 - pos1
          if (distance < 0) distance = -distance
          if (distance < 20000) short_cis++
          else long_cis++
        } else {
          trans++
        }
      }
      END {
        valid_rate = raw > 0 ? valid / raw : "NA"
        duplicate_rate = valid + total_dups > 0 \
          ? total_dups / (valid + total_dups) : "NA"
        total_contacts = cis + trans
        cis_fraction = total_contacts > 0 ? cis / total_contacts : "NA"
        trans_fraction = total_contacts > 0 ? trans / total_contacts : "NA"
        cis_trans_ratio = trans > 0 ? cis / trans : "Inf"
        long_short_ratio = short_cis > 0 ? long_cis / short_cis : "Inf"
        long_cis_fraction = cis > 0 ? long_cis / cis : "NA"
        print sample, group, raw, valid, dups, opt_dups, total_dups, \
          valid_rate, duplicate_rate, cis, trans, cis_fraction, \
          trans_fraction, cis_trans_ratio, short_cis, long_cis, \
          long_short_ratio, long_cis_fraction, has_hic, has_loops, \
          has_domains
      }
    ' "${merged}" >> "${SUMMARY}"

  awk -v sample="${sample_id}" -v group="${group}" '
    BEGIN { OFS = "\t" }
    $2 == $6 {
      distance = $7 - $3
      if (distance < 0) distance = -distance
      if (distance < 10000) bin = "0_10kb"
      else if (distance < 50000) bin = "10_50kb"
      else if (distance < 100000) bin = "50_100kb"
      else if (distance < 500000) bin = "100_500kb"
      else if (distance < 1000000) bin = "500kb_1Mb"
      else if (distance < 5000000) bin = "1_5Mb"
      else bin = "gt5Mb"
      count[bin]++
      distance_sum[bin] += distance
    }
    END {
      for (bin in count) {
        print sample, group, bin, count[bin], distance_sum[bin] / count[bin]
      }
    }
  ' "${merged}" >> "${DECAY}"
done < "${SAMPLE_FILE}"

echo "[DONE] Basic QC: ${SUMMARY}"
echo "[DONE] Distance decay: ${DECAY}"
