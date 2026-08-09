#!/usr/bin/env bash
#SBATCH -p smp
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${SCRIPT_DIR}/00_config.sh"; source "${SCRIPT_DIR}/00_lib.sh"
as_activate_conda "${ASE_ENV}"
sample_id="${1:-}"; [[ -n "${sample_id}" ]] || as_die "Usage: 15_prepare_ase_counts.sh SAMPLE_ID"
as_require_cmds python
gatk_table="${ASE_TABLE_DIR}/${sample_id}.output.table"; out="${ASE_COUNT_DIR}/${sample_id}.ase_counts.tsv"; as_require_file "${gatk_table}" "GATK ASE table"
python - "${gatk_table}" "${ASE_TEMP_DIR}/${sample_id}.snv.bed" <<'PY'
import csv,sys
src,dst=sys.argv[1:]
with open(src,newline='') as fh:
    rd=csv.DictReader(fh,delimiter='\t')
    if not rd.fieldnames: raise SystemExit('ASEReadCounter table has no header')
    aliases={'contig':['contig','CONTIG'],'position':['position','POSITION'],'ref':['refAllele','REF_ALLELE','ref'],'alt':['altAllele','ALT_ALLELE','alt'],'refCount':['refCount','REF_COUNT'],'altCount':['altCount','ALT_COUNT'],'totalCount':['totalCount','TOTAL_COUNT']}
    def choose(k):
        for x in aliases[k]:
            if x in rd.fieldnames: return x
        raise SystemExit(f'Missing required column {k}; columns={rd.fieldnames}')
    c={k:choose(k) for k in aliases}
    with open(dst,'w',newline='') as out:
        wr=csv.writer(out,delimiter='\t',lineterminator='\n')
        for row in rd:
            try: pos=int(row[c['position']])
            except (ValueError,TypeError): continue
            wr.writerow([row[c['contig']],pos-1,pos,row[c['ref']],row[c['alt']],row[c['refCount']],row[c['altCount']],row[c['totalCount']]])
PY
snv_bed="${ASE_TEMP_DIR}/${sample_id}.snv.bed"; [[ -s "${snv_bed}" ]] || as_die "No valid ASE rows parsed"
if as_sample_is_rnaseq "${sample_id}"; then
  # RNASeq has no peak stage: all valid ASEReadCounter sites continue to FDR.
  cp -f "${snv_bed}" "${out}"; as_info "RNASeq: no peak intersection applied"
else
  as_require_cmds bedtools; peak_file="${PEAK_DIR}/${sample_id}_Peaks.bed"; as_require_file "${peak_file}" "peak file"
  fixed_bed="${ASE_TEMP_DIR}/${sample_id}.snv.chrom_style_fixed.bed"
  peak_chr="$(awk 'NF>0 && $1 !~ /^#/ {print $1; exit}' "${peak_file}")"
  if [[ "${peak_chr}" == chr* ]]; then
    awk -F'\t' -v OFS='\t' '{if($1 !~ /^chr/) $1="chr"$1; print}' "${snv_bed}" > "${fixed_bed}"
  else
    awk -F'\t' -v OFS='\t' '{sub(/^chr/,"",$1); print}' "${snv_bed}" > "${fixed_bed}"
  fi
  bedtools intersect -a "${fixed_bed}" -b "${peak_file}" -u > "${out}"
  rm -f "${fixed_bed}"
fi
rm -f "${snv_bed}"
rows="$(wc -l < "${out}")"
(( rows > 0 )) || as_warn "Prepared ASE table is empty for ${sample_id}"
as_info "Prepared ${rows} ASE sites: ${out}"
