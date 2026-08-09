#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd bedtools
for f in "${GENE_BED}" "${TSS_2K_BED}" "${INTERGENIC_BED}" "${UTR5_BED}" \
         "${UTR3_BED}" "${EXON_BED}" "${INTRON_BED}" "${CDS_BED}"; do
  require_file "$f"
done

cd "${RESULT_DIR}"
out="${RESULT_DIR}/Peak_number_summary.csv"
printf 'sheep\tgenic\tTSS_proximal\tintergenic\t5UTR\t3UTR\texon\tCDS\tIntron\tall_peak\tgenic_pe\tTSS2k_pe\tintergenic_pe\tgenic_cov\tTSS2k_cov\tintergenic_cov\t5UTR_cov\t3UTR_cov\tExon_cov\tIntron_cov\tCDS_cov\n' > "${out}"

count_regions() {
  bedtools intersect -wa -u -a "$1" -b "$2" | awk 'END{print NR+0}'
}

# This function intentionally mirrors the coverage definition in the source
# Markdown: it sums the complete lengths of A records returned by intersect and
# divides by the merged annotation length. It is retained for result continuity.
coverage_ratio_source_definition() {
  local query="$1" annotation="$2"
  local numerator denominator
  numerator=$(bedtools intersect -a "${query}" -b <(sort -k1,1 -k2,2n "${annotation}" | bedtools merge -i stdin) \
    | awk '{s+=($3-$2+1)} END{print s+0}')
  denominator=$(sort -k1,1 -k2,2n "${annotation}" | bedtools merge -i stdin \
    | awk '{s+=($3-$2+1)} END{print s+0}')
  safe_ratio "${numerator}" "${denominator}" 10
}

for tissue in "${TISSUES[@]}"; do
  bed="${RESULT_DIR}/${tissue}_super_enhancer.bed"
  require_file "${bed}"

  genic=$(count_regions "${bed}" "${GENE_BED}")
  tss=$(count_regions "${bed}" "${TSS_2K_BED}")
  intergenic=$(count_regions "${bed}" "${INTERGENIC_BED}")
  utr5=$(count_regions "${bed}" "${UTR5_BED}")
  utr3=$(count_regions "${bed}" "${UTR3_BED}")
  exon=$(count_regions "${bed}" "${EXON_BED}")
  intron=$(count_regions "${bed}" "${INTRON_BED}")
  cds=$(count_regions "${bed}" "${CDS_BED}")
  total=$(awk 'END{print NR+0}' "${bed}")

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tissue}" "${genic}" "${tss}" "${intergenic}" "${utr5}" "${utr3}" "${exon}" "${cds}" "${intron}" "${total}" \
    "$(safe_ratio "${genic}" "${total}" 10)" \
    "$(safe_ratio "${tss}" "${total}" 10)" \
    "$(safe_ratio "${intergenic}" "${total}" 10)" \
    "$(coverage_ratio_source_definition "${bed}" "${GENE_BED}")" \
    "$(coverage_ratio_source_definition "${bed}" "${TSS_2K_BED}")" \
    "$(coverage_ratio_source_definition "${bed}" "${INTERGENIC_BED}")" \
    "$(coverage_ratio_source_definition "${bed}" "${UTR5_BED}")" \
    "$(coverage_ratio_source_definition "${bed}" "${UTR3_BED}")" \
    "$(coverage_ratio_source_definition "${bed}" "${EXON_BED}")" \
    "$(coverage_ratio_source_definition "${bed}" "${INTRON_BED}")" \
    "$(coverage_ratio_source_definition "${bed}" "${CDS_BED}")" >> "${out}"
done

log "Genomic annotation summary written to ${out}."
