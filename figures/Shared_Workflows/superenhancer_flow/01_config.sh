#!/usr/bin/env bash
# ==============================================================================
# Super-enhancer workflow: shared configuration
# ==============================================================================
# The workflow is anchored to the project directory requested by the user.
# All generated files and intermediate results are written below WORK_ROOT or
# directories derived from it, so the package itself can be unpacked anywhere.
#
# Paths under "External input and software paths" deliberately retain the
# original absolute locations recorded in the source Markdown. When migrating
# to a server with a different storage layout, update the variables in THIS file
# only; do not edit individual step scripts. This centralised configuration
# keeps path configuration consistent across workflow stages.
# ============================================================================

# ---- Anchored project paths -------------------------------------------------
WORK_ROOT="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer"
STATE_DIR="$(dirname "${WORK_ROOT}")"
RESULT_DIR="${WORK_ROOT}/AA_super_result_1"
STATE_ID_DIR="${STATE_DIR}/AA_state_id"
ENHANCER_COMBINE_DIR="${STATE_DIR}/AA_enhancer_combine"

# ---- External input and software paths --------------------------------------
# These defaults reproduce the source workflow. Replace only when the same
# underlying resource is stored elsewhere on the destination server.
ROSE_HOME="/vol2/mengzhu/soft/ROSE-1.3.2"
GTF_FILE="/vol2/mengzhu/genome/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.113.gtf.gz"
H3K27AC_BAM_DIR="/vol2/mengzhu/snakemake_sheep/clean"
H3K27AC_BW_DIR="/vol2/mengzhu/snakemake_sheep/clean/bw"
H3K27AC_LABELS_FILE="${H3K27AC_BW_DIR}/labels.txt"

GENOME_FASTA="/vol2/mengzhu/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna"
PROTEIN_TSS_BED="/vol2/mengzhu/genome/fold_enrich_region_protein/protein_TSS_esemble100_colin.bed"

GENE_BED="/vol2/mengzhu/genome/part_change_esemb100/Gene_esemble100_colin_ncbi.bed"
TSS_2K_BED="/vol2/mengzhu/genome/part_change_esemb100/TSS_esemble100_colin.bed_2k_new.bed"
INTERGENIC_BED="/vol2/mengzhu/genome/part_change_esemb100/Intergenic_esemble100_colin.bed"
UTR5_BED="/vol2/mengzhu/genome/part_change_esemb100/5UTR_esemble100_last_ucsc.bed"
UTR3_BED="/vol2/mengzhu/genome/part_change_esemb100/3UTR_esemble100_last_ucsc.bed"
EXON_BED="/vol2/mengzhu/genome/part_change_esemb100/exon_esemble100_colin.bed"
INTRON_BED="/vol2/mengzhu/genome/part_change_esemb100/Intron_esemble100_colin.bed"
CDS_BED="/vol2/mengzhu/genome/part_change_esemb100/CDS_esemble100_colin.bed"

EXPRESSION_TPM="/vol2/mengzhu/snakemake_sheep/expressiondir/all_tisssues_expression_tpm.csv"
EXPRESSION_TAU="/vol2/mengzhu/snakemake_sheep/expressiondir/all_tisssues_expression_tpm.median.tau.csv"
EXPRESSION_AVERAGE="/vol2/mengzhu/snakemake_sheep/expressiondir/Average/all_tisssues_expression_aveage_tpm.csv"
SHEEP_HUMAN_CONSERVATION="/vol2/mengzhu/genome/conservation_1/sheep_to_human_conservation_last.txt"
SHEEP_HUMAN_ID_MAP="/vol2/mengzhu/genome/conservation_1/sheep_human_esemble_ID.txt"

# ---- Workflow parameters ----------------------------------------------------
ROSE_GENOME_KEY="OVIARI2"
ROSE_ANNOTATION_NAME="oviAri2_refseq.ucsc"
ROSE_STITCH_DISTANCE=12500
ROSE_TSS_EXCLUSION=2500
GENOME_SIZE_BP=2478444698
THREADS="${THREADS:-24}"
HOMER_THREADS="${HOMER_THREADS:-50}"

# The expression table used in the source workflow contains metadata in columns
# 1-6 and 86 expression columns in columns 7-92.
EXPRESSION_FIRST_COLUMN=7
EXPRESSION_LAST_COLUMN=92

TISSUES=(
  abomasum adipose bone-marrow brainstem cecum cerebellum cerebral-cortex
  cervix colon cornua-uteri corpus-uteri duodenum epididymis heart hippocampus
  hypothalamus ileum jejunum kidney liver lung lymph-node mammary-gland
  medulla-oblongata midbrain muscle omasum optic-chiasm ovary oviduct pineal
  pituitary pons rectum reticulum rumen skin soft-horn spleen splenium testis
  thymus thyroid
)
STATES=(E6 E7 E8 E9 E10)
REPLICATES=(39 40)

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

die() {
  log "ERROR: $*"
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

require_file() {
  [[ -s "$1" ]] || die "Required file is missing or empty: $1"
}

safe_ratio() {
  # Usage: safe_ratio numerator denominator [scale]
  awk -v n="$1" -v d="$2" -v s="${3:-6}" 'BEGIN {
    if (d == 0 || d == "") print "NA";
    else printf "%.*f", s, n/d;
  }'
}
