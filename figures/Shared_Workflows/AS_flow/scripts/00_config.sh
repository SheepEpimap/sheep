#!/usr/bin/env bash
# ============================================================================
# AS pipeline central configuration
# ============================================================================
# Portability policy
#   1. External/raw inputs keep the original server paths as DEFAULTS.
#   2. Every default can be replaced either by editing this file or by exporting
#      an environment variable with the same name before running a script.
#   3. All generated files are written below this extracted package by default.
#
# Example override on another server:
#   export REFERENCE_FASTA=/new/path/Ovis_aries_v2.0.fasta
#   export INPUT_BAM_DIR=/new/path/bam
#   bash scripts/01_check_environment.sh
# ============================================================================

CONFIG_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
AS_PROJECT_ROOT="${AS_PROJECT_ROOT:-$(cd "${CONFIG_DIR}/.." && pwd)}"

# ---------- External software/environment locations: replaceable ----------
CONDA_ROOT="${CONDA_ROOT:-/vol2/wulingyun/miniconda3}"
ATAC_ENV="${ATAC_ENV:-atac_tools}"
RNASEQ_ENV="${RNASEQ_ENV:-rnaseq_tools}"
WASP_ENV="${WASP_ENV:-wasp_env}"
ASE_ENV="${ASE_ENV:-ase_env}"

# ---------- External/raw inputs: original paths retained as defaults ----------
PEAK_RAW_READS_DIR="${PEAK_RAW_READS_DIR:-/vol2/00_panlab_rawdata/01_gwasdata/WGS_235raw/35_high_depth}"
INPUT_BAM_DIR="${INPUT_BAM_DIR:-/vol2/mengzhu/snakemake_sheep/clean/bam1}"
RNASEQ_RAW_READS_DIR="${RNASEQ_RAW_READS_DIR:-/vol2/mengzhu/snakemake_sheep/Raw_Reads}"

REFERENCE_FASTA="${REFERENCE_FASTA:-/vol2/wulingyun/AS/data/Ovis_aries_v2.0.fasta}"
REFERENCE_GTF="${REFERENCE_GTF:-/vol2/wulingyun/AS/data/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf}"
REFERENCE_CHROM_SIZES_SOURCE="${REFERENCE_CHROM_SIZES_SOURCE:-/vol2/wulingyun/AS/data/Ovis_aries_v2.0.chrom.sizes}"
GENOME_SIZE="${GENOME_SIZE:-2628146905}"

PHASING_VCF="${PHASING_VCF:-/vol2/wulingyun/AS/vcf/35_sub_indv.vcf.gz}"
STAR_INDEX="${STAR_INDEX:-/vol2/wulingyun/AS/index/star_ovis_aries}"
SALMON_INDEX="${SALMON_INDEX:-/vol2/wulingyun/AS/index/salmon_ovis_aries}"

WASP_PATH="${WASP_PATH:-/vol2/wulingyun/AS/WASP/WASP-master}"
SPP_RUNNER="${SPP_RUNNER:-/vol2/wulingyun/AS/phantompeakqualtools/run_spp.R}"
SPP_R_LIB="${SPP_R_LIB:-/vol2/wulingyun/R/r_spp_lib}"

# ---------- Package-local resources and outputs ----------
RESOURCE_DIR="${RESOURCE_DIR:-${AS_PROJECT_ROOT}/resources}"
LOG_DIR="${LOG_DIR:-${AS_PROJECT_ROOT}/logs}"
WORK_DIR="${WORK_DIR:-${AS_PROJECT_ROOT}/work}"

SAMPLE_MAP="${SAMPLE_MAP:-${RESOURCE_DIR}/sample_to_individual.txt}"
CHROM_SIZES="${CHROM_SIZES:-${RESOURCE_DIR}/Ovis_aries_v2.0.chrom.sizes}"
CHROM_LIST="${CHROM_LIST:-${RESOURCE_DIR}/chroms.txt}"
MAPPABILITY_RESOURCE_DIR="${MAPPABILITY_RESOURCE_DIR:-${RESOURCE_DIR}/mappability}"
PHASING_MAP_DIR="${PHASING_MAP_DIR:-${RESOURCE_DIR}/beagle_maps}"

PEAK_WORK_DIR="${PEAK_WORK_DIR:-${WORK_DIR}/peak_pipeline}"
PEAK_TRIM_DIR="${PEAK_TRIM_DIR:-${PEAK_WORK_DIR}/01_Trimmed_Reads}"
MAPPABILITY_DIR="${MAPPABILITY_DIR:-${PEAK_WORK_DIR}/02_Mappability_filtered}"
PEAK_METRICS_DIR="${PEAK_METRICS_DIR:-${PEAK_WORK_DIR}/03_Metrics}"
PEAK_DIR="${PEAK_DIR:-${PEAK_WORK_DIR}/04_Peak_called}"
PEAK_TEMP_DIR="${PEAK_TEMP_DIR:-${PEAK_WORK_DIR}/Temp}"

RNASEQ_WORK_DIR="${RNASEQ_WORK_DIR:-${WORK_DIR}/rnaseq}"
RNASEQ_TRIM_DIR="${RNASEQ_TRIM_DIR:-${RNASEQ_WORK_DIR}/01_Trimmed_Reads}"
RNASEQ_ALIGN_DIR="${RNASEQ_ALIGN_DIR:-${RNASEQ_WORK_DIR}/02_Aligned_BAM}"
RNASEQ_QUANT_DIR="${RNASEQ_QUANT_DIR:-${RNASEQ_WORK_DIR}/03_Quantification}"
RNASEQ_METRICS_DIR="${RNASEQ_METRICS_DIR:-${RNASEQ_WORK_DIR}/04_Metrics}"
RNASEQ_TEMP_DIR="${RNASEQ_TEMP_DIR:-${RNASEQ_WORK_DIR}/Temp}"

PHASING_DIR="${PHASING_DIR:-${WORK_DIR}/phasing}"
WASP_HDF5_DIR="${WASP_HDF5_DIR:-${PHASING_DIR}/hdf5}"
WASP_HAPS_H5="${WASP_HAPS_H5:-${WASP_HDF5_DIR}/haps.h5}"
WASP_SNP_TAB_H5="${WASP_SNP_TAB_H5:-${WASP_HDF5_DIR}/snp_tab.h5}"
WASP_SNP_INDEX_H5="${WASP_SNP_INDEX_H5:-${WASP_HDF5_DIR}/snp_index.h5}"
WASP_WORK_DIR="${WASP_WORK_DIR:-${WORK_DIR}/wasp}"
WASP_TEMP_DIR="${WASP_TEMP_DIR:-${WASP_WORK_DIR}/temp}"
WASP_FIND_DIR="${WASP_FIND_DIR:-${WASP_WORK_DIR}/find_intersecting_snps}"
WASP_REMAP_DIR="${WASP_REMAP_DIR:-${WASP_WORK_DIR}/remap}"
WASP_FINAL_DIR="${WASP_FINAL_DIR:-${WASP_WORK_DIR}/filter_remapped_reads}"
WASP_SNP_TEXT_DIR="${WASP_SNP_TEXT_DIR:-${WASP_WORK_DIR}/output_snp}"

ASE_WORK_DIR="${ASE_WORK_DIR:-${WORK_DIR}/ase}"
ASE_TABLE_DIR="${ASE_TABLE_DIR:-${ASE_WORK_DIR}/ase_read_counter}"
ASE_COUNT_DIR="${ASE_COUNT_DIR:-${ASE_WORK_DIR}/prepared_counts}"
ASE_TEMP_DIR="${ASE_TEMP_DIR:-${ASE_WORK_DIR}/temp}"

STAT_WORK_DIR="${STAT_WORK_DIR:-${WORK_DIR}/statistics}"
BINOMIAL_DIR="${BINOMIAL_DIR:-${STAT_WORK_DIR}/binomial_results}"
FDR_DETAIL_DIR="${FDR_DETAIL_DIR:-${BINOMIAL_DIR}/fdr_details}"

# ---------- Analysis settings ----------
MAPQ="${MAPQ:-30}"
PEAK_THREADS="${PEAK_THREADS:-12}"
RNASEQ_THREADS="${RNASEQ_THREADS:-8}"
WASP_THREADS="${WASP_THREADS:-4}"
PHASING_THREADS="${PHASING_THREADS:-15}"

# Required rule:
ATAC_QVALUE="${ATAC_QVALUE:-0.01}"
OTHER_PEAK_QVALUE="${OTHER_PEAK_QVALUE:-0.05}"

# Input BAM suffixes used by the original data layout.
NON_ATAC_BAM_SUFFIX="${NON_ATAC_BAM_SUFFIX:-bowtie2.mapped.filtered.sort.bam}"
ATAC_BAM_SUFFIX="${ATAC_BAM_SUFFIX:-last.bam}"

# Create package-local directories. No external/raw directory is created or modified.
mkdir -p \
  "${RESOURCE_DIR}" "${LOG_DIR}" \
  "${PEAK_TRIM_DIR}" "${MAPPABILITY_DIR}" "${PEAK_METRICS_DIR}" "${PEAK_DIR}" "${PEAK_TEMP_DIR}" \
  "${RNASEQ_TRIM_DIR}" "${RNASEQ_ALIGN_DIR}" "${RNASEQ_QUANT_DIR}" "${RNASEQ_METRICS_DIR}" "${RNASEQ_TEMP_DIR}" \
  "${PHASING_DIR}" "${PHASING_MAP_DIR}" "${WASP_HDF5_DIR}" \
  "${WASP_TEMP_DIR}" "${WASP_FIND_DIR}" "${WASP_REMAP_DIR}" "${WASP_FINAL_DIR}" "${WASP_SNP_TEXT_DIR}" \
  "${ASE_TABLE_DIR}" "${ASE_COUNT_DIR}" "${ASE_TEMP_DIR}" \
  "${BINOMIAL_DIR}" "${FDR_DETAIL_DIR}"
