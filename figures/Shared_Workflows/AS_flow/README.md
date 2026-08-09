# Portable Allele-Specific Analysis Pipeline

A SLURM-oriented workflow for allele-specific analysis of ATAC-seq, histone-mark ChIP-seq, and RNA-seq data. The pipeline performs assay-specific preprocessing, read-mapping bias correction with WASP, allele counting with GATK ASEReadCounter, and multiple-testing correction.

## Key features

- Centralized, replaceable path and parameter configuration.
- Package-local output, temporary, and log directories.
- Separate preprocessing branches for peak-based assays and RNA-seq.
- Assay-aware MACS2 thresholds:
  - ATAC-seq: `q = 0.01`
  - Other peak-based assays: `q = 0.05`
  - RNA-seq: no peak calling
- WASP remapping-bias correction with merged retained and remapped reads.
- GATK ASEReadCounter-based allele counting.
- Fast analytical or optional simulation-based false discovery rate analysis.
- Static validation and software-version capture for reproducibility.

## Workflow overview

```text
Peak-based assays
existing BAM
  -> GenMap mappability filtering
  -> SPP fragment-length estimation
  -> MACS2 broad-peak calling
  -> WASP bias correction
  -> GATK ASEReadCounter
  -> retain ASE sites overlapping assay-specific peaks
  -> statistical testing and FDR correction

RNA-seq
paired-end FASTQ
  -> Trim Galore
  -> STAR alignment and MAPQ filtering
  -> StringTie and Salmon quantification
  -> WASP bias correction using STAR remapping
  -> GATK ASEReadCounter
  -> retain all valid ASEReadCounter sites
  -> statistical testing and FDR correction
```

## Repository structure

```text
.
├── README.md
├── CODE_AVAILABILITY.md
├── CONTRIBUTING.md
├── MANIFEST.sha256
├── envs/                  Conda environment definitions
├── resources/             Sample metadata and generated small reference files
├── scripts/               Numbered workflow scripts and shared configuration
├── work/                  Generated analysis outputs; created automatically
└── logs/                  SLURM, runtime, and reproducibility logs
```

Large sequencing files, reference genomes, indices, phased VCF files, and external software repositories are not bundled. Their locations are configured in `scripts/00_config.sh`.

## System requirements

The workflow is designed for a Linux high-performance computing environment with:

- Bash
- SLURM workload manager
- Conda or Mamba
- Sufficient memory and storage for genome-scale BAM, VCF, and FASTQ processing

The required bioinformatics tools are defined in:

```text
envs/atac_tools.yml
envs/rnaseq_tools.yml
envs/wasp_env.yml
envs/ase_env.yml
```

Several scripts contain default `#SBATCH` resource requests. Review the partition, CPU, memory, and wall-time directives before running the workflow on a different cluster.

## Installation

### 1. Extract the release

```bash
unzip AS_pipeline_portable_20260803.zip
cd AS_pipeline_portable_20260803
```

### 2. Run static validation

```bash
bash scripts/22_self_test.sh
```

This check validates Bash syntax, line endings, the assay-specific peak thresholds, and the absence of external absolute paths outside the central configuration file.

### 3. Create or update Conda environments

```bash
bash envs/01_create_environments.sh
```

Mamba can be used instead of Conda by creating the environments directly from the YAML files.

## Configuration

All replaceable external paths and major analysis parameters are defined in:

```text
scripts/00_config.sh
```

The default values may be changed by editing that file or by exporting environment variables before execution. Environment-variable overrides take precedence.

Example:

```bash
export CONDA_ROOT=/opt/miniconda3
export REFERENCE_FASTA=/data/reference/Ovis_aries_v2.0.fasta
export REFERENCE_GTF=/data/reference/Ovis_aries_v2.0.gtf
export PHASING_VCF=/data/genotypes/cohort.vcf.gz
export INPUT_BAM_DIR=/data/epigenome/bam
export RNASEQ_RAW_READS_DIR=/data/rnaseq/fastq
export STAR_INDEX=/data/index/star
export SALMON_INDEX=/data/index/salmon
export WASP_PATH=/opt/WASP

bash scripts/01_check_environment.sh
```

Generated files are written below `work/`, `resources/`, and `logs/` unless their corresponding variables are overridden.

## Required inputs

### Reference resources

Configure the following variables in `scripts/00_config.sh`:

| Variable | Required content |
|---|---|
| `REFERENCE_FASTA` | Reference genome FASTA |
| `REFERENCE_GTF` | Gene annotation GTF for RNA-seq |
| `GENOME_SIZE` | Effective genome size used by MACS2 |
| `PHASING_VCF` | Cohort genotype VCF compressed with bgzip |
| `STAR_INDEX` | STAR genome index |
| `SALMON_INDEX` | Salmon transcriptome index |
| `WASP_PATH` | WASP source directory |
| `SPP_RUNNER` | `run_spp.R` from phantompeakqualtools |

The reference FASTA must use the same chromosome identifiers as the VCF, BAM files, and chromosome-size files.

### Sample-to-individual map

Edit `resources/sample_to_individual.txt` before batch execution. The file must contain two tab- or whitespace-delimited columns without a header:

```text
ATAC_thymus_39      individual_001
RNASeq_liver_39     individual_001
H3K27ac_liver_39    individual_001
```

Column 1 is the sample identifier. Column 2 is the individual identifier represented in the phased genotype data.

### Peak-assay BAM naming

The default batch discovery rules are:

```text
ATAC*.last.bam
*.bowtie2.mapped.filtered.sort.bam
```

Files ending in `raw_sorted.bam` are excluded. These suffixes can be changed with `ATAC_BAM_SUFFIX` and `NON_ATAC_BAM_SUFFIX`.

### RNA-seq FASTQ naming

RNA-seq samples must use paired-end files named:

```text
RNASeq_<sample>_R1.fq.gz
RNASeq_<sample>_R2.fq.gz
```

The sample identifier is the filename prefix preceding `_R1.fq.gz` or `_R2.fq.gz`.

## Execution

### Step 1: Validate paths and software

```bash
bash scripts/01_check_environment.sh
```

Files reported as pending may be outputs of later preparation steps. Missing external inputs or executables must be resolved before the dependent step is submitted.

### Step 2: Prepare chromosome and mappability resources

```bash
sbatch scripts/02_prepare_reference_mappability.sh
```

This step creates chromosome-size files and GenMap resources for `k = 36`, `40`, and `72`.

### Step 3A: Process peak-based assays

Batch mode:

```bash
bash scripts/04_batch_submit_peak_pipeline.sh
```

Single-sample mode:

```bash
sbatch scripts/03_process_mappability_and_peaks.sh \
  ATAC_thymus_39 last.bam

sbatch scripts/03_process_mappability_and_peaks.sh \
  H3K27ac_liver_39 bowtie2.mapped.filtered.sort.bam
```

The script rejects sample identifiers beginning with `RNASeq` to prevent accidental peak calling on RNA-seq data.

### Step 3B: Process RNA-seq

Batch mode:

```bash
bash scripts/06_batch_submit_rnaseq.sh
```

Single-sample mode:

```bash
sbatch scripts/05_process_rnaseq.sh RNASeq_liver_39
```

RNA-seq processing includes read trimming, STAR alignment, MAPQ filtering, StringTie assembly and abundance estimation, and Salmon quantification. It does not run GenMap filtering, SPP, MACS2, or peak overlap.

### Step 4: Prepare and run genotype phasing

```bash
sbatch scripts/07_prepare_phasing_inputs.sh
sbatch scripts/08_run_phasing.sh
```

Review the generated Beagle map files in `resources/beagle_maps/`. A validated species-specific recombination map should be used when available.

### Step 5: Build WASP HDF5 resources

Compile the WASP `snp2h5` executable if it is not already available, then run:

```bash
sbatch scripts/09_prepare_wasp_hdf5.sh
```

Expected outputs:

```text
work/phasing/hdf5/haps.h5
work/phasing/hdf5/snp_tab.h5
work/phasing/hdf5/snp_index.h5
```

### Step 6: Run WASP mapping-bias correction

Batch mode:

```bash
bash scripts/11_batch_submit_wasp.sh
```

Single-sample mode:

```bash
sbatch scripts/10_run_wasp.sh ATAC_thymus_39
sbatch scripts/10_run_wasp.sh RNASeq_liver_39
```

An explicit input BAM can be supplied as the second argument:

```bash
sbatch scripts/10_run_wasp.sh SAMPLE_ID /path/to/input.bam
```

For RNA-seq, remapping is performed with STAR. For peak-based assays, remapping is performed with BWA. The final BAM combines reads that did not require remapping with reads that passed the WASP remapping test.

### Step 7: Prepare phased VCF files for GATK

```bash
sbatch scripts/12_prepare_phased_vcfs.sh
```

### Step 8: Count alleles with GATK ASEReadCounter

Batch mode:

```bash
bash scripts/14_batch_submit_gatk_ase.sh
```

Single-sample mode:

```bash
sbatch scripts/13_run_gatk_ase.sh SAMPLE_ID
```

### Step 9: Prepare assay-specific ASE count tables

Batch mode:

```bash
bash scripts/16_batch_prepare_ase_counts.sh
```

Single-sample mode:

```bash
sbatch scripts/15_prepare_ase_counts.sh SAMPLE_ID
```

Peak-based assays retain ASE sites overlapping the corresponding peak file. RNA-seq retains all valid ASEReadCounter sites because no peak set is generated for this assay.

### Step 10: Run statistical testing and FDR correction

The recommended analytical mode is:

```bash
bash scripts/20_batch_submit_fdr.sh resources/sample_to_individual.txt adjust
```

Single-sample execution:

```bash
sbatch scripts/19_run_fdr.sh SAMPLE_ID adjust
```

An optional simulation mode is available:

```bash
sbatch scripts/19_run_fdr.sh SAMPLE_ID simulation
```

The simulation count and number of cores can be overridden:

```bash
export FDR_SIMULATIONS=10000
export FDR_CORES=8
sbatch scripts/19_run_fdr.sh SAMPLE_ID simulation
```

### Step 11: Summarize retained sites

```bash
bash scripts/21_count_snps.sh
```

### Step 12: Record the software environment

Run this step on the compute environment used for the final analysis:

```bash
bash scripts/23_capture_software_versions.sh
```

The report is written to `logs/reproducibility/` and records the operating system, shell, Git commit when available, Conda environments, and versions of major external tools.

## Assay-specific parameters

The following defaults are enforced in `scripts/00_config.sh`:

```bash
ATAC_QVALUE=0.01
OTHER_PEAK_QVALUE=0.05
MAPQ=30
```

`03_process_mappability_and_peaks.sh` selects the MACS2 threshold from the sample identifier:

| Assay | Peak calling | Default threshold |
|---|---:|---:|
| ATAC-seq | Yes | `q = 0.01` |
| Histone-mark or other peak-based assay | Yes | `q = 0.05` |
| RNA-seq | No | Not applicable |

Any parameter changed for a manuscript analysis should be recorded in the Methods section and in the archived configuration used for that analysis.

## Main outputs

| Stage | Output directory | Principal files |
|---|---|---|
| Peak preprocessing | `work/peak_pipeline/` | mappability-filtered BAM, SPP metrics, peak BED files |
| RNA-seq preprocessing | `work/rnaseq/` | aligned and filtered BAM, StringTie tables, Salmon directories |
| Phasing | `work/phasing/` | per-chromosome phased VCF files and WASP HDF5 resources |
| WASP | `work/wasp/filter_remapped_reads/` | `<sample>.keep.merged.sorted.bam` and index |
| ASE counting | `work/ase/ase_read_counter/` | `<sample>.output.table` |
| ASE preparation | `work/ase/prepared_counts/` | `<sample>.ase_counts.tsv` |
| Statistical analysis | `work/statistics/` | binomial-test and FDR result tables |
| Logs | `logs/` | SLURM logs, runtime logs, and software-version reports |

## Restart behavior

Most computational steps check for non-empty final outputs and skip completed work. Before restarting a failed job:

1. Review the corresponding SLURM standard-output and standard-error files.
2. Remove only incomplete or corrupted output files from the failed step.
3. Preserve validated upstream outputs.
4. Resubmit the same numbered script.

Do not treat the existence of an empty file as evidence of successful completion.

## Validation

Run static validation after modifying scripts or configuration:

```bash
bash scripts/22_self_test.sh
```

This validation does not replace a biological end-to-end test. Before production analysis, run a representative sample through the complete workflow and verify:

- chromosome naming consistency;
- expected read counts before and after filtering;
- successful WASP read merging and BAM indexing;
- non-empty ASEReadCounter output;
- expected peak-overlap behavior for peak assays;
- absence of peak filtering for RNA-seq;
- stable results when the same release and configuration are rerun.

## Reproducibility and release archiving

For a manuscript-associated release:

1. Commit the exact scripts and configuration used for the reported analysis.
2. Record the Git commit hash and create an immutable release tag.
3. Run `scripts/23_capture_software_versions.sh` on the analysis system.
4. Preserve the resolved Conda environments or explicit package specifications.
5. Include all non-default parameters required to reproduce the reported results.
6. Archive the tagged release in a DOI-minting repository such as Zenodo.
7. Cite the archived DOI and exact release in the manuscript Code Availability section.
8. Provide editors and reviewers access to the code and required documentation during peer review.

A manuscript-ready statement template is provided in [`CODE_AVAILABILITY.md`](CODE_AVAILABILITY.md).

## Citation

Cite the immutable archived release used for the analysis rather than an unversioned repository branch. Before public release, add a valid `CITATION.cff` containing the authors, software title, version, repository URL, and archival DOI.

## License

No software license is assigned by this package. Before public distribution, the copyright holder should add an appropriate license and document any restrictions affecting code reuse. An Open Source Initiative-approved license is recommended when unrestricted reuse is intended.

## Contributing

Repository modification and validation requirements are described in [`CONTRIBUTING.md`](CONTRIBUTING.md).
