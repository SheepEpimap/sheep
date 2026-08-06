# Figure 2c - DNA methylation across chromatin states

This folder contains the complete RRBS/WGBS processing, validation, state-level
summarization and plotting workflow used for Figure 2c.

## Execution order

1. `01_download_ena_fastq.sh`
2. `02_verify_fastq_md5.sh`
3. `03_run_fastqc_multiqc.sh`
4. `04_prepare_bismark_references.sh`
5. `05_make_samplesheet.sh`
6. `06_run_bismark_one_sample.sh`
7. `07_run_bismark_array.sbatch`
8. `08_submit_bismark_array.sh`
9. `09_summarize_state_methylation.sh`
10. `10_plot_state_methylation.R`
11. `11_validate_methylation_outputs.sh`

The summarization step matches each methylation sample to the ChromHMM
annotation from the same tissue, intersects Bismark CpG coverage with E1-E11,
and calculates:

`weighted methylation (%) = 100 * methylated reads / total reads`.

The plotting step retains all eleven states and does not discard highly
methylated observations. The validation step checks per-sample completion,
coverage files, exactly one value per sample-state pair, complete E1-E11
representation, and a valid 0-100% range.

## Environment

```bash
mamba env create -f environment.yml
conda activate sheep_epimap_methylation
```

Copy the scripts to:

`/vol2/zhangshiwen/rrbs/wgbs_bismark/scripts/`

Then make shell scripts executable and submit the array:

```bash
chmod +x ./*.sh
./08_submit_bismark_array.sh
```

After all samples finish:

```bash
./09_summarize_state_methylation.sh
Rscript 10_plot_state_methylation.R
./11_validate_methylation_outputs.sh
```

Absolute paths preserve the original analysis layout and can be overridden by
the environment variables defined at the beginning of each script.
