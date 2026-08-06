import os

# --------------------
# Shared configuration
# --------------------
UNBLACKLIST_DIR = "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist"
PEAK_DIR        = "/vol2/mengzhu/snakemake_sheep/peak_unblacklist"
GENOME_SIZE     = "2628104905"

# Read the directory only once to avoid repeated I/O
bam_files = [f for f in os.listdir(UNBLACKLIST_DIR) if f.endswith(".bam")]

# Classify samples by prefix
SAMPLES_ATAC   = [f[:-4] for f in bam_files if f.startswith("ATAC")]
SAMPLES_BROAD  = [f[:-4] for f in bam_files if f.startswith(("H3K4me1", "H3K27me3"))]
SAMPLES_NARROW = [f[:-4] for f in bam_files if f.startswith(("H3K4me3", "H3K27ac"))]

print("ATAC:", SAMPLES_ATAC)
print("BROAD (H3K4me1/H3K27me3):", SAMPLES_BROAD)
print("NARROW (H3K4me3/H3K27ac):", SAMPLES_NARROW)


# --------------------
# Files required in the final output
# --------------------
rule all:
    input:
        # broadPeak: H3K4me1 / H3K27me3
        expand(f"{PEAK_DIR}" + "/{sample}_peaks.broadPeak",  sample=SAMPLES_BROAD),
        # narrowPeak: H3K4me3 / H3K27ac
        expand(f"{PEAK_DIR}" + "/{sample}_peaks.narrowPeak", sample=SAMPLES_NARROW),
        # narrowPeak: ATAC
        expand(f"{PEAK_DIR}" + "/{sample}_peaks.narrowPeak", sample=SAMPLES_ATAC)


# --------------------
# Rule 1: ATAC
# --------------------
rule macs2_callpeak_atac:
    # Allow only sample names starting with ATAC
    wildcard_constraints:
        sample = r"ATAC.*"
    input:
        bam = f"{UNBLACKLIST_DIR}" + "/{sample}.bam"
    output:
        peaks = f"{PEAK_DIR}" + "/{sample}_peaks.narrowPeak"
    conda:
        "/vol2/mengzhu/snakemake_sheep/Envs/macs.yaml"
    params:
        gsize = GENOME_SIZE
    shell:
        """
        macs2 callpeak \
            -t {input.bam} \
            -n {wildcards.sample} \
            -g {params.gsize} \
            -q 0.01 \
            -f BAMPE \
            --nomodel \
            --shift -100 \
            --extsize 200 \
            -B --keep-dup all --SPMR \
            --outdir {PEAK_DIR}
        """


# --------------------
# Rule 3: narrow peaks for H3K4me3/H3K27ac
# --------------------
rule macs2_callpeak_narrow:
    # Allow only sample names starting with H3K4me3 or H3K27ac
    wildcard_constraints:
        sample = r"(H3K4me3|H3K27ac).*"
    input:
        bam = f"{UNBLACKLIST_DIR}" + "/{sample}.bam"
    output:
        peaks = f"{PEAK_DIR}" + "/{sample}_peaks.narrowPeak"
    conda:
        "/vol2/mengzhu/snakemake_sheep/Envs/macs.yaml"
    params:
        genome_size = GENOME_SIZE
    shell:
        """
        macs2 callpeak \
            -t {input.bam} \
            -n {wildcards.sample} \
            -g {params.genome_size} \
            -q 0.01 \
            -f BAMPE \
            --fix-bimodal \
            --extsize 200 \
            -B --keep-dup all --SPMR \
            --to-large \
            --outdir {PEAK_DIR}
        """


# --------------------
# Rule 2: broad peaks for H3K4me1/H3K27me3
# --------------------
rule macs2_callpeak_broad:
    input:
        bam = f"{UNBLACKLIST_DIR}" + "/{sample}.bam"
    output:
        peaks = f"{PEAK_DIR}" + "/{sample}_peaks.broadPeak"
    conda:
        "/vol2/mengzhu/snakemake_sheep/Envs/macs.yaml"
    params:
        genome_size = GENOME_SIZE
    shell:
        """
        macs2 callpeak \
            -t {input.bam} \
            -n {wildcards.sample} \
            -g {params.genome_size} \
            -q 0.05 \
            -f BAMPE \
            --fix-bimodal \
            --extsize 200 \
            -B --keep-dup all --SPMR \
            --broad --to-large \
            --outdir {PEAK_DIR}
        """

