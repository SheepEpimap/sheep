#Normalize RNA-Seq data
import os

SAMPLES = [
    i[:-4] for i in os.listdir("Aligned_Reads")
    if i.endswith(".bam") and i.startswith("RNA")
]

print(SAMPLES)

GENOME_SIZE = "2628104905"
CHROMSIZES = "/vol2/mengzhu/genome/reference/sheep.size"  # Replace with the actual chrom.sizes file path

rule all:
    input:
        # Original RPKM bigWig
        expand(
            "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}.bw",
            sample=SAMPLES
        ),
        # ZScore bigwig
        expand(
            "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}_ZScores.bw",
            sample=SAMPLES
        )

# 1) BAM -> RPKM bigwig
rule rna_bam_coverage_bigwig:
    input:
        bam = "Aligned_Reads/{sample}.bam",
        bai = "Aligned_Reads/{sample}.bam.bai"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}.bw"
    conda:
        "/vol2/mengzhu/snakemake_sheep/Envs/deeptools.yaml"
    threads: 25
    params:
        gsize = GENOME_SIZE
    shell:
        "bamCoverage -p {threads} -b {input.bam} -o {output} "
        "--normalizeUsing RPKM --effectiveGenomeSize {params.gsize} -bs 100"

# 2) BAM -> RPKM bedGraph
rule rna_bam_coverage_bedgraph:
    input:
        bam = "Aligned_Reads/{sample}.bam",
        bai = "Aligned_Reads/{sample}.bam.bai"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}.bdg"
    conda:
        "/vol2/mengzhu/snakemake_sheep/Envs/deeptools.yaml"
    threads: 25
    params:
        gsize = GENOME_SIZE
    shell:
        "bamCoverage -p {threads} -b {input.bam} -o {output} -of bedgraph "
        "--normalizeUsing RPKM --effectiveGenomeSize {params.gsize} -bs 100"

# 3) Apply Z-score normalization to bedGraph
rule zscore_normalize_bedgraph:
    input:
        bedgraph   = "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}.bdg",
        chromsizes = CHROMSIZES    # Pass the chrom.sizes file to the script here
    output:
        bedgraph = "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}_ZScores.bdg"
    threads: 4
    conda:
        "/vol2/mengzhu/snakemake_sheep/Envs/scipy.yaml"
    script:
        "/vol2/mengzhu/snakemake_sheep/Scripts/ZScore_Normalize_BedGraph.py"

# 4) ZScore bedGraph -> ZScore bigwig
rule zscore_bedgraph_to_bigwig:
    input:
        bedgraph   = "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}_ZScores.bdg",
        chromsizes = CHROMSIZES
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/bw/{sample}_ZScores.bw"
    conda:
        "/vol2/mengzhu/snakemake_sheep/Envs/bdg2bw.yaml"
    shell:
        "bedGraphToBigWig {input.bedgraph} {input.chromsizes} {output}"
