import os
SAMPLES = [i[:-9] for i in os.listdir("/vol2/mengzhu/snakemake_sheep/H3") if i.endswith(".fq.gz") and i.startswith("H3")]
print(SAMPLES)
rule all:
    input:
         expand('/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.mapped.filtered.sort.bam.bai', sample=SAMPLES),  # Use expand to handle wildcards

rule trim_galore:
    input:
        fq1 = "/vol2/mengzhu/snakemake_sheep/Raw_Reads/{sample}_R1.fq.gz",
        fq2 = "/vol2/mengzhu/snakemake_sheep/Raw_Reads/{sample}_R2.fq.gz"
    output:
        clean_fq1 = "clean/{sample}_R1_val_1.fq.gz",
        clean_fq2 = "clean/{sample}_R2_val_2.fq.gz"
    conda:
        '../Envs/trim-galore.yaml'
    params:
        q = 25,
        phred = 33,
        length = 35,
        error_rate = 0.1,
        stringency = 4,
        outdir = "clean"
    log:
        "../logs/trim_galore_{sample}.log"
    threads: 4
    shell:
        """
        trim_galore -q {params.q} --phred{params.phred} --length {params.length} \
        -e {params.error_rate} --stringency {params.stringency} \
        --paired {input.fq1} {input.fq2} --gzip -o {params.outdir} &> {log}
        """

rule bowtie2_align:
    input:
        fq1 = "/vol2/mengzhu/snakemake_sheep/clean/{sample}_R1_val_1.fq.gz",
        fq2 = "/vol2/mengzhu/snakemake_sheep/clean/{sample}_R2_val_2.fq.gz",
    output:
        sam = "/vol2/mengzhu/snakemake_sheep/clean/{sample}_bowtie2.sam",
        log = "/vol2/mengzhu/snakemake_sheep/clean/{sample}_bowtie2.txt"
    conda:
        '../Envs/bowtie2.yaml'
    threads: 20
    shell:
        """
        bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
        -I 30 -X 700 -p {threads} \
        -x /public/home/mengzhu/reference/sheep \
        -1 {input.fq1} -2 {input.fq2} \
        -S {output.sam} 2> {output.log}
        """

# Rule to sort BAM file by coordinate
rule sort_sam:
    input:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}_bowtie2.sam"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.sorted.sam"
    conda:
        '../Envs/picard.yaml'
    shell:
        "picard SortSam "
        "INPUT={input} "
        "OUTPUT={output} "
        "SORT_ORDER=coordinate"

# Rule to mark duplicates
rule mark_duplicates:
    input:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.sorted.sam"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.sorted.dupMarked.sam"
    conda:
        '../Envs/picard.yaml'
    params:
        metrics_file="/vol2/mengzhu/snakemake_sheep/clean/{sample}.picard.dupMark.txt"
    shell:
        "picard MarkDuplicates "
        "INPUT={input} "
        "OUTPUT={output} "
        "METRICS_FILE={params.metrics_file}"

# Rule to remove duplicates
rule remove_duplicates:
    input:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.sorted.sam"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.sorted.rmdupMarked.sam"
    conda:
        '../Envs/picard.yaml'
    params:
        metrics_file="/vol2/mengzhu/snakemake_sheep/clean/{sample}.picard.rmDup.txt"
    shell:
        "picard MarkDuplicates "
        "INPUT={input} "
        "OUTPUT={output} "
        "REMOVE_DUPLICATES=true "
        "METRICS_FILE={params.metrics_file}"

# Rule to convert SAM to BAM and filter unmapped reads
rule sam_to_bam:
    input:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.sorted.rmdupMarked.sam"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.mapped.bam"
    conda:
        '../Envs/samtools.yaml'
    shell:
        "samtools view -q 2 -bS -F 0x04 {input} > {output}"

# Rule to sort BAM file
#rule sort_bam:
    #input:
       # "clean/{sample}.bowtie2.mapped.bam"
    #output:
        #"clean/{sample}.bowtie2.mapped.sort.bam"
    #conda:
        #'Envs/samtools.yaml'
    #params:
        #threads=5  # Set the number of threads
   # shell:
        #"samtools sort -@ {params.threads} -O bam -o {output} {input}"

# Rule to index BAM file
##rule index_bam:
    #input:
        #"clean/{sample}.bowtie2.mapped.sort.bam"
    #output:
       # "clean/{sample}.bowtie2.mapped.sort.bam.bai"
    #conda:
        #'Envs/samtools.yaml'
    #shell:
        #"samtools index {input}"


rule filter_alignments:
    input: 
        sam = '/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.sorted.rmdupMarked.sam'
    output: 
        bam = '/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.mapped.filtered.bam'
    #group: 'bwa_align'
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools view -h -F 1804 -q 30 {input.sam} | grep -v XA:Z | grep -v SA:Z | samtools view -S -b - > {output.bam}'

# Rule to sort BAM file
rule sort_bam:
    input:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.mapped.filtered.bam"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.mapped.filtered.sort.bam"
    conda:
        '../Envs/samtools.yaml'
    params:
        threads=5  # Set the number of threads
    shell:
        "samtools sort -@ {params.threads} -O bam -o {output} {input}"

# Rule to index BAM file
rule index_bam:
    input:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.mapped.filtered.sort.bam"
    output:
        "/vol2/mengzhu/snakemake_sheep/clean/{sample}.bowtie2.mapped.filtered.sort.bam.bai"
    conda:
        '../Envs/samtools.yaml'
    shell:
        "samtools index {input}"
