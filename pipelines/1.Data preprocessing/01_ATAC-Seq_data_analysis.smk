import os
SAMPLES = [i[:-11] for i in os.listdir("/vol2/zhangshiwen/atac_public/fastq/raw_data") if i.endswith("_1.fastq.gz")  and any(i.startswith(prefix) for prefix in ["DRR", "SRR"])]
print(SAMPLES)

ruleorder: samtools_index1 > samtools_index2
ruleorder: bedtools_bamtobed > bedtools_bamtobed_last

rule all:
    input:
        #expand("/vol2/zhangshiwen/atac_public/clean/{sample}.bw", sample=SAMPLES),
        #expand("/vol2/zhangshiwen/atac_public/peaks/{sample}_peaks.narrowPeak", sample=SAMPLES),
        expand("/vol2/zhangshiwen/atac_public/clean/{sample}.bam.bai", sample=SAMPLES),
        expand("/vol2/zhangshiwen/atac_public/clean/{sample}.last.bam.bai", sample=SAMPLES)

rule trim_galore:
    input:
        fq1 = "/vol2/zhangshiwen/atac_public/fastq/raw_data/{sample}_1.fastq.gz",
        fq2 = "/vol2/zhangshiwen/atac_public/fastq/raw_data/{sample}_2.fastq.gz"
    output:
        clean_fq1 = "/vol2/zhangshiwen/atac_public/clean/{sample}_1_val_1.fq.gz",
        clean_fq2 = "/vol2/zhangshiwen/atac_public/clean/{sample}_2_val_2.fq.gz"
    conda:
        '../Envs/trim-galore.yaml'
    params:
        q = 25,
        phred = 33,
        length = 35,
        error_rate = 0.1,
        stringency = 4,
        outdir = "/vol2/zhangshiwen/atac_public/clean"
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
        fq1 = "/vol2/zhangshiwen/atac_public/clean/{sample}_1_val_1.fq.gz",
        fq2 = "/vol2/zhangshiwen/atac_public/clean/{sample}_2_val_2.fq.gz"
    output:
        sam = "/vol2/zhangshiwen/atac_public/clean/{sample}.sam",
        log = "/vol2/zhangshiwen/atac_public/clean/{sample}_bowtie2.txt"
    conda:
        '../Envs/bowtie2.yaml'
    params:
        index = "/public/home/mengzhu/reference/sheep"
    threads: 5
    shell:
        """
        bowtie2 -p {threads} --very-sensitive -X 2000 -x {params.index} -1 {input.fq1} -2 {input.fq2} -S {output.sam} 2> {output.log}
        """

rule samtools_sort:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.sam"
    output:
        bam = "/vol2/zhangshiwen/atac_public/clean/{sample}.bam"
    conda:
        '../Envs/samtools.yaml'
    threads: 5
    shell:
        """
        samtools sort -@ {threads} -O bam -o {output.bam} {input}
        """

rule samtools_index1:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.bam"
    output:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.bam.bai"
    conda:
        '../Envs/samtools.yaml'
    shell:
        """
        samtools index {input}
        """

rule bedtools_bamtobed:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.bam"
    output:
        bed = "/vol2/zhangshiwen/atac_public/clean/{sample}.raw.bed"
    conda:
        '../Envs/bedtools.yaml'
    shell:
        """
        bedtools bamtobed -i {input} > {output}
        """

rule samtools_flagstat_raw:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.bam"
    output:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.raw.stat"
    conda:
        '../Envs/samtools.yaml'
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule sambamba_markdup:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.bam"
    output:
        bam = "/vol2/zhangshiwen/atac_public/clean/{sample}.rmdup.bam"
    conda:
        '../Envs/sambamba.yaml'
    shell:
        """
        sambamba markdup --overflow-list-size 600000 --tmpdir='./' -r {input} {output.bam}
        """

rule samtools_flagstat_rmdup:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.rmdup.bam"
    output:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.rmdup.stat"
    conda:
        '../Envs/samtools.yaml'
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule samtools_filter_mitochondria:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.rmdup.bam"
    output:
        bam = "/vol2/zhangshiwen/atac_public/clean/{sample}.last.bam"
    conda:
        '../Envs/samtools.yaml'
    threads: 5
    shell:
        """
        samtools view -h -f 2 -q 30 {input} | grep -v chrM | samtools sort -O bam -@ {threads} -o - > {output.bam}
        """

rule samtools_index2:
    input:
        bam = "/vol2/zhangshiwen/atac_public/clean/{sample}.last.bam",
        index1 = "/vol2/zhangshiwen/atac_public/clean/{sample}.bam.bai"  # Depends on the output of samtools_index1
    output:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.last.bam.bai"
    conda:
        '../Envs/samtools.yaml'
    shell:
        """
        samtools index {input.bam}
        """

rule samtools_flagstat_last:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.last.bam"
    output:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.last.stat"
    conda:
        '../Envs/samtools.yaml'
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule bedtools_bamtobed_last:
    input:
        bam = "/vol2/zhangshiwen/atac_public/clean/{sample}.last.bam",
        bed = "/vol2/zhangshiwen/atac_public/clean/{sample}.raw.bed"
    output:
        bed = "/vol2/zhangshiwen/atac_public/clean/{sample}.bed"
    conda:
        '../Envs/bedtools.yaml'
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output}
        """

rule bamCoverage:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.last.bam"
    output:
        bw = "/vol2/zhangshiwen/atac_public/clean/{sample}.bw"
    conda:
        '../Envs/bedtools.yaml'
    threads: 5
    shell:
        """
        bamCoverage -b {input} -o {output.bw} --binSize 20 --normalizeUsing RPKM --smoothLength 60 --extendReads 150 --centerReads -p {threads}
        """

rule macs2_callpeak:
    input:
        "/vol2/zhangshiwen/atac_public/clean/{sample}.bed"
    output:
        peaks = "/vol2/zhangshiwen/atac_public/peaks/{sample}_peaks.narrowPeak"
    conda:
        '../Envs/macs.yaml'
    params:
        gsize = "2628104905"
    shell:
        """
        macs2 callpeak -t {input} -g {params.gsize} --nomodel --shift -100 --extsize 200 -n {wildcards.sample} --outdir ['ATAC_hypothalamus_01', 'ATAC_spleen_02', 'ATAC_hypothalamus_02', 'ATAC_lung_02', 'ATAC_liver_02', 'ATAC_spleen_01', 'ATAC_hypothalamus_02', 'ATAC_adipose_01', 'ATAC_spleen_01', 'ATAC_liver_01', 'ATAC_duodenum_01', 'ATAC_duodenum_01', 'ATAC_lung_01', 'ATAC_duodenum_02', 'ATAC_liver_01', 'ATAC_heart_01', 'ATAC_heart_02', 'ATAC_lung_02', 'ATAC_hypothalamus_01', 'ATAC_adipose_01', 'ATAC_lung_01', 'ATAC_adipose_02', 'ATAC_duodenum_02', 'ATAC_adipose_02', 'ATAC_heart_01', 'ATAC_heart_02', 'ATAC_spleen_02', 'ATAC_liver_02']

        """
