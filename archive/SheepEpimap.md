# 1.Data preprocessing

Data preprocessing

## ATAC-Seq data analysis

```shell
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
```

## CUT&Tag data analysis

```shell
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
```



## **call peak**

```shell
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

```



## Data Quality control

## SPP

```R
# https://github.com/kundajelab/phantompeakqualtools?tab=readme-ov-file
#Install
conda install bioconda::r-spp

vim spp_cut.sh
#!/bin/bash
cat H3K4me1_39.txt | while read id; do
  echo $id
  Rscript /vol2/mengzhu/soft/phantompeakqualtools/run_spp.R \
    -c=$id \
    -savp \
    -out=$(basename $id ".bowtie2.mapped.filtered.sort.bam")_spp.txt
done

vim spp_atac.sh
#!/bin/bash
cat sample1.txt | while read id; do
  echo $id
  Rscript /vol2/mengzhu/soft/phantompeakqualtools/run_spp.R \
    -c=$id \
    -savp \
    -out=$(basename $id ".last.bam")_spp.txt
done

```

## FRIP

```sh
#!/bin/bash
rm narrow_frip.txt
ls *_peaks.narrowPeak | cut -d '_' -f 1-3 | while read id; do
    echo "$id"
    bed="/vol2/mengzhu/snakemake_sheep/clean/bam1/${id}.bed"
    peak="/vol2/mengzhu/snakemake_sheep/peak3/${id}_peaks.narrowPeak"
    Reads=$(bedtools intersect -a "$bed" -b "${peak}" | wc -l | awk '{print $1}')
    totalReads=$(wc -l "$bed" | awk '{print $1}')
    totalpeaks=$(wc -l "${peak}" | awk '{print $1}')
    FRIP=$(bc <<< "scale=2; 100 * $Reads / $totalReads")
    echo -e "$id\t$totalReads\t$totalpeaks\t$FRIP" >> "narrow_frip.txt"
done

vim broad_frip.sh
#!/bin/bash
ls *_peaks.broadPeak | cut -d '_' -f 1-3 | while read id; do
    echo "$id"
    bed="/vol2/mengzhu/snakemake_sheep/clean/bam1/${id}.bed"
    peak="/vol2/mengzhu/snakemake_sheep/peak3/${id}_peaks.broadPeak"
    Reads=$(bedtools intersect -a "$bed" -b "${peak}" | wc -l | awk '{print $1}')
    totalReads=$(wc -l "$bed" | awk '{print $1}')
    totalpeaks=$(wc -l "${peak}" | awk '{print $1}')
    FRIP=$(bc <<< "scale=2; 100 * $Reads / $totalReads")
    echo -e "$id\t$totalReads\t$totalpeaks\t$FRIP" >> "broad_frip.txt"
done
sbatch -c 10  -p low --mem 10G broad_frip.sh

#Calculate genome coverage
#!/bin/bash
ls ATAC*_peaks.narrowPeak | while read id;
do
  echo ${id}
  awk -F '\t' '!/^#/ {sum += $4} END {print "'${id}'", sum / 2628104905}' ${id} >> coverage.txt
done

#!/bin/bash

# Total genome length
genome=2628104905
ls ATAC*_peaks.narrowPeak | while read id; do
  sum=$(awk '{ sum += ($3 - $2) } END { print sum }' "${id}")
  # Calculate coverage
  coverage=$(awk -v s="${sum}" -v g="${genome}" 'END { printf("%.6f", s/g) }' <<< "")
  echo "${id}  ${sum}  ${coverage}" >> coverage.txt
done
```

## Heatmap Matrix(all score)

```shell
#Normalize CUT&Tag and ATAC-Seq data
import os
SAMPLES = [i[:-9] for i in os.listdir("/vol2/mengzhu/soft/ChromImpute/sheep2/imputed_bedgraph") if i.endswith(".bedGraph") and any(i.startswith(prefix) for prefix in ["H3", "ATAC"])]
print(SAMPLES)
rule all:
    input:
         expand('/vol2/mengzhu/soft/ChromImpute/sheep2/imputed_bedgraph/{sample}_ZScores.bw', sample=SAMPLES)  # Use expand to handle wildcards

#rule bam_coverage_bedgraph:
#    input:
#        bam = 'Aligned_Reads/{SAMPLES}.bam',
#        bai = 'Aligned_Reads/{SAMPLES}.bam.bai',
#    output: 'DeepTools/{SAMPLES}.bedGraph'
#    conda:
#        'Envs/deeptools.yaml'
#    threads: 2
#    shell: 'bamCoverage -b {input.bam} -o {output} -of=bedgraph -p={threads} --normalizeUsing RPKM --effectiveGenomeSize {config[genomesize]} --ignoreDuplicates --extendReads=200 -bs 100'


rule zscore_normalize_bedgraph:
    input:
        bedgraph = '/vol2/mengzhu/soft/ChromImpute/sheep2/imputed_bedgraph/{SAMPLES}.bedGraph',
        chromsizes = '/vol2/mengzhu/genome/reference/sheep1.size'
    output:
        bedgraph = '/vol2/mengzhu/soft/ChromImpute/sheep2/imputed_bedgraph/{SAMPLES}_ZScores.bdg'
    threads: 4
    conda:
        '../Envs/scipy.yaml'
    script:
        '/vol2/mengzhu/snakemake_sheep/Scripts/ZScore_Normalize_BedGraph.py'

rule zscore_bedgraph_to_bigwig:
    input:
        bedgraph = '/vol2/mengzhu/soft/ChromImpute/sheep2/imputed_bedgraph/{SAMPLES}_ZScores.bdg',
        chromsizes = '/vol2/mengzhu/genome/reference/sheep1.size'
    output:
        '/vol2/mengzhu/soft/ChromImpute/sheep2/imputed_bedgraph/{SAMPLES}_ZScores.bw'
    conda:
        '../Envs/bdg2bw.yaml'
    shell:
        'bedGraphToBigWig {input} {output}'
```



```shell
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
```

```shell
#!/bin/bash
cd /vol2/mengzhu/snakemake_sheep/clean/bw
rm  labels.txt
ls ATAC*ZScores.bw H3*ZScores.bw RNASeq*ZScores.bw | while read id;
do
   echo $(basename $id "_ZScores.bw") >> labels.txt
done

#!/bin/bash
# Read the label file and combine its contents into one line
ma=$(cat labels.txt | perl -p -e 's/\n/ /g')

# Check the number of matching BigWig files
bw_files=$(ls ATAC_*ZScores.bw H3K27ac_*ZScores.bw H3K27me3*ZScores.bw H3K4me1*ZScores.bw H3K4me3*ZScores.bw RNASeq*ZScores.bw)
num_bw_files=$(echo "$bw_files" | wc -w)
num_labels=$(echo "$ma" | wc -w)

echo "Number of BigWig files: $num_bw_files"
echo "Number of labels: $num_labels"

# Check whether the numbers of labels and files match
if [ "$num_bw_files" -ne "$num_labels" ]; then
    echo "Error: The number of labels does not match the number of BigWig files."
    exit 1
fi

# Print matching files and labels for debugging
echo "BigWig files: $bw_files"
echo "Labels: $ma"

# Run multiBigwigSummary
multiBigwigSummary bins -p=110 \
    -b $bw_files \
    -out all_MultiBigwigSummary.npz \
    --labels $ma \
    --binSize=1000 \
    --outRawCounts all_MultiBigwigSummary.txt

plotCorrelation -in /vol2/mengzhu/snakemake_sheep/clean/bw/all_MultiBigwigSummary.npz \
--corMethod pearson \
--skipZeros --plotTitle "Pearson Correlation of all Read Counts" \
--whatToPlot heatmap \
-o /vol2/mengzhu/snakemake_sheep/clean/bw/all_Pearson_Correlation.PDF \
--removeOutliers \
--colorMap bwr \
--outFileCorMatrix /vol2/mengzhu/snakemake_sheep/clean/bw/pearsonCorr_readCounts.tab

sed "s/'//g" pearsonCorr_readCounts.tab > 1.txt
cut -f 1 1.txt | cut -d "_" -f 1  > 2.txt
cut -f 1 1.txt | cut -d "_" -f 2 > 3.txt
cut -f 1 1.txt | cut -d "_" -f 3 > 4.txt
paste 1.txt 2.txt 3.txt 4.txt > pearsonCorr_readCounts_last.csv
sed -i '1d' pearsonCorr_readCounts_last.csv
cut -f 1 pearsonCorr_readCounts.tab | cut -d "_" -f 1 | sed "s/'//g"    #Remove single quotes

#!/bin/bash
cd /vol2/mengzhu/snakemake_sheep/clean/bw
sed '1d' pearsonCorr_readCounts.tab | sed "s/'//g" > 1.txt
cut -f 1 1.txt | cut -d "_" -f 1  > 2.txt
cut -f 1 1.txt | cut -d "_" -f 2 > 3.txt
cut -f 1 1.txt | cut -d "_" -f 3 > 4.txt

cd /vol2/mengzhu/snakemake_sheep/clean/bw
cp 3.txt 3_1.txt
cat /vol2/mengzhu/snakemake_sheep/Figures/tissue_list_sheep.txt |sed 's/\t/_/g' > 6.txt
paste /vol2/mengzhu/snakemake_sheep/Figures/tissue_list_sheep.txt 6.txt >tissue_list_1.txt 
sed 's/\r/\t/g' tissue_list_1.txt > tissue_list_1_modified.txt
cat tissue_list_1_modified.txt | while read line
do 
echo $line
arr=($line)
tissue=${arr[0]}
layer1=${arr[2]}
layer2=${arr[3]}
echo $tissue
sed -i "s/\<$tissue\>/$layer2/g" 3_1.txt 
done
cat 3_1.txt |sed 's/_/\t/g' > 3_2.txt
sed '1c Tissues Layer Bigtype' 3_2.txt |  sed 's/ /\t/g' > 3_3.txt
paste 1.txt 2.txt 3_3.txt 4.txt > PearsonCorr_readCounts_last.csv


```

## Singal around tss

Plot **Supplementary Fig. 4**

```shell
vim tss.sh
#!/bin/bash
cd /vol2/mengzhu/snakemake_sheep
for i in protein_coding lncRNA pseudogene miRNA snRNA snoRNA 
do
echo $i
grep $i /vol2/mengzhu/genome/Gene_esemble100_colin.bed > 2.bed
    for tissue in cerebellum_39;
    do
    computeMatrix scale-regions \
    -R 2.bed \
    -S "DeepTools/H3K27ac_"$tissue"_ZScores.bw" "DeepTools/H3K27me3_"$tissue"_ZScores.bw" "DeepTools/H3K4me1_"$tissue"_ZScores.bw" "DeepTools/H3K4me3_"$tissue"_ZScores.bw" "DeepTools/ATAC_"$tissue"_ZScores.bw" \
    -a 2500 -b 2500 \
    -out ComputeMatrix/${tissue}_${i}_region.mat.gz \
    -p 24 \
    --skipZeros

    plotProfile -m ComputeMatrix/${tissue}_${i}_region.mat.gz \
    -out Figures/${tissue}_${i}_region.PDF \
    --perGroup \
    --colors red gray yellow green pink \
    --samplesLabel H3K27ac H3K27me3 H3K4me1 H3K4me3 ATAC \
    -z ""
    done
done
```

## PCA for each assay



```shell
#!/bin/bash
for assay in H3K27ac H3K27me3 H3K4me1 H3K4me3 RNASeq
do

cd /vol2/mengzhu/snakemake_sheep/clean/bw

ma=$(grep ${assay} labels.txt|perl -p -e 's/\n/ /g')
echo ${ma}
multiBigwigSummary bins -p=5 \
-b ${assay}_*ZScores.bw \
-out ${assay}_MultiBigwigSummary.npz \
--labels ${ma} \
--binSize=1000

cd /vol2/mengzhu/snakemake_sheep
plotCorrelation -in ${assay}_MultiBigwigSummary.npz \
--corMethod pearson \
--skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap \
-o ${assay}_Pearson_Correlation.PDF \
--removeOutliers \
--colorMap bwr \
--outFileCorMatrix ${assay}_pearsonCorr_readCounts.tab


plotPCA -in ${assay}_MultiBigwigSummary.npz \
-o ${assay}_PCA_1_vs_2.pdf \
--plotTitle "PCA" \
--transpose --PCs 1 2 \
--outFileNameData ${assay}_PCA.tab
done


#!/bin/bash
for assay in ATAC
do
cd /vol2/mengzhu/snakemake_sheep/clean/bw

ma=$(grep ${assay} labels.txt|perl -p -e 's/\n/ /g')
echo ${ma}
multiBigwigSummary bins -p=24 \
-b  ATAC_*ZScores.bw \
-out ${assay}_MultiBigwigSummary.npz \
--labels ${ma} \
--binSize=1000 \
--outRawCounts ${assay}_MultiBigwigSummary.txt


plotCorrelation -in ${assay}_MultiBigwigSummary.npz \
--corMethod pearson \
--skipZeros --plotTitle "Pearson Correlation of Read Counts" \
--whatToPlot heatmap \
-o ${assay}_Pearson_Correlation.PDF \
--removeOutliers \
--colorMap bwr \
--outFileCorMatrix ${assay}_pearsonCorr_readCounts.tab


plotPCA -in ${assay}_MultiBigwigSummary.npz \
-o ${assay}_PCA_1_vs_2.pdf \
--plotTitle "PCA" \
--transpose --PCs 1 2 \
--outFileNameData ${assay}_PCA.tab

done
```



```shell
#Prepare _PCA_last.txt
sed '1d' H3K4me3_PCA.tab > H3K4me3_PCA_no_header.tab
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j]
        }
        print str
    }
}' H3K4me3_PCA_no_header.tab > H3K4me3_PCA_transposed.tab
sed -i '$d' H3K4me3_PCA_transposed.tab
cut -f1 H3K4me3_PCA_transposed.tab | cut -d '_' -f 2 >2
cut -f1 H3K4me3_PCA_transposed.tab | cut -d '_' -f 3 >3
paste H3K4me3_PCA_transposed.tab 2 3 >H3K4me3_PCA_last.txt
#Finally, Group and Rep still need to be updated
```

**Plot Supplementary Fig. 3**

```shell
#PCA
library(ggplot2)
setwd('/vol2/mengzhu/snakemake_sheep/Figures/') 
makers =c("RNASeq", "ATAC", "H3K27ac","H3K27me3","H3K4me3","H3K4me1")
pdf("PCA_all_Assay.pdf", width=7.4, height=5)
for(maker in makers)
{
data <- read.table(paste(maker, "_PCA_last.txt", sep=""), header = T) 
  # Convert the Rep column to a factor
  data$Rep <- as.factor(data$Rep)

print(ggplot(data)+geom_point(aes(x=X1,y=X2,color= group, shape = Rep),size=3)+
        theme(legend.title =element_blank())+labs(x="PCA1",y="PCA2")+labs(title = maker)+
        scale_shape_manual(values=c(16, 17, 12, 10))+
  scale_color_manual(values=c("abomasum" = "#f18264",
      "adipose" = "#ffc4f1",
      "bone-marrow" = "#d84b4b",
      "brainstem" = "#f4d578",
      "cecum" = "#daaa6c",
      "cerebellum" = "#efd80b",
      "cerebral-cortex" = "#dcd71a",
      "cervix" = "#b6d7a9",
      "colon" = "#f2c063",
      "cornua-uteri" = "#69d683",
      "corpus-uteri" = "#80d897",
      "duodenum" = "#eb9d63",
      "epididymis" = "#70d24b",
      "heart" = "#bc58e3",
      "hippocampus" = "#f1cb05",
      "hypothalamus" = "#825e19",
      "ileum" = "#ce9639",
      "jejunum" = "#eb951c",
      "kidney" = "#4f3136",
      "liver" = "#ad8c8b",
      "lung" = "#36b5f1",
      "lymph-node" = "#c33a11",
      "mammary-gland" = "#fed9d0",
      "medulla-oblongata" = "#d0b35b",
      "midbrain" = "#f9ed19",
      "muscle" = "#a180ca",
      "omasum" = "#f8ae81",
      "optic-chiasm" = "#fece01",
      "ovary" = "#69d28c",
      "oviduct" = "#79ffaa",
      "pineal" = "#807120",
      "pituitary" = "#f1d95d",
      "pons" = "#fcd222",
      "rectum" = "#efe0a8",
      "reticulum" = "#d97c68",
      "rumen" = "#fc9891",
      "skin" = "#d09dc5",
      "soft-horn" = "#a25d73",
      "spleen" = "#962932",
      "splenium" = "#7c6919",
      "testis" = "#7ef351",
      "thymus" = "#ff3a32",
      "thyroid" = "#f359d1"))+
        theme(axis.text.x=element_text(colour="black",family="Times",size=15), #Set the x-axis tick-label font, rotate it by 15 degrees, shift it down by 1 (hjust = 1), and use Times at size 20
              axis.text.y=element_text(family="Times",size=15,face="plain"), #Set the y-axis tick-label font family, size, and plain style
              axis.title.y=element_text(family="Times",size = 15,face="plain"),
              axis.title.x=element_text(family="Times",size = 15,face="plain"))+#Set the y-axis title font properties
        theme(legend.text=element_text(family="Times", colour="black",  #Set the legend-label font properties
                                       size=16))+
        theme(legend.title=element_text(family="Times", colour="black", #Set the legend-title font properties
                                        size=16))+theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_hline(yintercept = 0, colour="blue", linetype="dashed" )+
    geom_vline(xintercept = 0, colour="blue", linetype="dashed" )
  
)
}
dev.off() 
```



# 2.ChromImpute/ChromHMM



```shell
#!/bin/bash
#This step follows peak calling and produces *_treat_pileup.bdg and *_control_lambda.bdg files for each sample
#ls *_treat_pileup.bdg | cut -d '_' -f 1-3 |while read id; do
id=$1
echo $id
chipReads=$(cat /vol2/mengzhu/snakemake_sheep/clean/bam1/unblacklist/${id}.bed | wc -l | awk '{printf "%f", $1/1000000}')
macs2 bdgcmp -t ${id}_treat_pileup.bdg -c ${id}_control_lambda.bdg -o ${id}_ppois.bdg -m ppois -S ${chipReads}
slopBed -i ${id}_ppois.bdg -g /vol2/mengzhu/genome/reference/sheep.size -b 0 | bedClip stdin /vol2/mengzhu/genome/reference/sheep.size ${id}.pval.signal.bedgraph
awk 'BEGIN{OFS="\t"}
    NR==1 {
      last_chrom=$1; last_end=$3
      print $0
      next
    }
    {
      if ($1 == last_chrom && $2 < last_end) {
        $2 = last_end
      }
      if ($2 < $3) {
        print $0
        last_chrom = $1; last_end = $3
      }
    }' ${id}.pval.signal.bedgraph > ${id}.pval.signal_1.bedgraph
sort -k1,1 -k2,2n ${id}.pval.signal_1.bedgraph > ${id}.pval.signal.bedgraph.tmp
bedGraphToBigWig ${id}.pval.signal.bedgraph.tmp /vol2/mengzhu/genome/reference/sheep.size ${id}.pval.signal.bigwig
gzip -c -f ${id}.pval.signal.bedgraph.tmp > ${id}.pval.signal.bedgraph.gz
#done
#sbatch -p low -c 8 --mem=32G -t 10-0 01_pval_signal_bedgraph.sh
#for i in cat sample.tab; do sbatch -p low -w comput3 -c 4 --mem=8G -t 10-0 01_pval_signal_bedgraph.sh $i ; done

#0.Prepare Imputed_samples.tab
mkdir CHROMHMMDIR CONVERTEDDIR DISTANCEDIR IMPUTED INPUTDATADIR PREDICTORDIR TRAINDATA
for bdg in *.pval.signal.bedgraph.gz; do
  # Extract the tissue name
  tissue=$(echo "$bdg" | cut -d'_' -f 2-3 | cut -d '.' -f 1)
  # Extract the mark name (the prefix, such as ATAC or H3K27ac)
  mark=$(echo "$bdg" | cut -d'_' -f1)
  # Get the full path
  echo -e "${tissue}\t${mark}\t${bdg}" >> Imputed_samples.tab
done


#1.Convert data to the format required by ChromImpute
vim ChromImpute_Convert1.sh 
#!/bin/bash
java -jar -Xmx40G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar Convert INPUTDATADIR  Imputed_samples.tab /vol2/mengzhu/genome/reference/sheep1.size CONVERTEDDIR


#2.Generate the training dataset
vim ChromImpute_ComputeGlobalDist.sh
#!/bin/bash
java -jar -Xmx32G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ComputeGlobalDist -m ${mark} /vol2/mengzhu/soft/ChromImpute/sheep/CONVERTEDDIR ${sample_tab} /vol2/mengzhu/genome/reference/sheep.size /vol2/mengzhu/soft/ChromImpute/sheep/DISTANCEDIR
java -jar -Xmx32G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar GenerateTrainData /vol2/mengzhu/soft/ChromImpute/sheep/CONVERTEDDIR /vol2/mengzhu/soft/ChromImpute/sheep/DISTANCEDIR ${sample_tab} /vol2/mengzhu/genome/reference/sheep.size /vol2/mengzhu/soft/ChromImpute/sheep/TRAINDATA ${mark}
for i in H3K4me1 H3K4me3 H3K27ac H3K27me3 ATAC RNASeq; do sbatch -p smp -c 4 --mem=32G -t 10-0 ChromImpute_ComputeGlobalDist.sh $i imputed_samples.tab; done

#3.Start training and generate imputed results (cannot run in parallel)
vim ChromImpute_TrainApply.sh
#!/bin/bash
set -e
module load openjdk/16.0.2
EID=$1
assay=$2
sample_tab=$3
java -jar -Xmx8G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar Train TRAINDATA ${sample_tab} PREDICTORDIR ${EID} ${assay}
java -jar -Xmx8G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar Apply CONVERTEDDIR DISTANCEDIR PREDICTORDIR ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size IMPUTED ${EID} ${assay}
#for i in `cat Imputed_samples_39.tab | awk '{print $1}' | sort | uniq`; do for j in H3K4me3 H3K27ac H3K27me3; do sbatch -p smp -c 4 --mem=8G -t 10-0 -o Logs/GeneratePredictors.${i}_${j}.%j.out ChromImpute_TrainApply.sh $i $j Imputed_samples.tab; done; done
#for i in `cat Imputed_samples_40.tab | awk '{print $1}' | sort | uniq`; do for j in ATAC H3K4me1 H3K4me3 H3K27ac H3K27me3; do sbatch -p low -w comput3 -c 4 --mem=8G -t 10-0 -o Logs/GeneratePredictors.${i}_${j}.%j.out ChromImpute_TrainApply.sh $i $j Imputed_samples_40.tab; done; done

#4.Convert predictions to the input format required for ChromHMM learning
#The -g 2 option enables binarization; the exported files are already binarized and can be passed directly to ChromHMM LearnModel
mkdir ATAC H3K27ac H3K27me3 H3K4me1 H3K4me3
vim 04_ExportToChromHMM_ATAC.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/ATAC
#for i in {1..23}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_ATAC.${i}.out 04_ExportToChromHMM_ATAC.sh Temp/export2chromhmm_ATAC.${i}.tab; done
vim 04_ExportToChromHMM_H3K4me1.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K4me1
#for i in {1..3}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K4me1.${i}.out 04_ExportToChromHMM_H3K4me1.sh Temp/export2chromhmm_H3K4me1.${i}.tab; done
vim 04_ExportToChromHMM_H3K4me3.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K4me3
#for i in {1}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K4me3.1.out 04_ExportToChromHMM_H3K4me3.sh Temp/export2chromhmm_H3K4me3.1.tab; done
vim 04_ExportToChromHMM_H3K27ac.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K27ac
#for i in {1..20}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K27ac.${i}.out 04_ExportToChromHMM_H3K27ac.sh Temp/export2chromhmm_H3K27ac.${i}.tab; done
vim 04_ExportToChromHMM_H3K27me3.sh
#!/bin/bash
set -e
sample_tab=$1
java -jar -Xmx24G /vol2/mengzhu/soft/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /vol2/mengzhu/genome/reference/sheep1.size CHROMHMMDIR/H3K27me3
#for i in {1..26}; do sbatch -p low -c 2 --mem=4G -t 10-0 -o Logs/export2chromhmm_H3K27me3.${i}.out 04_ExportToChromHMM_H3K27me3.sh Temp/export2chromhmm_H3K27me3.${i}.tab; done


#5.Merge imputed data into the binarized files generated by ChromHMM
5.1 ATAC
#!/usr/bin/env bash
set -euo pipefail

ATAC_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/ATAC"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${ATAC_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    ATAC_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $ATAC_TXT"

    awk 'NR>2{print $1}' "$ATAC_TXT" > atac_col.tmp
    tail -n +3 "$MATRIX" | cut -f2- > other_marks.tmp
    paste atac_col.tmp other_marks.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm atac_col.tmp other_marks.tmp body.tmp
done

5.2 H3K27ac
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K27ac"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${H3_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    awk 'NR>2{print $1}' "$H3_TXT" > h27ac_col.tmp
    tail -n +3 "$MATRIX" > matrix_body.tmp

    cut -f1 matrix_body.tmp > col1.tmp
    cut -f3- matrix_body.tmp > others.tmp

    paste col1.tmp h27ac_col.tmp others.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h27ac_col.tmp matrix_body.tmp col1.tmp others.tmp body.tmp
done

5.3 H3K27me3
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K27me3"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${H3_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    # Extract the 0/1 column starting at row 3 from the single-column H3K27me3 file
    awk 'NR>2{print $1}' "$H3_TXT" > h27me3_col.tmp

    # Matrix body starting at row 3
    tail -n +3 "$MATRIX" > matrix_body.tmp

    # Retain columns 1-2 of the original matrix (ATAC, H3K27ac)
    cut -f1-2 matrix_body.tmp > col12.tmp
    # Retain column 4 onward from the original matrix (H3K4me1, H3K4me3)
    cut -f4-  matrix_body.tmp > others.tmp

    paste col12.tmp h27me3_col.tmp others.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h27me3_col.tmp matrix_body.tmp col12.tmp others.tmp body.tmp
done

5.4 H3K4me1
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K4me1"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_40"

for f in "${H3_DIR}"/*40*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    # Extract the 0/1 column starting at row 3 from the single-column H3K4me1 file
    awk 'NR>2{print $1}' "$H3_TXT" > h4me1_col.tmp

    # Matrix body starting at row 3
    tail -n +3 "$MATRIX" > matrix_body.tmp

    # Retain columns 1-3 of the original matrix (ATAC, H3K27ac, H3K27me3)
    cut -f1-3 matrix_body.tmp > col123.tmp
    # Retain column 5 onward from the original matrix (H3K4me3 etc.)
    cut -f5-  matrix_body.tmp > others.tmp

    # Reassemble: ATAC, H3K27ac, H3K27me3, new H3K4me1, H3K4me3...
    paste col123.tmp h4me1_col.tmp others.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h4me1_col.tmp matrix_body.tmp col123.tmp others.tmp body.tmp
done

5.5 H3K4me3
#!/usr/bin/env bash
set -euo pipefail

H3_DIR="/vol2/mengzhu/soft/ChromImpute/sheep_clean_no_blacklist_modif/CHROMHMMDIR/H3K4me3"
MATRIX_DIR="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/SAMPLEDATA_sheep_five_maker_39"

for f in "${H3_DIR}"/*39*_chr*_binary.txt; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    noext=${base%.txt}
    sample=${noext%_chr*_binary}
    chr=${noext#*_chr}
    chr=${chr%_binary}

    H3_TXT="$f"
    MATRIX="${MATRIX_DIR}/${sample}_chr${chr}_binary.txt"

    if [[ ! -s "$MATRIX" ]]; then
        echo "[WARN] skip: MATRIX not found or empty: $MATRIX"
        continue
    fi

    echo "Fixing $MATRIX using $H3_TXT"

    awk 'NR>2{print $1}' "$H3_TXT" > h4me3_col.tmp
    tail -n +3 "$MATRIX" > matrix_body.tmp

    cut -f1-4 matrix_body.tmp > col1234.tmp

    paste col1234.tmp h4me3_col.tmp > body.tmp

    {
        head -n2 "$MATRIX"
        cat body.tmp
    } > tmp && mv tmp "$MATRIX"

    rm h4me3_col.tmp matrix_body.tmp col1234.tmp body.tmp
done



#5.Use ChromHMM to predict chromatin states
vim ChromHMM_ModelOptm.sh
#!/bin/bash
# Script: ChromHMM_ModelOptm.sh
mkdir -p ChroHMM/LearnModel_${num_model}
num_model=$1
java -jar -Xmx64G /vol2/mengzhu/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /vol2/mengzhu/genome/reference/sheep1.size /vol2/mengzhu/soft/ChromImpute/sheep/CHROMHMMDIR ChroHMM/LearnModel_${num_model} ${num_model} Ramb2
# for i in {10..18}; do sbatch -p low -c 12 --mem=48G -t 10-0 -o LearnModel_${i}.%j.out ChromHMM_ModelOptm.sh $i; done

#6.Generate 2-20 chromatin-state solutions for Rep1 and Rep2
vim 02_LearnModel.sh
#!/bin/bash
num_model=$1
mkdir Rep1/LearnModel_${num_model}
mkdir Rep2/LearnModel_${num_model}
java -jar -Xmx22G /vol2/mengzhu/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /vol2/mengzhu/genome/reference/sheep1.size SAMPLEDATA_sheep_five_maker_39 Rep1/LearnModel_${num_model} ${num_model} Ramb2
java -jar -Xmx22G /vol2/mengzhu/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /vol2/mengzhu/genome/reference/sheep1.size SAMPLEDATA_sheep_five_maker_40 Rep2/LearnModel_${num_model} ${num_model} Ramb2
# for i in {2..20}; do sbatch -p low -c 11 --mem=22G -t 10-0 -o Logs/LearnModel_${i}.out 02_CompareModels.sh ${i}; done


#Plot correlations among the 2-20 chromatin-state solutions and select the optimal number of states
mkdir emissions_Rep1
cp LearnModel_*/emissions_*.txt emissions_Rep1
mkdir emissions_Rep2
cp LearnModel_*/emissions_*.txt emissions_Rep2
vim 03_CompareModels_Rep2.sh
#!/bin/bash
#four maker and 16 chromatin state 39
java -mx80000M -jar /vol2/mengzhu/ChromHMM/ChromHMM.jar CompareModels \
    -color 255,0,0 \
    /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/LearnModel_20/emissions_20.txt \
    /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/emissions_Rep2 \
    /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/CompareModels/Rep2
#sbatch -c 40  -p low --mem 80G 03_CompareModels_Rep2.sh
# Extract the first column, remove rows with the value "1", transpose, and save as the first row
cut -f1 Rep1.txt | grep -v "^1$" | tr '\n' '\t' > Rep1_average.txt
echo "" >> Rep1_average.txt  # Append a newline

# Calculate the mean of each column and insert "average" before the first column of the second row
awk '
NR > 1 {
    for(i=2; i<=NF; i++) {
        sum[i] += $i;
    }
    count++;
}
END {
    printf "average\t";  # Insert "average" before the first column of the second output row
    for(i=2; i<=NF; i++) {
        printf "%.2f\t", sum[i]/count;
    }
    print "";  # Newline
}' Rep1.txt >> Rep1_average.txt
```

# 3.TSR

```sh
#Split by chromatin state
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep1/LearnModel_11/new_LearnModel_11
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/LearnModel_11/new_LearnModel_11
ls *_11_segments.bed | while read id;
do
  echo $id
  b=$(basename $id "_11_segments.bed")
  echo $b
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11
  do
     grep -w $state $id | sort -k1,1 -k2,2n > "state_variability/"$b"_"$state".bed"
  done
done

#Merge Rep1 and Rep2
vim 01_merge.sh
#!/usr/bin/env bash
while read -r id; do
  echo "$id"
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11; do   # ← do is present
    cat \
      "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep1/LearnModel_11/new_LearnModel_11/state_variability/${id}_39_${state}.bed" \
      "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/LearnModel_11/new_LearnModel_11/state_variability/${id}_40_${state}.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -c 4 -o distinct \
      > "${id}_${state}.bed"
  done
done < sample.txt
#sbatch -c 5  -p low --mem 10G 01_merge.sh

#Gs
####get Gs: total regions of each state (Gs) across 14 tissues 
####Merge all tissues by chromatin state
vim 02_AAGs.sh
#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability
mkdir AAGs
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11
do
echo $i
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done
#sbatch -c 5  -p low --mem 10G 02_AAGs.sh

#Identify tissue-specific chromatin states
#Intersect each chromatin state from each tissue with the merged chromatin states; write the number of overlaps in column 4, or 0 when there is no overlap
vim 03_AARegulatory_module.sh
#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability
mkdir AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11
do
ls *${state}.bed | while read id;
do
echo $id
bedtools intersect -a <(sort -k1,1 -k2,2n /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs/${state}_Gs.bed) -b <(sort -k1,1 -k2,2n ${id}) -c -sorted >  "AARegulatory_module/"${id%%.*}"_Gs.bed" #-coutputs each record from A and appends a column containing the number of overlapping records in B (an integer; 0 indicates no overlap
done
done
#sbatch -c 5  -p low --mem 10G 03_AARegulatory_module.sh

cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability
rm 1.txt
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
echo chr start end > 1.txt
ls *_${state}.bed >> 1.txt
cat 1.txt | cut -d "_" -f 1 | perl -p -e 's/\n/ /g'| sed '$ s/.$/\n/' > AARegulatory_module/header.txt
#paste all tissues together
done

vim 01_state_merge.sh
#!/bin/bash
#Record the presence of all chromatin reads as 1 and their absence as 0
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
ls *${state}_Gs.bed  | while read id;
do
  echo $id
  cat ${id} | cut -f 4 | paste -s >> 2.txt
  done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 *_${state}_Gs.bed > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' >6.txt
cat header.txt 6.txt |sed 's/ /\t/g' > all_${state}_Gs.csv
done
#sbatch -c 5  -p low --mem 10G 01_state_merge.sh

#normalized the one count and count the number for each region
vim 02_normalization.sh
#!/bin/bash
#Normalize one count and count the number for each region
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
echo all_${state}_Gs.csv
cat all_${state}_Gs.csv | cut -f 4- |  sed 's/40/1/g' | sed 's/41/1/g'| sed 's/42/1/g'| sed 's/43/1/g'|  sed 's/30/1/g' | sed 's/31/1/g'| sed 's/32/1/g'| sed 's/33/1/g'|sed 's/34/1/g'| sed 's/35/1/g'| sed 's/36/1/g'| sed 's/37/1/g'| sed 's/38/1/g'| sed 's/39/1/g'| sed 's/20/1/g' | sed 's/21/1/g'| sed 's/22/1/g'| sed 's/23/1/g'|sed 's/24/1/g'| sed 's/25/1/g'| sed 's/26/1/g'| sed 's/27/1/g'| sed 's/28/1/g'| sed 's/29/1/g'| sed 's/10/1/g'| sed 's/11/1/g'| sed 's/12/1/g'| sed 's/13/1/g'| sed 's/14/1/g'| sed 's/15/1/g' | sed 's/16/1/g'| sed 's/17/1/g'|  sed 's/18/1/g'| sed 's/19/1/g'| sed 's/2/1/g'| sed 's/3/1/g'| sed 's/4/1/g'| sed 's/5/1/g' | sed 's/6/1/g'| sed 's/7/1/g'|  sed 's/8/1/g'| sed 's/9/1/g' > 3.txt
cat 3.txt | awk '{for(i=1;i<=NF;i++){a[NR]+=$i}print $0,a[NR]}' > 4.txt
cut  -f 1-3 all_${state}_Gs.csv > 5.txt
paste  5.txt 4.txt |sed 's/ /\t/g' > all_${state}_Gs_one_count.csv
done
#sbatch -c 5  -p low --mem 10G 02_normalization.sh

vim 03_TSR_tissue.sh
#!/bin/bash
#Because there are more than 30 tissues, tissues from the same system are treated as one tissue and are not excluded when identifying tissue-specific regulatory elements.
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
echo $state
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
mkdir AA_TSR_${state}
cat all_${state}_Gs_one_count.csv | awk '{if (($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30+$31+$32+$33+$34+$35+$36+$37+$38+$39+$40+$41+$42+$43+$44+$45+$46)==43) print $0}' > AA_TSR_${state}/TSR_All_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43)==12 \
&& ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_Nervous_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($39+$38+$30+$4)==4 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Stomach_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($15+$21+$20+$8+$12+$37)==6 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Digestive_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($11+$13+$14)==3 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Uterus_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($6+$25+$45+$46+$42)==5 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Immune_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($17+$29)==2 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Muscular_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($40+$41)==2 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5)==0) print $0}' \
> AA_TSR_${state}/TSR_Skin_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($32+$33+$26)==3 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Female-reproductive_common_${state}.txt

cat all_${state}_Gs_one_count.csv | awk '{if (($16+$44)==2 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$26+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Male-reproductive_common_${state}.txt

cat all_${state}_Gs_one_count.csv | awk '{if ($10==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_cerebral-cortex_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($28==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_midbrain_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($9==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_cerebellum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($7==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_brainstem_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($18==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_hippocampus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($19==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_hypothalamus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($27==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_medulla-oblongata_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($31==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_optic-chiasm_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($34==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_pineal_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($35==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_pituitary_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($36==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_pons_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($43==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_splenium_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($39==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_rumen_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($38==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_reticulum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($30==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_omasum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($4==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_abomasum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($15==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_duodenum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($21==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_jejunum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($20==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_ileum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($8==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_cecum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($12==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_colon_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($37==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_rectum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($11==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$13+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_cervix_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($13==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_cornua-uteri_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($14==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$13+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_corpus-uteri_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($32==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_ovary_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($33==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_oviduct_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($16==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$26+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_epididymis_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($44==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$26+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_testis_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($26==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_mammary-gland_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($6==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_bone-marrow_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($25==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_lymph-node_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($45==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_thymus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($46==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_thyroid_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($42==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_spleen_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($23==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_liver_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($22==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_kidney_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($24==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_lung_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($17==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_heart_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($29==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_muscle_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($5==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_adipose_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($40==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5)==0) print $0}' > AA_TSR_${state}/TSR_skin_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($41==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5)==0) print $0}' > AA_TSR_${state}/TSR_soft-horn_${state}.txt
done
#sbatch -c 100  -p low --mem 200G 03_TSR_tissue_10.sh
```

## TSR summary

```shell
vim 01_AA_Gs_summary.sh
#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AAGs
rm AA_Gs_summary.txt
ls *Gs.bed | while read id;
do
A=$(cat $id | wc -l | awk '{print $1}')
B=$(cat $id | awk '{print ($3-$2) }' | awk '{sum+=$1} END {print sum/NR}')
C=$(cat $id | awk '{print ($3-$2) }' | awk '{sum+=$1} END {print sum}')
echo $(basename $id "_Gs.bed") $A $B $C $(bc <<< "scale=3;($C/2628104905)") | sed 's/ /\t/g' >> AA_Gs_summary.txt
done
#sbatch -c 5  -p low --mem 10G AA_Gs_summary.sh

#!/bin/bash
cd /vol2/mengzhu/ChromHMM/new/OUTPUTSAMPLE_sheep_18_modify/state_variability/AARegulatory_module/all
for tissue in abomasum adipose bone-marrow brainstem cecum cerebellum cerebral-cortex cervix colon cornua-uteri corpus-uteri duodenum epididymis heart hippocampus hypothalamus ileum jejunum kidney liver lung lymph-node mammary-gland medulla-oblongata midbrain muscle omasum optic-chiasm ovary oviduct pineal pituitary pons rectum reticulum rumen skin soft-horn spleen splenium testis thymus thyroid Digestive_common Immune_common Nervous_common
do 
ls  TSR_${tissue}*_id.bed | while read id;
do
A=$(cat $id | wc -l | awk '{print $1}')
B=$(cat $id | awk '{print ($3-$2) }' | awk '{sum+=$1} END {print sum/NR}')
C=$(cat $id | awk '{print ($3-$2) }' | awk '{sum+=$1} END {print sum}')

echo $(basename $id "_id.bed") $A $B $C $(bc <<< "scale=3;($C/2628104905)") | sed 's/ /\t/g' >> 1.txt
done
echo State $tissue $tissue $tissue  $tissue > 2.txt
#paste 4.txt 3.txt | sed 's/ /\t/g' > 6.txt
cat 2.txt 1.txt |  sed 's/ /\t/g' > ${tissue}_Gs_summary.txt
rm 1.txt
done


for tissue in soft-born
do 
ls  TSR_${tissue}*_id.bed | while read id;
do
A=$(cat $id | wc -l | awk '{print $1}')
B=$(cat $id | awk '{print ($3-$2) }' | awk '{sum+=$1} END {print sum/NR}')
C=$(cat $id | awk '{print ($3-$2) }' | awk '{sum+=$1} END {print sum}')

echo $(basename $id "_id.bed") $A $B $C $(bc <<< "scale=3;($C/2628104905)") | sed 's/ /\t/g' >> 1.txt
done
echo State $tissue $tissue $tissue  $tissue > 2.txt
#paste 4.txt 3.txt | sed 's/ /\t/g' > 6.txt
cat 2.txt 1.txt |  sed 's/ /\t/g' > ${tissue}_Gs_summary.txt
rm 1.txt
done

for file in *_Gs_summary.txt; do
  cut -f2 "$file" >> number.txt
done

for file in *_Gs_summary.txt; do
  cut -f3 "$file" >> size.txt
done

for file in *_Gs_summary.txt; do
  cut -f5 "$file" >> genome_coverage.txt
done

files=(*_Gs_summary.txt)

# Get the number of lines in the file
num_lines=$(wc -l < "${files[0]}")

# Iterate over each line
for ((i=1; i<=num_lines; i++)); do
  # Extract row i from column 2 of each file
  line=""
  for file in "${files[@]}"; do
    line="$line$(sed -n "${i}p" "$file" | cut -f2)	"  # Use tabs as delimiters
  done
  # Remove the trailing tab from each line and write to merged_columns.txt
  echo -e "$line" >> number.txt
done
```

## TSR link to target gene



```shell
#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
    ls AA_TSR_${state}/TSR_*_${state}.txt | while read id;
    do
      echo $id

        cat $id|  cut -f 1-2 | sed 's/\t/:/g' > 1.txt
        cat $id | cut -f 3 | sed 's/^/-/g'  > 2.txt
        paste 1.txt 2.txt | sed 's/\t//g' > 3.txt
        cat $id|  cut -f 1-3 > 4.txt 
        paste 4.txt 3.txt > AA_TSR_${state}/$(basename $id ".txt")_id.bed
    done
done

vim 04_sheep_to_human.sh
#!/bin/bash
#Convert sheep tissue-specific enhancer regions to human coordinates
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11
do
  ls AA_TSR_${state}/TSR_*_${state}_id.bed | while read id;
    do
      echo $id
      mkdir -p AA_TSR_${state}/Lifttohuman

      out1=AA_TSR_${state}/Lifttohuman/$(basename $id ".bed")_id_hg38_lift.bed
      out2=AA_TSR_${state}/Lifttohuman/$(basename $id ".bed")_id_hg38_unlift.bed
      echo ${out1}
     liftOver -minMatch=0.1 $id /vol2/mengzhu/genome/GCF_016772045.2ToHg38_chr.over.chain.gz ${out1} ${out2}
  done
done
#sbatch -c 5  -p low --mem 10G 04_sheep_to_human.sh
```

## Promoter

```shell
#Try merging chromatin states E1, E2, E3, and E9 before identifying target genes (first inspect the target-gene results for each state)
cat E1_Gs.bed E2_Gs.bed E3_Gs.bed E3_Gs.bed |sort -k1,1 -k2,2n > promoter.bed
bedtools merge -i promoter.bed > promoter_merge.bed
cat E4_Gs.bed E5_Gs.bed E6_Gs.bed E7_Gs.bed |sort -k1,1 -k2,2n > enhancer.bed
bedtools merge -i enhancer.bed > enhancer_merge.bed
ls *_E1_Gs.bed | cut -d "_" -f 1 | while read id ; do
echo ${id}
cat AA_TSR_E1/TSR_${id}_E1_id.bed AA_TSR_E2/TSR_${id}_E2_id.bed AA_TSR_E3/TSR_${id}_E3_id.bed AA_TSR_E9/TSR_${id}_E9_id.bed |sort -k1,1 -k2,2n > promoter.bed
bedtools merge -i promoter.bed > promoter_${id}.bed
done

####Use bedtools intersect to identify the gene nearest each promoter as its target gene#####
mkdir -p AA_Target_gene/Sheep_TSSup2k
ls *_id.bed | while read id; do
    echo "$id"
    bedtools intersect -a "$id" \
        -b /vol2/mengzhu/genome/part_change_esemb100/TSS_esemble100_colin.bed_up2k.bed \
        -wa -wb > AA_Target_gene/Sheep_TSSup2k/"$(basename "$id" "_id.bed")"_gene_up2k.txt
done
cd AA_Target_gene/Sheep_TSSup2k
ls *_gene_up2k.txt | while read id;
do
echo $id
cat $id | cut -f 8 | sort > $(basename $id "_gene_up2k.txt")_gene.txt 
done

cp AA_TSR_E1/AA_Target_gene/Sheep_TSSup2k/*_gene.txt AA_TSR_E2/AA_Target_gene/Sheep_TSSup2k/*_gene.txt AA_TSR_E3/AA_Target_gene/Sheep_TSSup2k/*_gene.txt AA_TSR_E4/AA_Target_gene/Sheep_TSSup2k/*_gene.txt /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/promoter_TSSup2k

ls *_E1_gene.txt | while read id ; do
echo ${id}
cat $(basename $id "_E1_gene.txt")_E1_gene.txt $(basename $id "_E1_gene.txt")_E2_gene.txt $(basename $id "_E1_gene.txt")_E3_gene.txt $(basename $id "_E1_gene.txt")_E4_gene.txt \
  | sort -u > $(basename $id "_E1_gene.txt")_promoter.bed
done



```

## Enhancer

```shell
vim 01_Target_gene_TSR_summary.sh
#!/bin/bash
######Identify genes linked to tissue-specific enhancers by overlap#####
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5
mkdir -p AA_Target_gene/Human
rm 1.txt
ls *_id.bed | while read id;
do
echo $id
join -1 4 -2 1  <(sort -k 4 ${id}) <(sort -k 1 /vol2/zhangshiwen/sheep_cor/h3k27ac/H3K27ac_output_E5_confident.tsv) |sed 's/ /\t/g'  > AA_Target_gene/$(basename $id "_id.bed")_target_gene.txt
cat AA_Target_gene/$(basename $id "_id.bed")_target_gene.txt | cut -f 7 | awk '!seen[$1]++' > AA_Target_gene/$(basename $id "_id.bed")_gene_sheep.txt
join -1 1 -2 1 <(sort -k 1,1 AA_Target_gene/$(basename $id "_id.bed")_gene_sheep.txt)    <(sort -k 1,1 /vol2/mengzhu/genome/conservation_1/sheep_human_esemble_ID.txt) | sed 's/ /\t/g' > AA_Target_gene/$(basename $id "_id.bed")_gene_human.txt
cat AA_Target_gene/$(basename $id "_id.bed")_gene_human.txt | cut -f 2 > AA_Target_gene/Human/$(basename $id "_id.bed")_gene_human.txt 
A=$(cat AA_Target_gene/$(basename $id "_id.bed")_target_gene.txt |  wc -l | awk '{print $1}')
B=$(cat $id |  wc -l | awk '{print $1}')
C=$(cat AA_Target_gene/$(basename $id "_id.bed")_target_gene.txt | cut -f 1| awk '!seen[$1]++' | wc -l | awk '{print $1}')
D=$(cat AA_Target_gene/$(basename $id "_id.bed")_gene_sheep.txt  |  wc -l | awk '{print $1}')
E=$(cat AA_Target_gene/Human/$(basename $id "_id.bed")_gene_human.txt |  wc -l | awk '{print $1}')
echo $(basename $id "_id.bed") $A $B $C $(bc <<< "scale=10;($C/$B)") $D $E | sed 's/ /\t/g'>> 1.txt
done
echo sample pair origin_enhancer target_enhancer ratio target_gene change_to_human |  sed 's/ /\t/g' > 3.txt
cat 3.txt 1.txt > Target_gene_TSR_summary.csv
#sbatch -c 5  -p low --mem 10G 01_Target_gene_TSR_summary.sh

cd AA_Target_gene
ls TSR_*_gene_sheep.txt | while read id;
do
echo $id
sort $id > $(basename $id "_sheep.txt")_sheep_sorted.txt
done

cp AA_TSR_E5/AA_Target_gene/*_gene_sheep_sorted.txt AA_TSR_E6/AA_Target_gene/*_gene_sheep_sorted.txt AA_TSR_E7/AA_Target_gene/*_gene_sheep_sorted.txt AA_TSR_E8/AA_Target_gene/*_gene_sheep_sorted.txt enhancer

ls *_E5_gene_sheep_sorted.txt | while read id ; do
echo ${id}
cat $(basename $id "_E5_gene_sheep_sorted.txt")_E5_gene_sheep_sorted.txt $(basename $id "_E5_gene_sheep_sorted.txt")_E6_gene_sheep_sorted.txt $(basename $id "_E5_gene_sheep_sorted.txt")_E7_gene_sheep_sorted.txt $(basename $id "_E5_gene_sheep_sorted.txt")_E8_gene_sheep_sorted.txt \
  | sort -u > $(basename $id "_E5_gene_sheep_sorted.txt")_enhancer.bed
done
```

#### Enhancer/Promoter to ID

The ENCODE project defines a unique cCRE ID for each cCRE: EH38E1393970
EH38: is the reference-genome version
E: is shorthand for enhancer
1393970: is a seven-digit identifier ranging from 0000001 to 999999

Following this principle, sheep cCRE IDs can be written as follows (without the reference genome)

Ovis denotes Ovis aries

E5 denotes chromatin state 5
Identifiers range from 000001 to 999999

OvisE5000001-OvisE5999999

```shell
####Assign IDs to the merged chromatin states#####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs####
####Enhancer####
awk 'BEGIN{OFS="\t"}{
  region=$1":"$2"-"$3;
  id=sprintf("OvisE5%06d", ++i);
  print $1,$2,$3,region,id
}' E5_Gs.bed > E5_Gs_ID.bed

awk 'BEGIN{FS=OFS="\t"}
NR==FNR { id[$4]=$5; next }                    # Read E5_Gs_ID.bed first: key=$4, value=$5
FNR==1 {                                       # Process the header: insert the new column name after Enhancer
  printf "%s\t%s\t%s\tEnhancer_ID", $1,$2,$3
  for(i=4;i<=NF;i++) printf "\t%s", $i
  printf "\n"
  next
}
{
  eid = (($3 in id) ? id[$3] : "NA")           # Enhancer is in column 3
  printf "%s\t%s\t%s\t%s", $1,$2,$3,eid
  for(i=4;i<=NF;i++) printf "\t%s", $i
  printf "\n"
}' /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs/E5_Gs_ID.bed all_tables.with_header.tsv > all_tables.with_enhancerID.tsv

####Promoter####
awk 'BEGIN{OFS="\t"}{
  region=$1":"$2"-"$3;
  id=sprintf("OvisP1%06d", ++i);
  print $1,$2,$3,region,id
}' E1_Gs.bed > E1_Gs_ID.bed

awk 'BEGIN{FS=OFS="\t"}
NR==FNR { p[$4]=$5; next }   # Read E1_Gs_ID.bed first: key=$4, value=$5

FNR==1 {                    # Header: insert Promoter_ID after Promoter
  for(i=1;i<=5;i++) printf (i==1? "%s":"\t%s"), $i
  printf "\tPromoter_ID"
  for(i=6;i<=NF;i++) printf "\t%s", $i
  printf "\n"
  next
}

{
  pid = (($5 in p) ? p[$5] : "NA")   # Promoter is in column 5
  for(i=1;i<=5;i++) printf (i==1? "%s":"\t%s"), $i
  printf "\t%s", pid
  for(i=6;i<=NF;i++) printf "\t%s", $i
  printf "\n"
}' /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs/E1_Gs_ID.bed all_tables.with_enhancerID.tsv > all_tables.with_enhancer_promoterID.tsv
```

## GREAT predicts functions of cis-regulatory regions.

```shell
71_liftover_enhancer.sh
cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E4
ls *_E4.txt | while read id;
do
cat $id | cut -f 1-3 > 1.txt
cat $id | cut -f 1-3 | sed 's/\t/vsss/g' > 2.txt
paste 1.txt 2.txt > $(basename $id "_E4.txt").bed
done
###Convert sheep regions to human coordinates
vim 04_liftsheeptohuman.sh
#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E4
ls *.bed | while read id;
do
echo ${id}
liftOver -minMatch=0.1 ${id} /vol2/mengzhu/genome/GCF_016772045.2ToHg38_chr.over.chain.gz AA_liftsheeptohuman/lifted_$(basename $id ".bed").bed AA_liftsheeptohuman/$(basename $id ".bed")_unlifted.bed 
done
#sbatch -c 5  -p low --mem 10G 04_liftsheeptohuman.sh
```

## TSR GO

```shell
TSR GO
#Extract the Term ID, Term Name, and Binom FDR Q-Val columns
cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/GO_GREAT
ls shown-GOBiologicalProcess_*.tsv | while read id;
do
aa=${id#*_}
bb=${aa%%.*}_Gorich.txt

echo $aa $bb
cat $id |sed '1d'| cut -f 2 > 1.txt
cat $id |sed '1d'| cut -f 1,3 > 3.txt
paste 1.txt 3.txt > Gorich/$bb

done

#Merge all tissues and remove duplicates in the first column
cat abomasum_Gorich.txt adipose_Gorich.txt bone-marrow_Gorich.txt brainstem_Gorich.txt cecum_Gorich.txt cerebellum_Gorich.txt cerebral-cortex_Gorich.txt cervix_Gorich.txt colon_Gorich.txt cornua-uteri_Gorich.txt corpus-uteri_Gorich.txt duodenum_Gorich.txt epididymis_Gorich.txt heart_Gorich.txt hippocampus_Gorich.txt hypothalamus_Gorich.txt ileum_Gorich.txt jejunum_Gorich.txt kidney_Gorich.txt liver_Gorich.txt lung_Gorich.txt lymph-node_Gorich.txt mammary-gland_Gorich.txt medulla-oblongata_Gorich.txt midbrain_Gorich.txt muscle_Gorich.txt omasum_Gorich.txt optic-chiasm_Gorich.txt ovary_Gorich.txt oviduct_Gorich.txt pineal_Gorich.txt pituitary_Gorich.txt pons_Gorich.txt rectum_Gorich.txt reticulum_Gorich.txt rumen_Gorich.txt skin_Gorich.txt soft-horn_Gorich.txt spleen_Gorich.txt splenium_Gorich.txt testis_Gorich.txt thymus_Gorich.txt thyroid_Gorich.txt| awk '!seen[$1]++' > AA_total_Go_orgin.txt
cut -f 1 AA_total_Go_orgin.txt  | sed '1d' > AA_total_Go_list.txt


cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/Gorich
rm 5.txt
cat AA_total_Go_list.txt | while read state;
do
  echo $state
    ls *Gorich.txt | while read tissue;
    do
    if [ "`grep -w $state $tissue`" ]
    then
      grep -w $state $tissue > 1.txt
      count=$(cat 1.txt | awk -F"\t" '{a = -log($3)/log(10); printf("%0.4f\n",a)}')

      echo $state $tissue $count >> 5.txt
     else
     let count=0
     echo $state $tissue $count >> 5.txt
    fi
  done
done




for tissue in abomasum_Gorich.txt adipose_Gorich.txt bone-marrow_Gorich.txt brainstem_Gorich.txt cecum_Gorich.txt cerebellum_Gorich.txt cerebral-cortex_Gorich.txt cervix_Gorich.txt colon_Gorich.txt cornua-uteri_Gorich.txt corpus-uteri_Gorich.txt duodenum_Gorich.txt epididymis_Gorich.txt heart_Gorich.txt hippocampus_Gorich.txt hypothalamus_Gorich.txt ileum_Gorich.txt jejunum_Gorich.txt kidney_Gorich.txt liver_Gorich.txt lung_Gorich.txt lymph-node_Gorich.txt mammary-gland_Gorich.txt medulla-oblongata_Gorich.txt midbrain_Gorich.txt muscle_Gorich.txt omasum_Gorich.txt optic-chiasm_Gorich.txt ovary_Gorich.txt oviduct_Gorich.txt pineal_Gorich.txt pituitary_Gorich.txt pons_Gorich.txt rectum_Gorich.txt reticulum_Gorich.txt rumen_Gorich.txt skin_Gorich.txt soft-horn_Gorich.txt spleen_Gorich.txt splenium_Gorich.txt testis_Gorich.txt thymus_Gorich.txt thyroid_Gorich.txt
do
  echo $tissue
  grep -w $tissue 5.txt | cut -d " " -f 3 | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
grep -w thyroid_Gorich.txt 5.txt | cut -d " " -f 1 > 4.txt
echo Go abomasum adipose bone-marrow brainstem cecum cerebellum cerebral-cortex cervix colon cornua-uteri corpus-uteri duodenum epididymis heart hippocampus hypothalarmus ileum jejunum kidney liver lung lymph-node mammary-gland medulla-oblongata midbrain muscle omasum optic-chiasm ovary oviduct pineal pituitary pons rectum reticulum rumen skin soft-horn spleen splenium testis thymus thyroid > 7.txt
paste 4.txt 3.txt | sed 's/ /\t/g' > 6.txt
cat 7.txt 6.txt |  sed 's/ /\t/g'> 9.txt

cut -f 1-2 AA_total_Go_orgin.txt  > 8.txt
paste 8.txt 9.txt > TSR_go_enhancer_enrichment.csv


```

## Human Phenotype

```shell
cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman
ls shown-HumanPhenotypeOntology_*.tsv | while read id;
do
aa=${id#*_}
bb=${aa%%.*}_HPrich.txt

echo $aa $bb
cat $id |sed '1d'| cut -f 2 > 1.txt
cat $id |sed '1d'| cut -f 1,3 > 3.txt
paste 1.txt 3.txt > HPrich/$bb
done


cat abomasum_HPrich.txt adipose_HPrich.txt bone-marrow_HPrich.txt brainstem_HPrich.txt cecum_HPrich.txt cerebellum_HPrich.txt cerebral-cortex_HPrich.txt cervix_HPrich.txt colon_HPrich.txt cornua-uteri_HPrich.txt corpus-uteri_HPrich.txt duodenum_HPrich.txt epididymis_HPrich.txt heart_HPrich.txt hippocampus_HPrich.txt hypothalarmus_HPrich.txt ileum_HPrich.txt jejunum_HPrich.txt kidney_HPrich.txt liver_HPrich.txt lung_HPrich.txt lymph-node_HPrich.txt mammary-gland_HPrich.txt medulla-oblongata_HPrich.txt midbrain_HPrich.txt muscle_HPrich.txt omasum_HPrich.txt optic-chiasm_HPrich.txt ovary_HPrich.txt oviduct_HPrich.txt pineal_HPrich.txt pituitary_HPrich.txt pons_HPrich.txt rectum_HPrich.txt reticulum_HPrich.txt rumen_HPrich.txt skin_HPrich.txt spleen_HPrich.txt splenium_HPrich.txt testis_HPrich.txt thymus_HPrich.txt thyroid_HPrich.txt| awk '!seen[$1]++' > AA_total_Go_orgin_HP.txt
cut -f 1 AA_total_Go_orgin_HP.txt  | sed '1d'> AA_total_Go_list_HP.txt

#Extract human-phenotype Q-values for all tissues in order; replace missing values with 0
cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/HPrich
rm 5.txt
cat AA_total_Go_list_HP.txt | while read state;
do
  echo $state
    ls *HPrich.txt | while read tissue;
    do
    if [ "`grep -w $state $tissue`" ]
    then
      grep -w $state $tissue > 1.txt
      count=$(cat 1.txt | awk -F"\t" '{a = -log($3)/log(10); printf("%0.4f\n",a)}')

      echo $state $tissue $count >> 5.txt
     else
     let count=0
     echo $state $tissue $count >> 5.txt
    fi
  done
done

#Combine all Q-values into 2.txt in tissue order
rm 2.txt
for tissue in abomasum_HPrich.txt adipose_HPrich.txt bone-marrow_HPrich.txt brainstem_HPrich.txt cecum_HPrich.txt cerebellum_HPrich.txt cerebral-cortex_HPrich.txt cervix_HPrich.txt colon_HPrich.txt cornua-uteri_HPrich.txt corpus-uteri_HPrich.txt duodenum_HPrich.txt epididymis_HPrich.txt heart_HPrich.txt hippocampus_HPrich.txt hypothalarmus_HPrich.txt ileum_HPrich.txt jejunum_HPrich.txt kidney_HPrich.txt liver_HPrich.txt lung_HPrich.txt lymph-node_HPrich.txt mammary-gland_HPrich.txt medulla-oblongata_HPrich.txt midbrain_HPrich.txt muscle_HPrich.txt omasum_HPrich.txt optic-chiasm_HPrich.txt ovary_HPrich.txt oviduct_HPrich.txt pineal_HPrich.txt pituitary_HPrich.txt pons_HPrich.txt rectum_HPrich.txt reticulum_HPrich.txt rumen_HPrich.txt skin_HPrich.txt soft-horn_HPrich.txt spleen_HPrich.txt splenium_HPrich.txt testis_HPrich.txt thymus_HPrich.txt thyroid_HPrich.txt
do
  echo $tissue
  grep -w $tissue 5.txt | cut -d " " -f 3 | paste -s >> 2.txt
done
ls *_HPrich.txt | while read tissue; do
  echo $tissue
  grep -w $tissue 5.txt | cut -d " " -f 3 | paste -s >> 2.txt
done

#Transpose
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
grep -w ovary_HPrich.txt 5.txt | cut -d " " -f 1 > 4.txt
echo Go abomasum adipose bone-marrow brainstem cecum cerebellum cerebral-cortex cervix colon cornua-uteri corpus-uteri duodenum epididymis heart hippocampus hypothalarmus ileum jejunum kidney liver lung lymph-node mammary-gland medulla-oblongata midbrain muscle omasum optic-chiasm ovary oviduct pineal pituitary pons rectum reticulum rumen skin soft-horn spleen splenium testis thymus thyroid > 7.txt
paste 4.txt 3.txt | sed 's/ /\t/g' > 6.txt
cat 7.txt 6.txt |  sed 's/ /\t/g'> 9.txt

cut -f 1-2 AA_total_Go_orgin_HP.txt  > 8.txt
paste 8.txt 9.txt > TSR_go_HP_enrichment.csv

```

## mouse Phenotype

```shell
cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman
ls *MGIPhenotype* | while read id;
do
aa=${id#*_}
bb=${aa%%.*}_MPrich.txt

echo $aa $bb
cat $id |sed '1d'| cut -f 2 > 1.txt
cat $id |sed '1d'| cut -f 1,3 > 3.txt
paste 1.txt 3.txt > MPrich/$bb
done

cat abomasum_MPrich.txt adipose_MPrich.txt bone-marrow_MPrich.txt brainstem_MPrich.txt cecum_MPrich.txt cerebellum_MPrich.txt cerebral-cortex_MPrich.txt cervix_MPrich.txt colon_MPrich.txt cornua-uteri_MPrich.txt corpus-uteri_MPrich.txt duodenum_MPrich.txt epididymis_MPrich.txt heart_MPrich.txt hippocampus_MPrich.txt hypothalarmus_MPrich.txt ileum_MPrich.txt jejunum_MPrich.txt kidney_MPrich.txt liver_MPrich.txt lung_MPrich.txt lymph-node_MPrich.txt mammary-gland_MPrich.txt medulla-oblongata_MPrich.txt midbrain_MPrich.txt muscle_MPrich.txt omasum_MPrich.txt optic-chiasm_MPrich.txt ovary_MPrich.txt oviduct_MPrich.txt pineal_MPrich.txt pituitary_MPrich.txt pons_MPrich.txt rectum_MPrich.txt reticulum_MPrich.txt rumen_MPrich.txt skin_MPrich.txt soft-horn_MPrich.txt spleen_MPrich.txt splenium_MPrich.txt testis_MPrich.txt thymus_MPrich.txt thyroid_MPrich.txt | awk '!seen[$1]++' > AA_total_Go_orgin_MP.txt
cut -f 1 AA_total_Go_orgin_MP.txt  | sed '1d'> AA_total_Go_list_MP.txt


cd /vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/MPrich
rm 5.txt
cat AA_total_Go_list_MP.txt | while read state;
do
  echo $state
    for tissue in abomasum_MPrich.txt adipose_MPrich.txt bone-marrow_MPrich.txt brainstem_MPrich.txt cecum_MPrich.txt cerebellum_MPrich.txt cerebral-cortex_MPrich.txt cervix_MPrich.txt colon_MPrich.txt cornua-uteri_MPrich.txt corpus-uteri_MPrich.txt duodenum_MPrich.txt epididymis_MPrich.txt heart_MPrich.txt hippocampus_MPrich.txt hypothalarmus_MPrich.txt ileum_MPrich.txt jejunum_MPrich.txt kidney_MPrich.txt liver_MPrich.txt lung_MPrich.txt lymph-node_MPrich.txt mammary-gland_MPrich.txt medulla-oblongata_MPrich.txt midbrain_MPrich.txt muscle_MPrich.txt omasum_MPrich.txt optic-chiasm_MPrich.txt ovary_MPrich.txt oviduct_MPrich.txt pineal_MPrich.txt pituitary_MPrich.txt pons_MPrich.txt rectum_MPrich.txt reticulum_MPrich.txt rumen_MPrich.txt skin_MPrich.txt soft-horn_MPrich.txt spleen_MPrich.txt splenium_MPrich.txt testis_MPrich.txt thymus_MPrich.txt thyroid_MPrich.txt
    do
    if [ "`grep -w $state $tissue`" ]
    then
      grep -w $state $tissue > 1.txt
      count=$(cat 1.txt | awk -F"\t" '{a = -log($3)/log(10); printf("%0.4f\n",a)}')

      echo $state $tissue $count >> 5.txt
     else
     let count=0
     echo $state $tissue $count >> 5.txt
    fi
  done
done



for tissue in abomasum_MPrich.txt adipose_MPrich.txt bone-marrow_MPrich.txt brainstem_MPrich.txt cecum_MPrich.txt cerebellum_MPrich.txt cerebral-cortex_MPrich.txt cervix_MPrich.txt colon_MPrich.txt cornua-uteri_MPrich.txt corpus-uteri_MPrich.txt duodenum_MPrich.txt epididymis_MPrich.txt heart_MPrich.txt hippocampus_MPrich.txt hypothalarmus_MPrich.txt ileum_MPrich.txt jejunum_MPrich.txt kidney_MPrich.txt liver_MPrich.txt lung_MPrich.txt lymph-node_MPrich.txt mammary-gland_MPrich.txt medulla-oblongata_MPrich.txt midbrain_MPrich.txt muscle_MPrich.txt omasum_MPrich.txt optic-chiasm_MPrich.txt ovary_MPrich.txt oviduct_MPrich.txt pineal_MPrich.txt pituitary_MPrich.txt pons_MPrich.txt rectum_MPrich.txt reticulum_MPrich.txt rumen_MPrich.txt skin_MPrich.txt soft-horn_MPrich.txt spleen_MPrich.txt splenium_MPrich.txt testis_MPrich.txt thymus_MPrich.txt thyroid_MPrich.txt
do
  echo $tissue
  grep -w $tissue 5.txt | cut -d " " -f 3 | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
grep -w oviduct_MPrich.txt 5.txt | cut -d " " -f 1 > 4.txt
echo Go abomasum adipose bone-marrow brainstem cecum cerebellum cerebral-cortex cervix colon cornua-uteri corpus-uteri duodenum epididymis heart hippocampus hypothalarmus ileum jejunum kidney liver lung lymph-node mammary-gland medulla-oblongata midbrain muscle omasum optic-chiasm ovary oviduct pineal pituitary pons rectum reticulum rumen skin soft-horn spleen splenium testis thymus thyroid > 7.txt
paste 4.txt 3.txt | sed 's/ /\t/g' > 6.txt
cat 7.txt 6.txt |  sed 's/ /\t/g'> 9.txt
sed -i 's/inf/100/g' 9.txt
cut -f 1-2 AA_total_Go_orgin_MP.txt  > 8.txt
paste 8.txt 9.txt > TSR_go_MP_enrichment.csv
```

#Plot Figure 2g/Figure 2h/Figure 2i for TSR GO/HP/MP

```python
library(ComplexHeatmap) 
library(circlize)
setwd('/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/GO_GREAT/Gorich/') #mac
#setwd('/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/GO_GREAT/HPrich/') #mac
#setwd('/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_liftsheeptohuman/GO_GREAT/MPrich/') #mac
data <- read.csv('TSR_go_enhancer_enrichment.csv', sep = "\t", header = T) 
#data <- read.csv('TSR_go_HP_enrichment.csv', sep = "\t", header = T) 
#data <- read.csv('TSR_go_MP_enrichment.csv', sep = "\t", header = T) 
colnames(data) <- gsub("\\.", "-", colnames(data))  
data1 = data[,c(4:46)]
rownames(data1)=data[,2] 
data1 <- data1[rowSums(data1)>0,]
split = factor(colnames(data1), levels= c("abomasum", "adipose", "bone-marrow", "brainstem", "cecum", "cerebellum", "cerebral-cortex", "cervix", "colon", "cornua-uteri", "corpus-uteri", "duodenum", "epididymis", "heart", "hippocampus", "hypothalamus", "ileum", "jejunum", "kidney", "liver", "lung", "lymph-node", "mammary-gland", "medulla-oblongata", "midbrain", "muscle", "omasum", "optic-chiasm", "ovary", "oviduct", "pineal", "pituitary", "pons", "rectum", "reticulum", "rumen", "skin", "soft-horn", "spleen", "splenium", "testis", "thymus", "thyroid"))
colnames(data1)= factor(colnames(data1), levels=c("abomasum", "adipose", "bone-marrow", "brainstem", "cecum", "cerebellum", "cerebral-cortex", "cervix", "colon", "cornua-uteri", "corpus-uteri", "duodenum", "epididymis", "heart", "hippocampus", "hypothalamus", "ileum", "jejunum", "kidney", "liver", "lung", "lymph-node", "mammary-gland", "medulla-oblongata", "midbrain", "muscle", "omasum", "optic-chiasm", "ovary", "oviduct", "pineal", "pituitary", "pons", "rectum", "reticulum", "rumen", "skin", "soft-horn", "spleen", "splenium", "testis", "thymus", "thyroid"))

type = colnames(data1)
ha = HeatmapAnnotation(tissue = type,annotation_name_side = "right",annotation_legend_param = list(at = colnames(data1), labels= c("abomasum", "adipose", "bone-marrow", "brainstem", "cecum", "cerebellum", "cerebral-cortex", "cervix", "colon", "cornua-uteri", "corpus-uteri", "duodenum", "epididymis", "heart", "hippocampus", "hypothalamus", "ileum", "jejunum", "kidney", "liver", "lung", "lymph-node", "mammary-gland", "medulla-oblongata", "midbrain", "muscle", "omasum", "optic-chiasm", "ovary", "oviduct", "pineal", "pituitary", "pons", "rectum", "reticulum", "rumen", "skin", "soft-horn", "spleen", "splenium", "testis", "thymus", "thyroid")),
                       col = list(tissue =c("abomasum" = "#f18264",
      "adipose" = "#ffc4f1",
      "bone-marrow" = "#d84b4b",
      "brainstem" = "#f4d578",
      "cecum" = "#daaa6c",
      "cerebellum" = "#efd80b",
      "cerebral-cortex" = "#dcd71a",
      "cervix" = "#b6d7a9",
      "colon" = "#f2c063",
      "cornua-uteri" = "#69d683",
      "corpus-uteri" = "#80d897",
      "duodenum" = "#eb9d63",
      "epididymis" = "#70d24b",
      "heart" = "#bc58e3",
      "hippocampus" = "#f1cb05",
      "hypothalamus" = "#825e19",
      "ileum" = "#ce9639",
      "jejunum" = "#eb951c",
      "kidney" = "#4f3136",
      "liver" = "#ad8c8b",
      "lung" = "#36b5f1",
      "lymph-node" = "#c33a11",
      "mammary-gland" = "#fed9d0",
      "medulla-oblongata" = "#d0b35b",
      "midbrain" = "#f9ed19",
      "muscle" = "#a180ca",
      "omasum" = "#f8ae81",
      "optic-chiasm" = "#fece01",
      "ovary" = "#69d28c",
      "oviduct" = "#79ffaa",
      "pineal" = "#807120",
      "pituitary" = "#f1d95d",
      "pons" = "#fcd222",
      "rectum" = "#efe0a8",
      "reticulum" = "#d97c68",
      "rumen" = "#fc9891",
      "skin" = "#d09dc5",
      "soft-horn" = "#a25d73",
      "spleen" = "#962932",
      "splenium" = "#7c6919",
      "testis" = "#7ef351",
      "thymus" = "#ff3a32",
      "thyroid" = "#f359d1")))
pdf("/vol2/mengzhu/jupyter/figure/TSR_go_GO_enrichment.pdf", width = 12, height = 8)
#pdf("/vol2/mengzhu/jupyter/figure/TSR_go_HP_enrichment.pdf", width = 12, height = 8)
#pdf("/vol2/mengzhu/jupyter/figure/TSR_go_MP_enrichment.pdf", width = 12, height = 8)
Heatmap(data1, border = TRUE, show_column_names = F,show_row_names = F, column_gap = unit(0, "mm"), cluster_column_slices = FALSE, column_title =NULL, column_split =split, bottom_annotation = ha, row_names_gp = gpar(fontsize = 7), cluster_rows = FALSE, cluster_columns = FALSE,col = colorRamp2(c(0, 0, 20),c("white", "white", "red"))) 
dev.off()
```

## TSR MOTIF

### Hormer

```shell
###########Final code#############
###########Use HOMER to identify motifs#############
###########Remove only sequences from the corresponding tissue from the background###########
##/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/Hormer_motif_remove1##
vim 03_TSR_remove1.sh
#!/bin/bash
set -e
id=$1
findMotifsGenome.pl ${id}_id.bed \
/vol2/mengzhu/genome/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna \
Hormer_motif_remove1/${id}_remove1 \
-bg background_$(basename $id "_E5").bed \
-len 8,10,12 -size 200 -mask -p 5
# The script processes one ID supplied as $1; the unmatched source `done` was removed.
#for i in `cat sample.txt`; do sbatch -c 5  -p low --mem 10G 03_TSR_remove1.sh $i; done
```

### Meme fimo

Use MEME FIMO to find motifs in known sequences (including sequence-position information)

```shell

fasta-get-markov -m 1 TSR_adipose_E5_enhancer.fa TSR_adipose_E5_enhancer.bg
fimo --bgfile enhancers.bg --qv-thresh --thresh 0.01 motifs.meme enhancers.fa


#Generate sample.txt
for f in *_E5_id.bed; do
  echo "${f%_id.bed}"
done > sample.txt

########Use MEME FIMO to identify TFBSs##########
##/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5##
vim 02_fimo_motif_p0.001.sh
#!/bin/bash
set -euo pipefail

id="$1"
#fasta-get-markov -m 1 motif/${id}_enhancer.fa motif/${id}_enhancer.bg
motif_file="/vol2/mengzhu/SheepFANNG/03_chrmatin_noblacklist/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E4/motif/jaspar_sheep_core_meme/JASPAR2024_CORE_vertebrates_redundant_pfms_meme.txt"
seq_file="motif/${id}_enhancer.fa"
bg_file="motif/${id}_enhancer.bg"
outdir="motif/fimo_vertebrates_bg0.001/${id}_enhancer"

mkdir -p "$(dirname "$outdir")"

# To use q-value (FDR) as the threshold (recommended for enhancer scans), uncomment the following line:
# Common q-value thresholds are 0.01 or 0.05; use 0.001 for greater stringency
fimo --oc "$outdir" \
     --max-stored-scores 5000000 \
     --thresh 0.001 \
     --bfile "$bg_file" \
     "$motif_file" \
     "$seq_file"
#for i in `cat sample1.txt`; do sbatch -c 40  -p low --mem 80G 02_fimo_motif_p0.001.sh $i; done
sbatch -c 5  -p smp --mem 40G 02_fimo_motif_p0.001.sh TSR_adipose_E5z
```

## TSE

```shell
#Remove rows with a maximum value below 0.1
####/vol2/mengzhu/snakemake_sheep/expressiondir/Average####
awk -F'\t' '
  NR==1 {print; next}
  {
    max = 0
    for (i=7; i<=NF; i++) {
      v = $i + 0
      if (v > max) max = v
    }
    if (max >= 0.1) print
  }' all_tisssues_expression_aveage_tpm.csv > filtered.tsv
cut -f1,7-49 filtered.tsv > expr.tsv

# SPM: one 0-1 score per tissue
vim 02_tspex_spm.sh
#!/bin/bash
tspex --log expr.tsv tspex_spm.tsv spm
#sbatch -p low -c 4 --mem=8G  02_tspex_spm.sh

# Tau: one 0-1 score per gene
vim 03_tspex_tau.sh
#!/bin/bash
tspex expr.tsv tspex_tau.tsv tau
#sbatch -p low -c 4 --mem=8G 03_tspex_tau.sh

# Tsi: one 0-1 score per gene
cat 04_tspex_tsi.sh
#!/bin/bash
tspex --log expr.tsv tspex_tsi.tsv tsi
#sbatch -p low -c 4 --mem=8G 04_tspex_tsi.sh

```

# 4.Network

 

```shell
vim 01_natwork.sh
#!/bin/bash
########1.First identify enhancers
#######2.Then identify enhancer target genes#######
#######3.Then identify enhancer motifs using MEME FIMO#######
#######3.Then add the tissue-specific results#######
#######3.Then link promoters within 2 kb upstream of target genes#######
#######4.Then link TFs#######
#######/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/motif/fimo_vertebrates_bg0.001/natwork#####
cat sample1.txt | while read -r Tissue id Motif
do
    echo "$id"
    echo "$Motif"
    echo "$Tissue"
####Identify enhancers potentially bound by motifs from the FIMO results#####
    fimo="../TSR_${Tissue}_E5_enhancer/fimo.tsv"
    out="${Tissue}_${Motif}_0.bed"

    awk -F'\t' -v OFS='\t' -v id="$id" -v pthr="2.6e-05" '
        BEGIN{
            print "#chrom","start","end","motif_alt_id","score","p-value","q-value","strand"
        }
        FNR>1 && index($1, id)>0 && ($8+0) < pthr {
            print $3, $4, $5, $2, $7, $8, $9, $6
        }
    ' "$fimo" > "$out"
bedtools intersect -wo -a ${Tissue}_${Motif}_0.bed -b /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/TSR_${Tissue}_E5_id.bed  > ${Tissue}_${Motif}_enhancer.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $4,$12,$5,"'"${Tissue}"'"}' ${Tissue}_${Motif}_enhancer.txt > ${Tissue}_${Motif}_TF_Enhn_node.txt
awk -F'\t' 'BEGIN{OFS="\t"} FNR>1{print $12,"Enhancer","'"${Tissue}"'"; print $4,"Motif","'"${Tissue}"'"}' ${Tissue}_${Motif}_enhancer.txt | awk '!seen[$0]++' > ${Tissue}_${Motif}_TF_Enhn_edge.txt
#####Determine motif-enhancer-gene relationships from enhancer-gene pairs########
awk 'NR==FNR { keep[$2]=1; next } ($1 in keep)' \
    ${Tissue}_${Motif}_TF_Enhn_node.txt /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_Target_gene/TSR_${Tissue}_E5_target_gene.txt > TSR_${Tissue}_E5_target_gene.intersect.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$7,$5,"'"${Tissue}"'"}' TSR_${Tissue}_E5_target_gene.intersect.txt > ${Tissue}_${Motif}_Enhn_Gene_node.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,"Enhancer","'"${Tissue}"'"; print $7,"Gene","'"${Tissue}"'"}' TSR_${Tissue}_E5_target_gene.intersect.txt | awk '!seen[$0]++' > ${Tissue}_${Motif}_Enhn_Gene_edge.txt
######Identify gene promoters from the 2-kb upstream regions of genes
awk 'NR==FNR { keep[$7]=1; next } ($4 in keep)' \
    TSR_${Tissue}_E5_target_gene.intersect.txt /vol2/mengzhu/genome/part_change_esemb100/TSS_esemble100_colin.bed_up2k.bed > TSR_${Tissue}_E5_gene_tss2k.txt
bedtools intersect -wo -a TSR_${Tissue}_E5_gene_tss2k.txt -b /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E1/TSR_${Tissue}_E1_id.bed  > TSR_${Tissue}_E1_gene.txt
#awk -F'\t' 'BEGIN{OFS="\t"} {print $11,$4,"1","'"${Tissue}"'"}' TSR_${Tissue}_E1_gene.txt > ${Tissue}_${Motif}_Promo_Gene_node.txt
awk -F'\t' -v OFS='\t' -v tissue="$Tissue" '
  $4 !~ /^LOC/ { print $11, $4, "1", tissue }
' TSR_${Tissue}_E1_gene.txt > ${Tissue}_${Motif}_Promo_Gene_node.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,"Promoter","'"${Tissue}"'"; print $2,"Gene","'"${Tissue}"'"}' ${Tissue}_${Motif}_Promo_Gene_node.txt | awk '!seen[$0]++' > ${Tissue}_${Motif}_Promo_Gene_edge.txt
#####Filter genes in ${Tissue}_${Motif}_Enhn_Gene_node.txt by promoter presence, retaining only genes with promoters
#####Extract enhancers#####
awk -F'\t' '
  NR==FNR {
    if ($2 !~ /^LOC/) keep[$2]=1
    next
  }
  ($2 in keep)
' ${Tissue}_${Motif}_Promo_Gene_node.txt ${Tissue}_${Motif}_Enhn_Gene_node.txt \
> ${Tissue}_${Motif}_Enhn_Gene_node_2.txt
#awk 'NR==FNR { keep[$2]=1; next } ($2 in keep)' ${Tissue}_${Motif}_Promo_Gene_node.txt ${Tissue}_${Motif}_Enhn_Gene_node.txt > ${Tissue}_${Motif}_Enhn_Gene_node_2.txt
awk 'NR==FNR { keep[$1]=1; next } ($1 in keep)' ${Tissue}_${Motif}_Enhn_Gene_node_2.txt ${Tissue}_${Motif}_Enhn_Gene_edge.txt > ${Tissue}_${Motif}_Enhn_Gene_edge_2.txt
#####Filter enhancers in ${Tissue}_${Motif}_TF_Enhn_node.txt by gene links, retaining only enhancers linked to genes
#####Extract motifs#####
awk 'NR==FNR { keep[$1]=1; next } ($2 in keep)' ${Tissue}_${Motif}_Enhn_Gene_node_2.txt ${Tissue}_${Motif}_TF_Enhn_node.txt > ${Tissue}_${Motif}_TF_Enhn_node_2.txt
awk 'NR==FNR { keep[$1]=1; next } ($1 in keep)' ${Tissue}_${Motif}_TF_Enhn_node_2.txt ${Tissue}_${Motif}_TF_Enhn_edge.txt > ${Tissue}_${Motif}_TF_Enhn_edge_2.txt
#####Merge to generate node and edge files######
cat ${Tissue}_${Motif}_TF_Enhn_node_2.txt ${Tissue}_${Motif}_Enhn_Gene_node_2.txt ${Tissue}_${Motif}_Promo_Gene_node.txt > node/${Tissue}_${Motif}_node.txt
cat ${Tissue}_${Motif}_TF_Enhn_edge_2.txt ${Tissue}_${Motif}_Enhn_Gene_edge_2.txt ${Tissue}_${Motif}_Promo_Gene_edge.txt > node/${Tissue}_${Motif}_edge.txt
done
#sbatch -c 50  -p low --mem 100G 01_natwork.sh
```

```shell
cat 02_data_statistics.sh
#!/bin/bash
ls *_TF_Enhn_node_2.txt |cut -d '_' -f 1-2 | while read sample;
do
echo ${sample}
cat ${sample}_Enhn_Gene_node.txt | cut -f 1 | awk '!seen[$1]++' > ${sample}_enhancer_region.txt
cat ${sample}_Enhn_Gene_node.txt | cut -f 2 | awk '!seen[$1]++' > ${sample}_enhancer_gene.txt
all=$(cat ${sample}_node.txt |  wc -l | awk '{print $1}')
promoter=$(cat ${sample}_Promo_Gene_node.txt | cut -f 1 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
promoterGene=$(cat ${sample}_Promo_Gene_node.txt | cut -f 2 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
enhancer=$(cat ${sample}_Enhn_Gene_node.txt | cut -f 1 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
enhancerGene=$(cat ${sample}_Enhn_Gene_node.txt | cut -f 2 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
rm ${sample}_region_count.txt
cat ${sample}_enhancer_region.txt | while read id;
do
Reads1=$(grep -w $id ${sample}_Enhn_Gene_node.txt | wc -l | awk '{print $1}')
echo $id $Reads1 >> ${sample}_region_count.txt
done
rm ${sample}_gene_count.txt
cat ${sample}_enhancer_gene.txt | while read id;
do
Reads2=$(grep -w $id ${sample}_Enhn_Gene_node.txt | wc -l | awk '{print $1}')
echo $id $Reads2 >> ${sample}_gene_count.txt
done

  read enhancer_avg enhancer_med < <(
    cut -d' ' -f2 ${sample}_region_count.txt \
    | awk '/^[0-9]+(\.[0-9]+)?$/' \
    | sort -n \
    | awk '{x[NR]=$1; s+=$1}
           END{
              if (NR==0){print "NA\tNA"; exit}
              if (NR%2){m=x[(NR+1)/2]} else {m=(x[NR/2]+x[NR/2+1])/2}
              print s/NR "\t" m
           }'
  )

  read Gene_avg Gene_med < <(
    cut -d' ' -f2 ${sample}_gene_count.txt \
    | awk '/^[0-9]+(\.[0-9]+)?$/' \
    | sort -n \
    | awk '{x[NR]=$1; s+=$1}
           END{
              if (NR==0){print "NA\tNA"; exit}
              if (NR%2){m=x[(NR+1)/2]} else {m=(x[NR/2]+x[NR/2+1])/2}
              print s/NR "\t" m
           }'
  )

  # Write the current-sample summary (adding four columns: region_mean, region_median, gene_mean, and gene_median)
  echo "$all" "$sample" "$enhancer" "$enhancerGene" "$enhancer_avg" "$enhancer_med" "$Gene_avg" "$Gene_med" "$promoter" "$promoterGene" \
    > ${sample}_summary.txt

  # Append to the combined table
  cat ${sample}_summary.txt >> all_target_gene_summary.txt
done

# Convert to tab-delimited format
sed 's/ /\t/g' all_target_gene_summary.txt > all_target_gene_summary.csv
#sbatch -c 50  -p low --mem 100G 02_data_statistics.sh

```

# 5.Selection signatures



```shell
#Modern Asian versus modern European SNPs: /storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/06.population/01.admixture/v1.all/vcf/chrAuto.vcf.gz
#Ancient Asian and European DNA: /storage/public/home/2021050411/20.phase/v3/info0.80/

#ancient asia 2222
#/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/fst.ancient_CEA/ancientasia_4000y.10000_10000.windowed.weir.fst.filter.t0.01.annotation

#ancient europe 2330
#/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/fst.ancient_EUR/ancienteurope_4000y.10000_10000.windowed.weir.fst.filter.t0.01.annotation

#moderneurope_modernasia 2134
#/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/fst_CEA_EUR/moderneurope_modernasia.10000_10000.windowed.weir.fst.filter.t0.01.annotation

#Demestication 300
#/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/Demestication_domestic/Demestication_domestic_sheep.bed

#Locations of tissue-specific regulatory elements
#/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E1/TSR_thyroid_E1_id.bed

#Merge chr*.1000.maf0.01.filter.vcf.gz files
vim vcf_concat.sh
#!/bin/bash
ls -1v chr*.1000.maf0.01.filter.vcf.gz > vcflist
bcftools concat --file-list vcflist -Oz -o chrAuto_ancient.vcf.gz --threads 8
bcftools index --threads 8 chrAuto_ancient.vcf.gz
#sbatch -c 8  -p low --mem 16G vcf_concat.sh

```

#Calculate the proportion of selected regions in noncoding regions

```shell
#!/bin/bash
set -euo pipefail

# =========================
# 1. Input file
# =========================
fst_file="Demestication_domestic/Demestication_domestic_sheep.bed"
cds_file="/vol2/mengzhu/genome/part_change_esemb100/CDS_esemble100_colin.bed"

# Output directory
outdir="fst_cds_stat"
mkdir -p "${outdir}"

# =========================
# 2. Convert the FST window file to BED
#    Original file:
#    CHROM BIN_START BIN_END ...
#    BED requires:
#    chrom start end
#
#    Assumptions:
#    BIN_START/BIN_END are commonly used interval boundaries,
#    convert to BED with start = BIN_START - 1 and end = BIN_END
# =========================
#awk 'BEGIN{FS=OFS="\t"}
#NR==1 {next}
#{
#    start = $2 - 1
#    if (start < 0) start = 0
    # Retain the original information for downstream tracking
#    print $1, start, $3, $1 ":" $2 "-" $3, $2, $3, $4, $5, $6, $7, $8
#}' "${fst_file}" \
#| sort -k1,1 -k2,2n -k3,3n \
#> "${outdir}/Demestication_domestic_windows.sorted.bed"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {next}
NF>=3 && $1!="" && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ {
    start = $2 - 1
    if (start < 0) start = 0
    print $1, start, $3
}' "${fst_file}" \
| sed 's/\r$//' \
| sort -k1,1 -k2,2n -k3,3n \
> "${outdir}/Demestication_domestic_windows.sorted.bed"

# Field descriptions:
# col1  chrom
# col2  bed_start
# col3  bed_end
# col4  window_id
# col5  original_BIN_START
# col6  original_BIN_END
# col7  N_VARIANTS
# col8  WEIGHTED_FST
# col9  MEAN_FST
# col10 type
# col11 ann

# =========================
# 3. Preprocess the CDS file
#    Keep only the first three columns: chr/start/end
#    Then sort and merge to prevent duplicated CDS segments from affecting the statistics
# =========================
cut -f1-3 "${cds_file}" \
| awk 'BEGIN{FS=OFS="\t"} $2 < $3 {print $1, $2, $3}' \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> "${outdir}/cds.merged.bed"

# =========================
# 4. Count windows that overlap CDS regions
#    -u: Output a window from A once if it has any overlap with B
# =========================
bedtools intersect \
    -a "${outdir}/Demestication_domestic_windows.sorted.bed" \
    -b "${outdir}/cds.merged.bed" \
    -sorted \
    -u \
> "${outdir}/Demestication_domestic_windows.in_CDS.bed"

# =========================
# 5. Count windows with no CDS overlap
#    -v: Output only windows from A that have no overlap at all
# =========================
bedtools intersect \
    -a "${outdir}/Demestication_domestic_windows.sorted.bed" \
    -b "${outdir}/cds.merged.bed" \
    -sorted \
    -v \
> "${outdir}/Demestication_domestic_windows.in_nonCDS.bed"

# =========================
# 6. Label each window as CDS or nonCDS
#    -c: Count the number of CDS overlaps for each window
# =========================
bedtools intersect \
    -a "${outdir}/Demestication_domestic_windows.sorted.bed" \
    -b "${outdir}/cds.merged.bed" \
    -sorted \
    -c \
| awk 'BEGIN{FS=OFS="\t"; print "chrom","bed_start","bed_end","window_id","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST","type","ann","cds_overlap_count","class"}
{
    cls = ($NF > 0 ? "CDS" : "nonCDS")
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,cls
}' \
> "${outdir}/Demestication_domestic_windows.CDS_classification.tsv"

# =========================
# 7. Calculate totals and proportions
# =========================
total_windows=$(wc -l < "${outdir}/Demestication_domestic_windows.sorted.bed")
cds_windows=$(wc -l < "${outdir}/Demestication_domestic_windows.in_CDS.bed")
noncds_windows=$(wc -l < "${outdir}/Demestication_domestic_windows.in_nonCDS.bed")

cds_ratio=$(awk -v a="${cds_windows}" -v b="${total_windows}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')
noncds_ratio=$(awk -v a="${noncds_windows}" -v b="${total_windows}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')

# Percentage format
cds_pct=$(awk -v x="${cds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')
noncds_pct=$(awk -v x="${noncds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')

# Write the summary table
{
    echo -e "category\tcount\tratio\tpercent"
    echo -e "CDS\t${cds_windows}\t${cds_ratio}\t${cds_pct}"
    echo -e "nonCDS\t${noncds_windows}\t${noncds_ratio}\t${noncds_pct}"
    echo -e "total\t${total_windows}\t1.000000\t100.00%"
} > "${outdir}/Demestication_domestic_summary_CDS_vs_nonCDS.tsv"

# Also print to the screen
echo "===== Summary ====="
cat "${outdir}/Demestication_domestic_summary_CDS_vs_nonCDS.tsv"

echo
echo "Output files:"
echo "${outdir}/Demestication_domestic_windows.sorted.bed"
echo "${outdir}/cds.merged.bed"
echo "${outdir}/Demestication_domestic_windows.in_CDS.bed"
echo "${outdir}/Demestication_domestic_windows.in_nonCDS.bed"
echo "${outdir}/Demestication_domestic_windows.CDS_classification.tsv"
echo "${outdir}/Demestication_domestic_summary_CDS_vs_nonCDS.tsv"

cut -f1-3 "Demestication_domestic_windows.sorted.bed" \
| sed 's/\r$//' \
| awk 'BEGIN{FS=OFS="\t"} NF>=3 && $2 < $3 {print $1,$2,$3}' \
| sort -k1,1 -k2,2n -k3,3n \
> "Demestication_domestic_windows_3col.sorted.bed"

bedtools coverage \
    -a "Demestication_domestic_windows_3col.sorted.bed" \
    -b "cds.merged.bed" \
> "Demestication_domestic.coverage_by_CDS.tsv"

read total_bp cds_bp < <(
    awk 'BEGIN{total=0; cds=0}
    {
        # coverage output:
        # col1 chr
        # col2 start
        # col3 end
        # col4 number of overlapping B features
        # col5 number of covered bp in A
        # col6 length of the A interval
        # col7 coverage fraction
        cds += $5
        total += $6
    }
    END{
        print total, cds
    }' "Demestication_domestic.coverage_by_CDS.tsv"
)

noncds_bp=$(( total_bp - cds_bp ))

cds_ratio=$(awk -v a="${cds_bp}" -v b="${total_bp}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')
noncds_ratio=$(awk -v a="${noncds_bp}" -v b="${total_bp}" 'BEGIN{if(b==0){printf "%.6f",0}else{printf "%.6f",a/b}}')

cds_pct=$(awk -v x="${cds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')
noncds_pct=$(awk -v x="${noncds_ratio}" 'BEGIN{printf "%.2f%%", x*100}')

# =========================
# 7. Write results
# =========================
{
    echo -e "category\tbp\tratio\tpercent"
    echo -e "CDS\t${cds_bp}\t${cds_ratio}\t${cds_pct}"
    echo -e "nonCDS\t${noncds_bp}\t${noncds_ratio}\t${noncds_pct}"
    echo -e "total_selected_region\t${total_bp}\t1.000000\t100.00%"
} > "Demestication_domestic_selected_region_CDS_percent.tsv"


for f in *_selected_region_CDS_percent.tsv; do
    sample=${f%_selected_region_CDS_percent.tsv}
    cds_pct=$(awk '$1=="CDS" {print $4}' "$f")
    noncds_pct=$(awk '$1=="nonCDS" {print $4}' "$f")
    echo -e "${sample}\t${cds_pct}\t${noncds_pct}"
done > all_selected_region_CDS_percent_summary.tsv

```

### Top 10% fold enrichment

Install: https://bioconda.github.io/recipes/gat/README.html

```sh
# 1) Create the environment and install dependencies (mamba is faster; conda works as well)
mamba create -n gat -c conda-forge -c bioconda python=2.7 gat=1.3.6
# 2) Activate the environment
conda activate gat
# 3) Verify the installation
gat-run.py --help
```

Run: https://gat.readthedocs.io/en/latest/tutorialIntervalOverlap.html

```shell
####0. Generate the genome-size file####
####/vol2/mengzhu/genome/reference####
awk 'BEGIN{OFS="\t"}
{
  chr=$1
  if(chr !~ /^chr/) chr="chr"chr
  print chr,0,$2,"ws"
}' sheep1.size \
| sort -k1,1 -k2,2n \
| gzip -c > sheep_workspace_allchr.bed.gz

####1. Extract the top 10% of selection-signal regions####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature####
in="moderneurope_modernasia.10000_10000.windowed.weir.fst.filter.t0.01.annotation"
out="fst.europe_asia_10000.decile.top10pct.tsv"

n=$(( $(wc -l < "$in") - 1 ))      # Number of data rows (excluding the header)
k=$(( (n + 9) / 10 ))              # Top 10%: ceil(n/10)
(( k < 1 )) && k=1                 # Keep at least one row (for very small files)

{
  head -n 1 "$in"
  tail -n +2 "$in" | sort -t $'\t' -k5,5gr | head -n "$k"
} > "$out"

####2. Generate fst.ancient_CEA.top10pct.bed.gz (process selection signatures)####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature####
awk 'BEGIN{FS=OFS="\t"}
NR==1{next}
{
  chr=$1;
  if(chr !~ /^chr/) chr="chr"chr;   # 13 -> chr13 (leave unchanged if already chr13)
  start=$2-1;                       # 1-based -> BED 0-based
  end=$3;
  if(start<0) start=0;
  print chr,start,end,"selection"   # Use one track name in column 4 so intervals are not treated as separate tracks
}' ancientasia_4000y.10000_10000.windowed.weir.fst.filter.t0.01.annotation \
| sort -k1,1 -k2,2n \
| gzip -c > GAT_fst.ancientasia_10000.decile.top10pct_1.bed.gz

####3. Generate GAT_TSR_*_E5.bed.gz files (process TSR enhancers)####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5####
outdir="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature"
mkdir -p "$outdir"

shopt -s nullglob  # If no files match, do not pass the wildcard literally

for inbed in TSR_*_E5_id.bed; do
  state="${inbed%_id.bed}"   # For example: TSR_thymus_E5
  out="$outdir/GAT_${state}.bed.gz"

  awk -v OFS='\t' -v s="$state" '{print $1,$2,$3,s}' "$inbed" \
  | sort -k1,1 -k2,2n \
  | gzip -c > "$out"

  echo "[OK] $inbed -> $out"
done

####4. Run####
####Use GAT to calculate fold enrichment between selection signals and tissue-specific enhancers
vim 02_GAT_europe_asia.sh
#!/bin/bash
set -euo pipefail

SEG="GAT_fst.europe_asia_10000.decile.top10pct.bed.gz"
WS="/vol2/mengzhu/genome/reference/sheep_workspace_allchr.bed.gz"

OUTDIR="GAT_results_ancient_10000_europe_asia"
LOGDIR="GAT_logs_ancient_10000_europe_asia"
mkdir -p "$OUTDIR" "$LOGDIR"

shopt -s nullglob

for ANN in GAT_*.bed.gz; do
  # Skip the segment file itself (it also matches GAT_*.bed.gz)
  [[ "$ANN" == "$SEG" ]] && continue

  base="$(basename "$ANN" .bed.gz)"   # e.g. GAT_TSR_thymus_E5
  tag="${base#GAT_}"                  # e.g. TSR_thymus_E5

  echo "[RUN] annotation=$ANN -> ${OUTDIR}/${tag}.tsv"

  gat-run.py \
    --segment-file="$SEG" \
    --annotation-file="$ANN" \
    --workspace-file="$WS" \
    --counter=nucleotide-overlap \
    --num-samples=10000 \
    --qvalue-method=BH \
    --num-threads=10 \
    --nbuckets=100001 \
    --log="${LOGDIR}/${tag}.log" \
    > "${OUTDIR}/${tag}.europe_asia.tsv"
done
#sbatch -c 10  -p low --mem 20G 02_GAT_europe_asia.sh

###Domestication analysis###
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "selection"}' Demestication_domestic_sheep.bed > Demestication_domestic_sheep.selection.bed
gzip -c Demestication_domestic_sheep.selection.bed >Demestication_domestic_sheep.selection.bed.gz

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "selection"}' Demestication_mouflon.bed > Demestication_mouflon_sheep.selection.bed
gzip -c Demestication_mouflon_sheep.selection.bed >Demestication_mouflon_sheep.selection.bed.gz

vim 02_GAT_Demestication_mouflon.sh
#!/bin/bash
set -euo pipefail

SEG="Demestication_mouflon_sheep.selection.bed.gz"
WS="/vol2/mengzhu/genome/reference/sheep_workspace_allchr.bed.gz"

OUTDIR="GAT_results_Demestication_mouflon"
LOGDIR="GAT_logs_Demestication_mouflon"
mkdir -p "$OUTDIR" "$LOGDIR"

shopt -s nullglob

for ANN in GAT_*.bed.gz; do
  # Skip the segment file itself (it also matches GAT_*.bed.gz)
  [[ "$ANN" == "$SEG" ]] && continue

  base="$(basename "$ANN" .bed.gz)"   # e.g. GAT_TSR_thymus_E5
  tag="${base#GAT_}"                  # e.g. TSR_thymus_E5

  echo "[RUN] annotation=$ANN -> ${OUTDIR}/${tag}.tsv"

  gat-run.py \
    --segment-file="$SEG" \
    --annotation-file="$ANN" \
    --workspace-file="$WS" \
    --counter=nucleotide-overlap \
    --num-samples=10000 \
    --qvalue-method=BH \
    --num-threads=10 \
    --nbuckets=100001 \
    --log="${LOGDIR}/${tag}.log" \
    > "${OUTDIR}/${tag}.Demestication_mouflon.tsv"
done
#sbatch -c 10  -p low --mem 20G 02_GAT_Demestication_mouflon.sh

#Prepare data and generate the plotting input file
awk 'FNR>1 || NR==1' TSR_*_E5.ancient_CEA_1.tsv > All_TSR_E5.ancient_CEA_1.tsv
awk 'FNR>1 || NR==1' TSR_*_E5.ancient_EUR_1.tsv > All_TSR_E5.ancient_EUR_1.tsv
awk 'FNR>1 || NR==1' TSR_*_E5.europe_asia_1.tsv > All_TSR_E5.europe_asia_1.tsv
awk 'FNR>1 || NR==1' TSR_*.Demestication_domestic.tsv > All_TSR_E5.Demestication_domestic.tsv
awk 'BEGIN{
  OFS="\t";
  print "Demestication_domestic_annotation","Demestication_domestic_fold","Demestication_domestic_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.Demestication_domestic.tsv > All_TSR_E5.Demestication_domestic_1.tsv

awk 'FNR>1 || NR==1' TSR_*.Demestication_mouflon.tsv > All_TSR_E5.Demestication_mouflon.tsv
awk 'BEGIN{
  OFS="\t";
  print "Demestication_mouflon_annotation","Demestication_mouflon_fold","Demestication_mouflon_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.Demestication_mouflon.tsv > All_TSR_E5.Demestication_mouflon_1.tsv

awk 'BEGIN{
  OFS="\t";
  print "ancient_EUR_1000_annotation","ancient_EUR_1000_fold","ancient_EUR_1000_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.ancient_EUR_1.tsv > All_TSR_E5.ancient_EUR_2.tsv

awk 'BEGIN{
  OFS="\t";
  print "europe_asia_1000_annotation","europe_asia_1000_fold","europe_asia_1000_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.europe_asia_1.tsv > All_TSR_E5.europe_asia_2.tsv

awk 'BEGIN{
  OFS="\t";
  print "ancient_CEA_1000_annotation","ancient_CEA_1000_fold","ancient_CEA_1000_qvalue"
}
NR>1{
  print $2,$8,$11
}' All_TSR_E5.ancient_CEA_1.tsv > All_TSR_E5.ancient_CEA_2.tsv

paste \
  All_TSR_E5.ancient_CEA_2.tsv \
  <(cut -f2-3 All_TSR_E5.ancient_EUR_2.tsv) \
  <(cut -f2-3 All_TSR_E5.europe_asia_2.tsv) \
  <(cut -f2-3 All_TSR_E5.Demestication_domestic_1.tsv) \
> All_TSR_E5_Selection_signatures_1.tsv

```

#Code for plotting Figure 7c

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam

# =========================
# 1. Input files and parameters
# =========================
vcf_file = "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/vcf/chrAuto_ancient.vcf.gz"
fst_file = "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/vcf/chrAuto_neweurope_oldeurope_point.windowed.weir.fst"
region = "chr8:90340000-90380000"
out_prefix = "/vol2/mengzhu/jupyter/figure/Selection_signature"

# Smoothing window: rolling mean by number of SNPs
smooth_window = 15

# Label the positions with the highest FST values
top_n = 5


# =========================
# 2. Utility functions
# =========================
def normalize_chrom(chrom: str) -> str:
    """
    Normalize chromosome names by removing the chr prefix
    chr8 -> 8
    8    -> 8
    """
    chrom = str(chrom).strip()
    chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    return chrom


def parse_region(region_str: str):
    """
    Parse a region string
    Example: chr8:90340000-90380000
    Return: chrom, start, end
    """
    region_str = region_str.replace(",", "").strip()
    chrom_part, pos_part = region_str.split(":")
    start_str, end_str = pos_part.split("-")
    chrom = normalize_chrom(chrom_part)
    start = int(start_str)
    end = int(end_str)
    return chrom, start, end


def safe_float(x):
    """
    Safely convert a value to a floating-point number
    """
    if x is None:
        return np.nan
    if isinstance(x, (tuple, list)):
        if len(x) == 0:
            return np.nan
        x = x[0]
    try:
        return float(x)
    except Exception:
        return np.nan


# =========================
# 3. Read regional FST data
# =========================
def read_fst_region(fst_path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """
    Read the FST file and filter the target region
    File format:
    CHROM   POS     WEIR_AND_COCKERHAM_FST
    """
    df = pd.read_csv(
        fst_path,
        sep=r"\s+",
        engine="python",
        dtype={"CHROM": str}
    )

    required_cols = ["CHROM", "POS", "WEIR_AND_COCKERHAM_FST"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"The FST file is missing column: {col}")

    df["CHROM_norm"] = df["CHROM"].map(normalize_chrom)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["FST"] = pd.to_numeric(df["WEIR_AND_COCKERHAM_FST"], errors="coerce")

    region_df = df[
        (df["CHROM_norm"] == chrom) &
        (df["POS"] >= start) &
        (df["POS"] <= end)
    ][["CHROM_norm", "POS", "FST"]].copy()

    region_df.columns = ["CHROM", "POS", "FST"]
    region_df = region_df.sort_values("POS").reset_index(drop=True)
    return region_df


# =========================
# 4. Read regional VCF data using the .csi index
# =========================
def read_vcf_region(vcf_path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """
    Read the VCF by region with pysam
    Support both contig naming conventions: 8 and chr8
    """
    rows = []
    vcf = pysam.VariantFile(vcf_path)

    # Read contig names from the VCF header to determine whether names use 8 or chr8
    contigs_in_vcf = set(vcf.header.contigs.keys())

    candidates = [chrom, f"chr{chrom}"]
    fetch_contig = None
    for c in candidates:
        if c in contigs_in_vcf:
            fetch_contig = c
            break

    if fetch_contig is None:
        raise ValueError(
            f"The VCF does not contain chromosome {chrom} or chr{chrom}."
            f" Example contigs from the header: {list(contigs_in_vcf)[:10]}"
        )

    # pysam fetch: start is 0-based and stop defines a half-open interval
    for rec in vcf.fetch(fetch_contig, start - 1, end):
        rows.append({
            "CHROM": normalize_chrom(rec.chrom),
            "POS": rec.pos,
            "REF": rec.ref,
            "ALT": ",".join(rec.alts) if rec.alts else ".",
            "AF": safe_float(rec.info.get("AF")),
            "RAF": safe_float(rec.info.get("RAF")),
            "INFO_SCORE": safe_float(rec.info.get("INFO")),
            "MAF": safe_float(rec.info.get("MAF"))
        })

    vcf.close()

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("POS").reset_index(drop=True)
    return df


# =========================
# 5. Plotting
# =========================
def plot_region(df: pd.DataFrame, region_str: str, out_prefix: str,
                smooth_window: int = 15, top_n: int = 5):
    """
    Three-panel plot:
    1) FST scatter plot and smoothed curve
    2) INFO scatter plot
    3) MAF scatter plot
    """
    if df.empty:
        raise ValueError("No data are available for plotting.")

    df = df.sort_values("POS").copy()

    win = min(smooth_window, max(3, len(df)))
    df["FST_SMOOTH"] = df["FST"].rolling(
        window=win, center=True, min_periods=1
    ).mean()

    label_df = df.dropna(subset=["FST"]).nlargest(min(top_n, len(df)), "FST").copy()

    fig, axes = plt.subplots(
        3, 1,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3.5, 1.2, 1.2]}
    )

    ax1, ax2, ax3 = axes

    # ---------- panel 1: FST ----------
    ax1.scatter(df["POS"], df["FST"], s=22, alpha=0.8, label="Per-SNP FST")
    ax1.plot(df["POS"], df["FST_SMOOTH"], linewidth=2, label=f"Rolling mean (n={win})")
    ax1.axhline(0, linestyle="--", linewidth=1)
    ax1.set_ylabel("Weir & Cockerham FST")
    ax1.set_title(f"Local fine-scale selection signal: {region_str}")

    for _, row in label_df.iterrows():
        ax1.scatter(row["POS"], row["FST"], s=40, zorder=3)
        ax1.annotate(
            f"{int(row['POS'])}\n{row['FST']:.3f}",
            xy=(row["POS"], row["FST"]),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8
        )

    ax1.legend(frameon=False, loc="upper right")

    # ---------- panel 2: INFO ----------
    if "INFO_SCORE" in df.columns and df["INFO_SCORE"].notna().any():
        ax2.scatter(df["POS"], df["INFO_SCORE"], s=18, alpha=0.8)
        ax2.set_ylabel("INFO")
        ymin = max(0, np.nanmin(df["INFO_SCORE"]) - 0.05)
        ymax = min(1.05, np.nanmax(df["INFO_SCORE"]) + 0.05)
        if np.isfinite(ymin) and np.isfinite(ymax) and ymin < ymax:
            ax2.set_ylim(ymin, ymax)
    else:
        ax2.text(0.5, 0.5, "No INFO field", transform=ax2.transAxes,
                 ha="center", va="center")
        ax2.set_ylabel("INFO")

    # ---------- panel 3: MAF ----------
    if "MAF" in df.columns and df["MAF"].notna().any():
        ax3.scatter(df["POS"], df["MAF"], s=18, alpha=0.8)
        ax3.set_ylabel("MAF")
        ax3.set_ylim(0, 0.55)
    else:
        ax3.text(0.5, 0.5, "No MAF field", transform=ax3.transAxes,
                 ha="center", va="center")
        ax3.set_ylabel("MAF")

    chrom, start, end = parse_region(region_str)
    ax3.set_xlim(start, end)
    ax3.set_xlabel("Genomic position")

    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    plt.close()


# =========================
# 6. Main workflow
# =========================
def main():
    chrom, start, end = parse_region(region)

    print("Reading FST...")
    fst_df = read_fst_region(fst_file, chrom, start, end)
    if fst_df.empty:
        raise ValueError(f"Region  {region}  has no sites in the FST file.")

    print("Reading VCF with CSI index...")
    vcf_df = read_vcf_region(vcf_file, chrom, start, end)

    print("Merging...")
    merged = fst_df.merge(vcf_df, on=["CHROM", "POS"], how="left")
    merged = merged.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    # Save the regional data table
    merged.to_csv(f"{out_prefix}.region_data.tsv", sep="\t", index=False)

    print("Plotting...")
    plot_region(
        df=merged,
        region_str=region,
        out_prefix=out_prefix,
        smooth_window=smooth_window,
        top_n=top_n
    )

    print("Done!")
    print(f"Output table : {out_prefix}.region_data.tsv")
    print(f"Output figure: {out_prefix}.png")
    print(f"Output figure: {out_prefix}.pdf")


if __name__ == "__main__":
    main()
```

### Calculate point FST

```sh
####Calculate neweurope_oldeurope_point####
####Calculate newasia_oldasia_point####
vim 01.fst_point_1.sh
#!/bin/bash
# run fst between two populations
mkdir -p neweurope_newasia_point neweurope_oldasia_point newasia_oldeurope_point oldasia_oldeurope_point
for chr in chr{1..26} 
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J neweurope_oldeurope_${chr}_fst.point -e neweurope_oldeurope_point/${chr}.point.%J.log -o neweurope_oldeurope_point/${chr}.point.%J.log "bash fst_point.sh $chr neweurope oldeurope"
         jsub -q normal -n 1 -R "span[hosts=1]" -J neweurope_newasia_${chr}_fst.point -e neweurope_newasia_point/${chr}.point.%J.log -o neweurope_newasia_point/${chr}.point.%J.log "bash fst_point.sh $chr neweurope newasia"
         jsub -q normal -n 1 -R "span[hosts=1]" -J neweurope_oldasiae_${chr}_fst.point -e neweurope_oldasia_point/${chr}.point.%J.log -o neweurope_oldasia_point/${chr}.point.%J.log "bash fst_point.sh $chr neweurope oldasia"
        jsub -q normal -n 1 -R "span[hosts=1]" -J newasia_oldasia_${chr}_fst.point -e newasia_oldasia_point/${chr}.point.%J.log -o newasia_oldasia_point/${chr}.point.%J.log "bash fst_point.sh $chr newasia oldasia"
        jsub -q normal -n 1 -R "span[hosts=1]" -J newasia_oldeurope_${chr}_fst.point -e newasia_oldeurope_point/${chr}.point.%J.log -o newasia_oldeurope_point/${chr}.point.%J.log "bash fst_point.sh $chr newasia oldeurope"
        jsub -q normal -n 1 -R "span[hosts=1]" -J oldasia_oldeurope_${chr}_fst.point -e oldasia_oldeurope_point/${chr}.point.%J.log -o oldasia_oldeurope_point/${chr}.point.%J.log "bash fst_point.sh $chr oldasia oldeurope"
        
done

# combine and filter by top 1%
bash fst.combine_point.sh neweurope oldeurope
bash fst.combine_point.sh newasia oldasia

cat fst.combine_point.sh
#!/bin/bash
pop1=$1
pop2=$2
win=$3
step=$4

for chr in chr{1..26} ; do sed '1d' ${pop1}_${pop2}/${chr}.${win}_${step}.windowed.weir.fst | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t%.6f\t%.6f\n",$5,$6}' ; done | sed '1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST' > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst
awk '$4>=40' ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter

# filter and extract top 1% region
lines=`cat ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter | wc -l`
top_lines=$(($lines/100+1))
csvtk sort -t -k 5:Nr ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter | head -$top_lines > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01
```

```sh
####Calculate eur_cea_point####
vim 01.fst_point.sh
#!/bin/bash
# run fst between two populations
mkdir -p eur_cea_point
ln -s eur_cea cea_eur
for chr in chr{1..26}
do

        jsub -q normal -n 1 -R "span[hosts=1]" -J eur_cea_${chr}_fst_point -e eur_cea_point/${chr}_point.%J.log -o eur_cea_point/${chr}_point.%J.log "bash fst_point.sh $chr eur cea"
done

# combine and filter by top 1%
bash fst.combine.sh eur cea 50000 10000
bash fst.combine.sh eur cea 10000 10000
bash fst.combine.sh eur mou 50000 10000
bash fst.combine.sh cea mou 50000 10000
```

```sh
vim 01.merge.sh
#!/bin/bash
{ 
    head -n 1 chr1.point.weir.fst
    awk 'FNR > 1' chr*.point.weir.fst | sort -k1,1n -k2,2n
} > chrAuto_${1}.windowed.weir.fst
#jsub -q normal -n 8 -o output.%J -e error.%J bash 01.merge.sh newasia_oldeurope_point
```



```shell
#
```



```shell
#
```


























































































































