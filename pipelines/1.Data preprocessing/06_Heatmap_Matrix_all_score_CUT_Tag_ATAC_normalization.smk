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
