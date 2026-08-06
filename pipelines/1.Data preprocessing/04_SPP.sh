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

