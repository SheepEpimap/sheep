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
