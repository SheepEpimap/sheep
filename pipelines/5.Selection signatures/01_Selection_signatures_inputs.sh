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

