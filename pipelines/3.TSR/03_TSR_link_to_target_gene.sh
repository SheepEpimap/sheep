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
