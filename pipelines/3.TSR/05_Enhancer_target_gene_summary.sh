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
