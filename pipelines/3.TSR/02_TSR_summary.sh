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
