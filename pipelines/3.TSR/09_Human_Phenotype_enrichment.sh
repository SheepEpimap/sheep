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

