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


