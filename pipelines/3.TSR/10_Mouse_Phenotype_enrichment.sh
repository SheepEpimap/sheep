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
