#Split by chromatin state
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep1/LearnModel_11/new_LearnModel_11
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/LearnModel_11/new_LearnModel_11
ls *_11_segments.bed | while read id;
do
  echo $id
  b=$(basename $id "_11_segments.bed")
  echo $b
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11
  do
     grep -w $state $id | sort -k1,1 -k2,2n > "state_variability/"$b"_"$state".bed"
  done
done

#Merge Rep1 and Rep2
vim 01_merge.sh
#!/usr/bin/env bash
while read -r id; do
  echo "$id"
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11; do   # ← do is present
    cat \
      "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep1/LearnModel_11/new_LearnModel_11/state_variability/${id}_39_${state}.bed" \
      "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Rep2/LearnModel_11/new_LearnModel_11/state_variability/${id}_40_${state}.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -c 4 -o distinct \
      > "${id}_${state}.bed"
  done
done < sample.txt
#sbatch -c 5  -p low --mem 10G 01_merge.sh

#Gs
####get Gs: total regions of each state (Gs) across 14 tissues 
####Merge all tissues by chromatin state
vim 02_AAGs.sh
#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability
mkdir AAGs
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11
do
echo $i
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done
#sbatch -c 5  -p low --mem 10G 02_AAGs.sh

#Identify tissue-specific chromatin states
#Intersect each chromatin state from each tissue with the merged chromatin states; write the number of overlaps in column 4, or 0 when there is no overlap
vim 03_AARegulatory_module.sh
#!/bin/bash
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability
mkdir AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11
do
ls *${state}.bed | while read id;
do
echo $id
bedtools intersect -a <(sort -k1,1 -k2,2n /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs/${state}_Gs.bed) -b <(sort -k1,1 -k2,2n ${id}) -c -sorted >  "AARegulatory_module/"${id%%.*}"_Gs.bed" #-coutputs each record from A and appends a column containing the number of overlapping records in B (an integer; 0 indicates no overlap
done
done
#sbatch -c 5  -p low --mem 10G 03_AARegulatory_module.sh

cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability
rm 1.txt
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
echo chr start end > 1.txt
ls *_${state}.bed >> 1.txt
cat 1.txt | cut -d "_" -f 1 | perl -p -e 's/\n/ /g'| sed '$ s/.$/\n/' > AARegulatory_module/header.txt
#paste all tissues together
done

vim 01_state_merge.sh
#!/bin/bash
#Record the presence of all chromatin reads as 1 and their absence as 0
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
ls *${state}_Gs.bed  | while read id;
do
  echo $id
  cat ${id} | cut -f 4 | paste -s >> 2.txt
  done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 *_${state}_Gs.bed > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' >6.txt
cat header.txt 6.txt |sed 's/ /\t/g' > all_${state}_Gs.csv
done
#sbatch -c 5  -p low --mem 10G 01_state_merge.sh

#normalized the one count and count the number for each region
vim 02_normalization.sh
#!/bin/bash
#Normalize one count and count the number for each region
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
echo all_${state}_Gs.csv
cat all_${state}_Gs.csv | cut -f 4- |  sed 's/40/1/g' | sed 's/41/1/g'| sed 's/42/1/g'| sed 's/43/1/g'|  sed 's/30/1/g' | sed 's/31/1/g'| sed 's/32/1/g'| sed 's/33/1/g'|sed 's/34/1/g'| sed 's/35/1/g'| sed 's/36/1/g'| sed 's/37/1/g'| sed 's/38/1/g'| sed 's/39/1/g'| sed 's/20/1/g' | sed 's/21/1/g'| sed 's/22/1/g'| sed 's/23/1/g'|sed 's/24/1/g'| sed 's/25/1/g'| sed 's/26/1/g'| sed 's/27/1/g'| sed 's/28/1/g'| sed 's/29/1/g'| sed 's/10/1/g'| sed 's/11/1/g'| sed 's/12/1/g'| sed 's/13/1/g'| sed 's/14/1/g'| sed 's/15/1/g' | sed 's/16/1/g'| sed 's/17/1/g'|  sed 's/18/1/g'| sed 's/19/1/g'| sed 's/2/1/g'| sed 's/3/1/g'| sed 's/4/1/g'| sed 's/5/1/g' | sed 's/6/1/g'| sed 's/7/1/g'|  sed 's/8/1/g'| sed 's/9/1/g' > 3.txt
cat 3.txt | awk '{for(i=1;i<=NF;i++){a[NR]+=$i}print $0,a[NR]}' > 4.txt
cut  -f 1-3 all_${state}_Gs.csv > 5.txt
paste  5.txt 4.txt |sed 's/ /\t/g' > all_${state}_Gs_one_count.csv
done
#sbatch -c 5  -p low --mem 10G 02_normalization.sh

vim 03_TSR_tissue.sh
#!/bin/bash
#Because there are more than 30 tissues, tissues from the same system are treated as one tissue and are not excluded when identifying tissue-specific regulatory elements.
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10
do
echo $state
cd /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module
mkdir AA_TSR_${state}
cat all_${state}_Gs_one_count.csv | awk '{if (($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30+$31+$32+$33+$34+$35+$36+$37+$38+$39+$40+$41+$42+$43+$44+$45+$46)==43) print $0}' > AA_TSR_${state}/TSR_All_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43)==12 \
&& ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_Nervous_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($39+$38+$30+$4)==4 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Stomach_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($15+$21+$20+$8+$12+$37)==6 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Digestive_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($11+$13+$14)==3 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Uterus_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($6+$25+$45+$46+$42)==5 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Immune_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($17+$29)==2 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Muscular_System_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($40+$41)==2 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5)==0) print $0}' \
> AA_TSR_${state}/TSR_Skin_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($32+$33+$26)==3 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Female-reproductive_common_${state}.txt

cat all_${state}_Gs_one_count.csv | awk '{if (($16+$44)==2 \
&& ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$26+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' \
> AA_TSR_${state}/TSR_Male-reproductive_common_${state}.txt

cat all_${state}_Gs_one_count.csv | awk '{if ($10==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_cerebral-cortex_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($28==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_midbrain_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($9==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_cerebellum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($7==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_brainstem_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($18==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_hippocampus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($19==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_hypothalamus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($27==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_medulla-oblongata_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($31==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_optic-chiasm_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($34==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_pineal_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($35==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_pituitary_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($36==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_pons_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($43==1 && ($4+$5+$6+$8+$11+$12+$13+$14+$15+$16+$17+$20+$21+$22+$23+$24+$25+$26+$29+$30+$32+$33+$37+$38+$39+$40+$41+$42+$44+$45+$46)==0) print $0}' > AA_TSR_${state}/TSR_splenium_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($39==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_rumen_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($38==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_reticulum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($30==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_omasum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($4==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_abomasum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($15==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_duodenum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($21==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_jejunum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($20==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_ileum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($8==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_cecum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($12==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_colon_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($37==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_rectum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($11==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$13+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_cervix_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($13==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_cornua-uteri_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($14==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$13+$37+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_corpus-uteri_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($32==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_ovary_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($33==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_oviduct_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($16==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$26+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_epididymis_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($44==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$26+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_testis_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($26==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$16+$44+$6+$25+$45+$46+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_mammary-gland_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($6==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_bone-marrow_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($25==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_lymph-node_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($45==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_thymus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($46==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_thyroid_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($42==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$22+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_spleen_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($23==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$22+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_liver_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($22==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$24+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_kidney_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($24==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$17+$29+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_lung_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($17==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_heart_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($29==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$5+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_muscle_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($5==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$40+$41)==0) print $0}' > AA_TSR_${state}/TSR_adipose_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($40==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5)==0) print $0}' > AA_TSR_${state}/TSR_skin_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($41==1 && ($10+$28+$9+$7+$18+$19+$27+$31+$34+$35+$36+$43+$39+$38+$30+$4+$15+$21+$20+$8+$12+$37+$11+$13+$14+$32+$33+$16+$44+$26+$6+$25+$45+$46+$42+$23+$22+$24+$17+$29+$5)==0) print $0}' > AA_TSR_${state}/TSR_soft-horn_${state}.txt
done
#sbatch -c 100  -p low --mem 200G 03_TSR_tissue_10.sh
