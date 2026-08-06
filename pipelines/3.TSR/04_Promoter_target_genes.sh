#Try merging chromatin states E1, E2, E3, and E9 before identifying target genes (first inspect the target-gene results for each state)
cat E1_Gs.bed E2_Gs.bed E3_Gs.bed E3_Gs.bed |sort -k1,1 -k2,2n > promoter.bed
bedtools merge -i promoter.bed > promoter_merge.bed
cat E4_Gs.bed E5_Gs.bed E6_Gs.bed E7_Gs.bed |sort -k1,1 -k2,2n > enhancer.bed
bedtools merge -i enhancer.bed > enhancer_merge.bed
ls *_E1_Gs.bed | cut -d "_" -f 1 | while read id ; do
echo ${id}
cat AA_TSR_E1/TSR_${id}_E1_id.bed AA_TSR_E2/TSR_${id}_E2_id.bed AA_TSR_E3/TSR_${id}_E3_id.bed AA_TSR_E9/TSR_${id}_E9_id.bed |sort -k1,1 -k2,2n > promoter.bed
bedtools merge -i promoter.bed > promoter_${id}.bed
done

####Use bedtools intersect to identify the gene nearest each promoter as its target gene#####
mkdir -p AA_Target_gene/Sheep_TSSup2k
ls *_id.bed | while read id; do
    echo "$id"
    bedtools intersect -a "$id" \
        -b /vol2/mengzhu/genome/part_change_esemb100/TSS_esemble100_colin.bed_up2k.bed \
        -wa -wb > AA_Target_gene/Sheep_TSSup2k/"$(basename "$id" "_id.bed")"_gene_up2k.txt
done
cd AA_Target_gene/Sheep_TSSup2k
ls *_gene_up2k.txt | while read id;
do
echo $id
cat $id | cut -f 8 | sort > $(basename $id "_gene_up2k.txt")_gene.txt 
done

cp AA_TSR_E1/AA_Target_gene/Sheep_TSSup2k/*_gene.txt AA_TSR_E2/AA_Target_gene/Sheep_TSSup2k/*_gene.txt AA_TSR_E3/AA_Target_gene/Sheep_TSSup2k/*_gene.txt AA_TSR_E4/AA_Target_gene/Sheep_TSSup2k/*_gene.txt /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/promoter_TSSup2k

ls *_E1_gene.txt | while read id ; do
echo ${id}
cat $(basename $id "_E1_gene.txt")_E1_gene.txt $(basename $id "_E1_gene.txt")_E2_gene.txt $(basename $id "_E1_gene.txt")_E3_gene.txt $(basename $id "_E1_gene.txt")_E4_gene.txt \
  | sort -u > $(basename $id "_E1_gene.txt")_promoter.bed
done



