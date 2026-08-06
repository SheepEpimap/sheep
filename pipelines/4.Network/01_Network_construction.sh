vim 01_natwork.sh
#!/bin/bash
########1.First identify enhancers
#######2.Then identify enhancer target genes#######
#######3.Then identify enhancer motifs using MEME FIMO#######
#######3.Then add the tissue-specific results#######
#######3.Then link promoters within 2 kb upstream of target genes#######
#######4.Then link TFs#######
#######/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/motif/fimo_vertebrates_bg0.001/natwork#####
cat sample1.txt | while read -r Tissue id Motif
do
    echo "$id"
    echo "$Motif"
    echo "$Tissue"
####Identify enhancers potentially bound by motifs from the FIMO results#####
    fimo="../TSR_${Tissue}_E5_enhancer/fimo.tsv"
    out="${Tissue}_${Motif}_0.bed"

    awk -F'\t' -v OFS='\t' -v id="$id" -v pthr="2.6e-05" '
        BEGIN{
            print "#chrom","start","end","motif_alt_id","score","p-value","q-value","strand"
        }
        FNR>1 && index($1, id)>0 && ($8+0) < pthr {
            print $3, $4, $5, $2, $7, $8, $9, $6
        }
    ' "$fimo" > "$out"
bedtools intersect -wo -a ${Tissue}_${Motif}_0.bed -b /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/TSR_${Tissue}_E5_id.bed  > ${Tissue}_${Motif}_enhancer.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $4,$12,$5,"'"${Tissue}"'"}' ${Tissue}_${Motif}_enhancer.txt > ${Tissue}_${Motif}_TF_Enhn_node.txt
awk -F'\t' 'BEGIN{OFS="\t"} FNR>1{print $12,"Enhancer","'"${Tissue}"'"; print $4,"Motif","'"${Tissue}"'"}' ${Tissue}_${Motif}_enhancer.txt | awk '!seen[$0]++' > ${Tissue}_${Motif}_TF_Enhn_edge.txt
#####Determine motif-enhancer-gene relationships from enhancer-gene pairs########
awk 'NR==FNR { keep[$2]=1; next } ($1 in keep)' \
    ${Tissue}_${Motif}_TF_Enhn_node.txt /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/AA_Target_gene/TSR_${Tissue}_E5_target_gene.txt > TSR_${Tissue}_E5_target_gene.intersect.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$7,$5,"'"${Tissue}"'"}' TSR_${Tissue}_E5_target_gene.intersect.txt > ${Tissue}_${Motif}_Enhn_Gene_node.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,"Enhancer","'"${Tissue}"'"; print $7,"Gene","'"${Tissue}"'"}' TSR_${Tissue}_E5_target_gene.intersect.txt | awk '!seen[$0]++' > ${Tissue}_${Motif}_Enhn_Gene_edge.txt
######Identify gene promoters from the 2-kb upstream regions of genes
awk 'NR==FNR { keep[$7]=1; next } ($4 in keep)' \
    TSR_${Tissue}_E5_target_gene.intersect.txt /vol2/mengzhu/genome/part_change_esemb100/TSS_esemble100_colin.bed_up2k.bed > TSR_${Tissue}_E5_gene_tss2k.txt
bedtools intersect -wo -a TSR_${Tissue}_E5_gene_tss2k.txt -b /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E1/TSR_${Tissue}_E1_id.bed  > TSR_${Tissue}_E1_gene.txt
#awk -F'\t' 'BEGIN{OFS="\t"} {print $11,$4,"1","'"${Tissue}"'"}' TSR_${Tissue}_E1_gene.txt > ${Tissue}_${Motif}_Promo_Gene_node.txt
awk -F'\t' -v OFS='\t' -v tissue="$Tissue" '
  $4 !~ /^LOC/ { print $11, $4, "1", tissue }
' TSR_${Tissue}_E1_gene.txt > ${Tissue}_${Motif}_Promo_Gene_node.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,"Promoter","'"${Tissue}"'"; print $2,"Gene","'"${Tissue}"'"}' ${Tissue}_${Motif}_Promo_Gene_node.txt | awk '!seen[$0]++' > ${Tissue}_${Motif}_Promo_Gene_edge.txt
#####Filter genes in ${Tissue}_${Motif}_Enhn_Gene_node.txt by promoter presence, retaining only genes with promoters
#####Extract enhancers#####
awk -F'\t' '
  NR==FNR {
    if ($2 !~ /^LOC/) keep[$2]=1
    next
  }
  ($2 in keep)
' ${Tissue}_${Motif}_Promo_Gene_node.txt ${Tissue}_${Motif}_Enhn_Gene_node.txt \
> ${Tissue}_${Motif}_Enhn_Gene_node_2.txt
#awk 'NR==FNR { keep[$2]=1; next } ($2 in keep)' ${Tissue}_${Motif}_Promo_Gene_node.txt ${Tissue}_${Motif}_Enhn_Gene_node.txt > ${Tissue}_${Motif}_Enhn_Gene_node_2.txt
awk 'NR==FNR { keep[$1]=1; next } ($1 in keep)' ${Tissue}_${Motif}_Enhn_Gene_node_2.txt ${Tissue}_${Motif}_Enhn_Gene_edge.txt > ${Tissue}_${Motif}_Enhn_Gene_edge_2.txt
#####Filter enhancers in ${Tissue}_${Motif}_TF_Enhn_node.txt by gene links, retaining only enhancers linked to genes
#####Extract motifs#####
awk 'NR==FNR { keep[$1]=1; next } ($2 in keep)' ${Tissue}_${Motif}_Enhn_Gene_node_2.txt ${Tissue}_${Motif}_TF_Enhn_node.txt > ${Tissue}_${Motif}_TF_Enhn_node_2.txt
awk 'NR==FNR { keep[$1]=1; next } ($1 in keep)' ${Tissue}_${Motif}_TF_Enhn_node_2.txt ${Tissue}_${Motif}_TF_Enhn_edge.txt > ${Tissue}_${Motif}_TF_Enhn_edge_2.txt
#####Merge to generate node and edge files######
cat ${Tissue}_${Motif}_TF_Enhn_node_2.txt ${Tissue}_${Motif}_Enhn_Gene_node_2.txt ${Tissue}_${Motif}_Promo_Gene_node.txt > node/${Tissue}_${Motif}_node.txt
cat ${Tissue}_${Motif}_TF_Enhn_edge_2.txt ${Tissue}_${Motif}_Enhn_Gene_edge_2.txt ${Tissue}_${Motif}_Promo_Gene_edge.txt > node/${Tissue}_${Motif}_edge.txt
done
#sbatch -c 50  -p low --mem 100G 01_natwork.sh
