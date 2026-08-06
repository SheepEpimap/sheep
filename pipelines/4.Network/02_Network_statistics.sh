cat 02_data_statistics.sh
#!/bin/bash
ls *_TF_Enhn_node_2.txt |cut -d '_' -f 1-2 | while read sample;
do
echo ${sample}
cat ${sample}_Enhn_Gene_node.txt | cut -f 1 | awk '!seen[$1]++' > ${sample}_enhancer_region.txt
cat ${sample}_Enhn_Gene_node.txt | cut -f 2 | awk '!seen[$1]++' > ${sample}_enhancer_gene.txt
all=$(cat ${sample}_node.txt |  wc -l | awk '{print $1}')
promoter=$(cat ${sample}_Promo_Gene_node.txt | cut -f 1 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
promoterGene=$(cat ${sample}_Promo_Gene_node.txt | cut -f 2 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
enhancer=$(cat ${sample}_Enhn_Gene_node.txt | cut -f 1 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
enhancerGene=$(cat ${sample}_Enhn_Gene_node.txt | cut -f 2 | awk '!seen[$1]++' |  wc -l | awk '{print $1}')
rm ${sample}_region_count.txt
cat ${sample}_enhancer_region.txt | while read id;
do
Reads1=$(grep -w $id ${sample}_Enhn_Gene_node.txt | wc -l | awk '{print $1}')
echo $id $Reads1 >> ${sample}_region_count.txt
done
rm ${sample}_gene_count.txt
cat ${sample}_enhancer_gene.txt | while read id;
do
Reads2=$(grep -w $id ${sample}_Enhn_Gene_node.txt | wc -l | awk '{print $1}')
echo $id $Reads2 >> ${sample}_gene_count.txt
done

  read enhancer_avg enhancer_med < <(
    cut -d' ' -f2 ${sample}_region_count.txt \
    | awk '/^[0-9]+(\.[0-9]+)?$/' \
    | sort -n \
    | awk '{x[NR]=$1; s+=$1}
           END{
              if (NR==0){print "NA\tNA"; exit}
              if (NR%2){m=x[(NR+1)/2]} else {m=(x[NR/2]+x[NR/2+1])/2}
              print s/NR "\t" m
           }'
  )

  read Gene_avg Gene_med < <(
    cut -d' ' -f2 ${sample}_gene_count.txt \
    | awk '/^[0-9]+(\.[0-9]+)?$/' \
    | sort -n \
    | awk '{x[NR]=$1; s+=$1}
           END{
              if (NR==0){print "NA\tNA"; exit}
              if (NR%2){m=x[(NR+1)/2]} else {m=(x[NR/2]+x[NR/2+1])/2}
              print s/NR "\t" m
           }'
  )

  # Write the current-sample summary (adding four columns: region_mean, region_median, gene_mean, and gene_median)
  echo "$all" "$sample" "$enhancer" "$enhancerGene" "$enhancer_avg" "$enhancer_med" "$Gene_avg" "$Gene_med" "$promoter" "$promoterGene" \
    > ${sample}_summary.txt

  # Append to the combined table
  cat ${sample}_summary.txt >> all_target_gene_summary.txt
done

# Convert to tab-delimited format
sed 's/ /\t/g' all_target_gene_summary.txt > all_target_gene_summary.csv
#sbatch -c 50  -p low --mem 100G 02_data_statistics.sh

