####Assign IDs to the merged chromatin states#####
####/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs####
####Enhancer####
awk 'BEGIN{OFS="\t"}{
  region=$1":"$2"-"$3;
  id=sprintf("OvisE5%06d", ++i);
  print $1,$2,$3,region,id
}' E5_Gs.bed > E5_Gs_ID.bed

awk 'BEGIN{FS=OFS="\t"}
NR==FNR { id[$4]=$5; next }                    # Read E5_Gs_ID.bed first: key=$4, value=$5
FNR==1 {                                       # Process the header: insert the new column name after Enhancer
  printf "%s\t%s\t%s\tEnhancer_ID", $1,$2,$3
  for(i=4;i<=NF;i++) printf "\t%s", $i
  printf "\n"
  next
}
{
  eid = (($3 in id) ? id[$3] : "NA")           # Enhancer is in column 3
  printf "%s\t%s\t%s\t%s", $1,$2,$3,eid
  for(i=4;i<=NF;i++) printf "\t%s", $i
  printf "\n"
}' /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs/E5_Gs_ID.bed all_tables.with_header.tsv > all_tables.with_enhancerID.tsv

####Promoter####
awk 'BEGIN{OFS="\t"}{
  region=$1":"$2"-"$3;
  id=sprintf("OvisP1%06d", ++i);
  print $1,$2,$3,region,id
}' E1_Gs.bed > E1_Gs_ID.bed

awk 'BEGIN{FS=OFS="\t"}
NR==FNR { p[$4]=$5; next }   # Read E1_Gs_ID.bed first: key=$4, value=$5

FNR==1 {                    # Header: insert Promoter_ID after Promoter
  for(i=1;i<=5;i++) printf (i==1? "%s":"\t%s"), $i
  printf "\tPromoter_ID"
  for(i=6;i<=NF;i++) printf "\t%s", $i
  printf "\n"
  next
}

{
  pid = (($5 in p) ? p[$5] : "NA")   # Promoter is in column 5
  for(i=1;i<=5;i++) printf (i==1? "%s":"\t%s"), $i
  printf "\t%s", pid
  for(i=6;i<=NF;i++) printf "\t%s", $i
  printf "\n"
}' /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AAGs/E1_Gs_ID.bed all_tables.with_enhancerID.tsv > all_tables.with_enhancer_promoterID.tsv
