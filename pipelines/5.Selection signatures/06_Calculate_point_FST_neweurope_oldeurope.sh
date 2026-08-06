####Calculate neweurope_oldeurope_point####
####Calculate newasia_oldasia_point####
vim 01.fst_point_1.sh
#!/bin/bash
# run fst between two populations
mkdir -p neweurope_newasia_point neweurope_oldasia_point newasia_oldeurope_point oldasia_oldeurope_point
for chr in chr{1..26} 
do
        jsub -q normal -n 1 -R "span[hosts=1]" -J neweurope_oldeurope_${chr}_fst.point -e neweurope_oldeurope_point/${chr}.point.%J.log -o neweurope_oldeurope_point/${chr}.point.%J.log "bash fst_point.sh $chr neweurope oldeurope"
         jsub -q normal -n 1 -R "span[hosts=1]" -J neweurope_newasia_${chr}_fst.point -e neweurope_newasia_point/${chr}.point.%J.log -o neweurope_newasia_point/${chr}.point.%J.log "bash fst_point.sh $chr neweurope newasia"
         jsub -q normal -n 1 -R "span[hosts=1]" -J neweurope_oldasiae_${chr}_fst.point -e neweurope_oldasia_point/${chr}.point.%J.log -o neweurope_oldasia_point/${chr}.point.%J.log "bash fst_point.sh $chr neweurope oldasia"
        jsub -q normal -n 1 -R "span[hosts=1]" -J newasia_oldasia_${chr}_fst.point -e newasia_oldasia_point/${chr}.point.%J.log -o newasia_oldasia_point/${chr}.point.%J.log "bash fst_point.sh $chr newasia oldasia"
        jsub -q normal -n 1 -R "span[hosts=1]" -J newasia_oldeurope_${chr}_fst.point -e newasia_oldeurope_point/${chr}.point.%J.log -o newasia_oldeurope_point/${chr}.point.%J.log "bash fst_point.sh $chr newasia oldeurope"
        jsub -q normal -n 1 -R "span[hosts=1]" -J oldasia_oldeurope_${chr}_fst.point -e oldasia_oldeurope_point/${chr}.point.%J.log -o oldasia_oldeurope_point/${chr}.point.%J.log "bash fst_point.sh $chr oldasia oldeurope"
        
done

# combine and filter by top 1%
bash fst.combine_point.sh neweurope oldeurope
bash fst.combine_point.sh newasia oldasia

cat fst.combine_point.sh
#!/bin/bash
pop1=$1
pop2=$2
win=$3
step=$4

for chr in chr{1..26} ; do sed '1d' ${pop1}_${pop2}/${chr}.${win}_${step}.windowed.weir.fst | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t%.6f\t%.6f\n",$5,$6}' ; done | sed '1iCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST' > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst
awk '$4>=40' ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter

# filter and extract top 1% region
lines=`cat ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter | wc -l`
top_lines=$(($lines/100+1))
csvtk sort -t -k 5:Nr ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter | head -$top_lines > ${pop1}_${pop2}/chrAuto.${win}_${step}.windowed.weir.fst.filter.t0.01
