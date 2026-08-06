####Calculate eur_cea_point####
vim 01.fst_point.sh
#!/bin/bash
# run fst between two populations
mkdir -p eur_cea_point
ln -s eur_cea cea_eur
for chr in chr{1..26}
do

        jsub -q normal -n 1 -R "span[hosts=1]" -J eur_cea_${chr}_fst_point -e eur_cea_point/${chr}_point.%J.log -o eur_cea_point/${chr}_point.%J.log "bash fst_point.sh $chr eur cea"
done

# combine and filter by top 1%
bash fst.combine.sh eur cea 50000 10000
bash fst.combine.sh eur cea 10000 10000
bash fst.combine.sh eur mou 50000 10000
bash fst.combine.sh cea mou 50000 10000
