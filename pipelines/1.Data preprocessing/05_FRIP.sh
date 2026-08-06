#!/bin/bash
rm narrow_frip.txt
ls *_peaks.narrowPeak | cut -d '_' -f 1-3 | while read id; do
    echo "$id"
    bed="/vol2/mengzhu/snakemake_sheep/clean/bam1/${id}.bed"
    peak="/vol2/mengzhu/snakemake_sheep/peak3/${id}_peaks.narrowPeak"
    Reads=$(bedtools intersect -a "$bed" -b "${peak}" | wc -l | awk '{print $1}')
    totalReads=$(wc -l "$bed" | awk '{print $1}')
    totalpeaks=$(wc -l "${peak}" | awk '{print $1}')
    FRIP=$(bc <<< "scale=2; 100 * $Reads / $totalReads")
    echo -e "$id\t$totalReads\t$totalpeaks\t$FRIP" >> "narrow_frip.txt"
done

vim broad_frip.sh
#!/bin/bash
ls *_peaks.broadPeak | cut -d '_' -f 1-3 | while read id; do
    echo "$id"
    bed="/vol2/mengzhu/snakemake_sheep/clean/bam1/${id}.bed"
    peak="/vol2/mengzhu/snakemake_sheep/peak3/${id}_peaks.broadPeak"
    Reads=$(bedtools intersect -a "$bed" -b "${peak}" | wc -l | awk '{print $1}')
    totalReads=$(wc -l "$bed" | awk '{print $1}')
    totalpeaks=$(wc -l "${peak}" | awk '{print $1}')
    FRIP=$(bc <<< "scale=2; 100 * $Reads / $totalReads")
    echo -e "$id\t$totalReads\t$totalpeaks\t$FRIP" >> "broad_frip.txt"
done
sbatch -c 10  -p low --mem 10G broad_frip.sh

#Calculate genome coverage
#!/bin/bash
ls ATAC*_peaks.narrowPeak | while read id;
do
  echo ${id}
  awk -F '\t' '!/^#/ {sum += $4} END {print "'${id}'", sum / 2628104905}' ${id} >> coverage.txt
done

#!/bin/bash

# Total genome length
genome=2628104905
ls ATAC*_peaks.narrowPeak | while read id; do
  sum=$(awk '{ sum += ($3 - $2) } END { print sum }' "${id}")
  # Calculate coverage
  coverage=$(awk -v s="${sum}" -v g="${genome}" 'END { printf("%.6f", s/g) }' <<< "")
  echo "${id}  ${sum}  ${coverage}" >> coverage.txt
done
