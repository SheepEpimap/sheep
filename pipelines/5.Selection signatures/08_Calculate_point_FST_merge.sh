vim 01.merge.sh
#!/bin/bash
{ 
    head -n 1 chr1.point.weir.fst
    awk 'FNR > 1' chr*.point.weir.fst | sort -k1,1n -k2,2n
} > chrAuto_${1}.windowed.weir.fst
#jsub -q normal -n 8 -o output.%J -e error.%J bash 01.merge.sh newasia_oldeurope_point
