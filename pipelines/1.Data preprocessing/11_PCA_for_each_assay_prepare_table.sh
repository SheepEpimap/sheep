#Prepare _PCA_last.txt
sed '1d' H3K4me3_PCA.tab > H3K4me3_PCA_no_header.tab
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j]
        }
        print str
    }
}' H3K4me3_PCA_no_header.tab > H3K4me3_PCA_transposed.tab
sed -i '$d' H3K4me3_PCA_transposed.tab
cut -f1 H3K4me3_PCA_transposed.tab | cut -d '_' -f 2 >2
cut -f1 H3K4me3_PCA_transposed.tab | cut -d '_' -f 3 >3
paste H3K4me3_PCA_transposed.tab 2 3 >H3K4me3_PCA_last.txt
#Finally, Group and Rep still need to be updated
