

### Scripts to perform small tasks, but can be generalized and be very useful

1) Awk script to remove Zeros before plotting into multiheatmapper
awk '{IFS="\t"; if (/^@/) print $0; else Count = 0; for (i = 10; i <= NF; i++) Count = Count+$i; if(Count != 0) print $0 }' input > outputfilename_zeroremoved.tab

# Eg. for zipped .tab file
awk '{IFS="\t"; if (/^@/) print $0; else Count = 0; for (i = 10; i <= NF; i++) Count = Count+$i; if(Count != 0) print $0 }' <(gzip -dc multideeptools_compositematrix_AllChip_allcells.tab.gz | gzip) > multideeptools_compositematrix_AllChip_allcells_zeroremoved.tab.gz

2) Open a file and add a delimiter to a column then grep it from another file 

awk '{print $1"\t"}' file1 | grep -f - file2 > outputfile

3) convert vcf to bed 

sed -e 's/chr//' file.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > out.bed

