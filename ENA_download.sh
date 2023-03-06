#!/bin/sh
input=$1
njobs=$2
num=$(head -1 $input | tr '\t' '\n' | grep -n "fastq_ftp" | cut -d: -f1)
name=$(head -1 $input | tr '\t' '\n' | grep -n "sample_title" | cut -d: -f1)
sed 1d $input | cut -f$num | cut -d ';' -f1 > fname1
sed 1d $input | cut -f$num | cut -d ';' -f2 > fname2
sed 1d $input | cut -f$name > fname3
sed -i 's/ /_/g' fname3
## Fetch read1
echo "Fetching Read 1 and unpaired reads with suffix: _R1.fastq.gz"
paste fname1 fname3 | awk 'OFS="\t"{print $2"_R1.fastq.gz", $1}' | xargs -P $njobs -L 1 -n 2 wget -O
## fetch read2
echo "Fetching Read 2 for paired reads with suffix: _R2.fastq.gz"
paste fname2 fname3 | awk '/_2.fastq.gz/ {print $2"_R2.fastq.gz", "\t", $1}' | xargs -P $njobs -L 1 -n 2 wget -O
