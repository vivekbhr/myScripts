#!/bin/bash
tsv=$1
dist=$2
genome=$3
sed 1d $tsv | awk -F '\t' -v X=$dist '{mid=(int($2)+int($3))/2;printf("%s\t%d\t%d\t%s\n", chr$1, (mid-X<0?0:mid-X), mid+X, $13);}' | awk '{print >$4}'
mkdir clusters_meme
for c in cluster_*;
do mv $c clusters_meme/$c.bed
bedtools getfasta -fi $genome -fo clusters_meme/$c.fa -bed clusters_meme/$c.bed -name
meme -oc clusters_meme/$c -dna -seed 2020 -maxsize 100000000 -minsites 10 -minw 4 -maxw 20 clusters_meme/$c.fa -nmotifs 10 -revcomp
done
