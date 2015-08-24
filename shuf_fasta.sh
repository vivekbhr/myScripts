
echo 'script to shuffle and print a given number of fasta entries from a file. 
eg: shuf_fasta.sh <inputfasta> <num_toSample> <seed_seq> > out.fasta
Note: you cant fix a seed '



inputfa=${1}

cat ${inputfa} |\
awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' |\
shuf |\
head -n ${2} |\
awk 'BEGIN{FS="\t"} {printf("%s\n%s\n",$1,$2)}'
