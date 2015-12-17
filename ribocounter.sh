#!/bin/bash

echo "counts total,mapped,primary and primary ribosomal RNA reads in a sample. Output = rDNA_counts_excludeMulti.txt"
folder=${1}
echo "add ribocounter.R to your PATH"

cd ${folder}
# make headers
echo -e "Filename" >> Filename.txt
echo -e "Total.alignments" >> total_alignments.txt
echo -e "Primary.alignments" >> primary_alignments.txt
echo -e "Primary.rDNA.alignments" >> Primary_rDNA_alignments.txt
echo -e "Total.reads" >> Total_reads.txt
echo -e "Mapped.reads" >> Mapped_reads.txt
echo -e "Primary.Mapped.reads" >> Primary-Mapped_reads.txt
echo -e "rDNA.reads" >> rDNA_reads.txt
echo -e "rDNA.reads.primary" >> rDNA_reads_primary.txt

## count and output
for file in $(ls | grep .bam$)
do
echo ${file} >> Filename.txt
samtools view -c ${file} >> total_alignments.txt
samtools view -cF260 ${file} >> primary_alignments.txt
samtools view -cF260 ${file} rDNA >> Primary_rDNA_alignments.txt
## Now count reads
samtools bam2fq ${file} | grep "@M" | wc -l >> Total_reads.txt
samtools view -bF4 ${file} | samtools bam2fq - | grep "@M" | wc -l >> Mapped_reads.txt
samtools view -bF260 ${file} | samtools bam2fq - | grep "@M" | wc -l >> Primary-Mapped_reads.txt

#samtools view ${file} | grep "NH:i:[0-9]" | samtools bam2fq - | grep "/1" | wc -l >> Multi-Mapped_reads.txt
#samtools view -bF260 ${file} | samtools bam2fq - | grep "/1" | wc -l >> Primary_Multi-Mapped_reads.txt
samtools view -bF4 ${file} rDNA | samtools bam2fq - | grep "@M" | wc -l >> rDNA_reads.txt
samtools view -bF260 ${file} rDNA | samtools bam2fq - | grep "@M" | wc -l >> rDNA_reads_primary.txt
done

# paste results and remove other files
paste Filename.txt total_alignments.txt primary_alignments.txt Primary_rDNA_alignments.txt Total_reads.txt Mapped_reads.txt Primary-Mapped_reads.txt rDNA_reads.txt rDNA_reads_primary.txt > rDNA_counts_excludeMulti.txt

rm Filename.txt total_alignments.txt primary_alignments.txt Primary_rDNA_alignments.txt Total_reads.txt rDNA_reads.txt Mapped_reads.txt Primary-Mapped_reads.txt rDNA_reads_primary.txt

## Run the Rscript

ribocounter.R $(pwd)
