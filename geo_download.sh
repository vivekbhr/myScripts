#!/usr/bin/env bash
# geo_download.sh
#
# Example
# -------
# To download all of the read archive files for GEO12001:
#   geo_download GEO012001
#
# For 'esearch', 'efetch', 'xtract', you must install Entrez Direct:
#   http://www.ncbi.nlm.nih.gov/books/NBK179288/
#
# For 'fastq-dump', you must install SRA Tools:
#   https://github.com/ncbi/sra-tools


GEO=${1}
echo "expecting esearch installed in /hpc/hub_oudenaarden/vbhardwaj/programs"
echo "expecting parallel fastq-dump under module files"

export export PATH=/hpc/hub_oudenaarden/vbhardwaj/programs/edirect/:$PATH 
# get SRA from geo
SRA=$(esearch -db gds -query "${GEO}" | efetch -format docsum | xtract -pattern ExtRelation -element RelationType,TargetObject | cut -f2 | grep SRP)

# get SRR from SRA
esearch -db sra -query $SRA | \
  efetch -format docsum | \
  xtract -pattern DocumentSummary -element Run@acc | \
  tr '\t' '\n' > ${SRA}.txt

# use Experiment@name Run@acc : for experiment name as well
# module load parallel-fastq-dump
# get fastq from SRRs
while read -r LINE; do
    parallel-fastq-dump -t 10 --split-3 --gzip -s "$LINE" &
done < ${SRA}.txt
