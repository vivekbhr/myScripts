#!/usr/bin/env bash
# geo_download.sh <GSE_ID>
#
# Example
# -------
# To download all of the read archive files for GEO12001:
#   geo_download GEO012001


GEO=${1}
echo "expecting esearch and parallel-fastq-dump installed"
echo "if not, install via conda install -c bioconda parallel-fastq-dump entrez-direct"

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
     echo $LINE
    parallel-fastq-dump -t 10 --split-3 --gzip -s "$LINE";
done < ${SRA}.txt
