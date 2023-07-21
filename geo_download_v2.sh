#!/usr/bin/bash
# to create the appropriate env: conda create -n sra_tools -c bioconda -c conda-forge python>=3.5 sra-tools==3.0.0 entrez-direct
TMPDIR=/hpc/hub_oudenaarden/vbhardwaj/tempdir
SLURM=/hpc/hub_oudenaarden/vbhardwaj/programs/myScripts/SlurmEasy
THREADS=20
# get SRA from geo
GEO=$1
SRA=$(esearch -db gds -query "${GEO}" | efetch -format docsum | xtract -pattern ExtRelation -element RelationType,TargetObject | cut -f2 | grep "SRP\|SRX")
echo 'Fetched SRP ID: ' ${SRA}

# getsampleInfo
echo 'Dumping sample info to samples.txt..'
touch samples.txt
for f in $SRA;
do esearch -db sra -query $f | \
    efetch -format docsum | \
    xtract -pattern DocumentSummary -element Title Run@acc | \
    tr ' ' '_' | tr ':_' '_' | tr '; ' '_' | sed 's/__/_/g' | tr -d '(|)' >> samples.txt
done

# make one directory per sample and dump the SRA IDs in it
echo 'Creating directories and dumping fastq..'

while read LINE ;
do dir=$(echo $LINE | awk '{print $1}');
samples=$(echo $LINE | awk 'OFS=" " {for (i=2; i<=NF; i++) print $i}');
  echo $dir #&& mkdir -p ${dir} && echo $samples | tr '\t' '\n' > ${dir}/samples.txt
## prep fastq dump script
  echo -e '#!/usr/bin/env bash \n' > fq_cmd.sh && \
  echo -e "prefetch --force all --max-size 10t -O ${dir} ${samples} && " | tr '\n' ' ' >> fq_cmd.sh && \
  echo -e "cd ${dir} && vdb-validate ${samples} &&
  fasterq-dump -m 10000MB -b 100MB -c 1000MB --temp $TMPDIR -e $THREADS -p --split-files --include-technical -O . ${samples}" | tr '\n' ' ' >> fq_cmd.sh && \
  echo -e '\n' >> fq_cmd.sh && \
  chmod 777 fq_cmd.sh && \
  $SLURM -t $THREADS -m 2G fq_cmd.sh
done < samples.txt
# echo -e "export HOME=$TMPDIR \n" >> fq_cmd.sh && \

# for each SRA ID, do fastq dump
#while read -r LINE;
#  do name=$(echo $LINE | cut -d ' ' -f2) && \
#  echo -e '#!/usr/bin/env bash' > fq_cmd.sh && \
#  echo -e "export HOME=$TMPDIR" >> fq_cmd.sh && \
#  echo -e "prefetch --max-size 10t -o $name $(echo $name | cut -d ' ' -f2 | sed 's/.sra//g')" >> fq_cmd.sh && \
#  echo -e "vdb-validate $name" >> fq_cmd.sh && \
#  echo -e "fasterq-dump -m 10000MB -b 100MB -c 1000MB --temp $TMPDIR -e 20 -p --split-files --include-technical -o "$LINE"" >> fq_cmd.sh && \
#  chmod 777 fq_cmd.sh && \
#  $SLURM -t 20 --no_log fq_cmd.sh;
#done < $1
