

echo 'This scripts takes input fastq dir, output (links) dir and samplesheet.csv file, and makes the links to fastq dir with proper names.'

seqdir=${1}
mydir=${2}
samplesheet=${3}

awk 'FS="," {print "ln -s '${seqdir}/'"$3"_R1.fastq.gz" " " "'${mydir}/'"$4"_R1.fastq.gz"} ' ${samplesheet} | sed 1d > links.sh
awk 'FS="," {print "ln -s '${seqdir}/'"$3"_R2.fastq.gz" " " "'${mydir}/'"$4"_R2.fastq.gz"} ' ${samplesheet} | sed 1d >> links.sh
chmod ug+rwx links.sh
./links.sh
rm links.sh


