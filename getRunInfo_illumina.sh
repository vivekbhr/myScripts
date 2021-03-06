

echo "given input folder with illumina fastq and output file name. This script spits out the run info in a table"
infolder=${1}
outfile=${2}

for file in $(ls ${infolder}| grep R1 | grep "fastq");
do basename ${infolder}/${file} _R1.fastq.gz >> ${infolder}/filenames.txt
zcat ${file} | head -1 | tr ":" "\t" | cut -f1-5 | sed 's/@//g' >> ${infolder}/info.txt
done
echo -e "sample"'\t'"machine_id"'\t'"run_number"'\t'"flowcell_ID"'\t'"lane"'\t'"tile" > ${outfile}
paste ${infolder}/filenames.txt ${infolder}/info.txt >> ${outfile}
rm ${infolder}/filenames.txt ${infolder}/info.txt
