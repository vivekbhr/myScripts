
## The script makes links of all fastq files from a given dir in akhtar/sequencing_data

seqdata=${1}

for dir in $(ls $seqdata); 
do for file in $(ls $seqdata/$dir | grep _R1.fastq.gz); 
do base=$(basename $seqdata/$dir/$file _R1.fastq.gz);
echo "Linking sample : ${base}" 
ln -s ${seqdata}/${dir}/${base}_R1.fastq.gz . ; 
ln -s ${seqdata}/${dir}/${base}_R2.fastq.gz . ; 
done; done
