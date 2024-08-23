indir=$1
regex=$2
version=$3

if [ -z "$3" ]; then
    echo "fetching fastqs"
    for name in $(ls $indir | grep $regex | grep _R1 | sed 's/_[A-Z0-9_]*_S[0-9]*_L[0-9]*//' \
    | sed s/_R1_001.fastq.gz// | sort | uniq)
    do 
    cat $indir/${name}_*_R1_001.fastq.gz > ${name}_R1.fastq.gz && \
    cat $indir/${name}_*_R2_001.fastq.gz > ${name}_R2.fastq.gz && \
    echo $indir/${name}_*_R1_001.fastq.gz $name "Done" >> rename.txt &
    done

elif [ $version=="v1" ]; then
    echo "fetching fastqs from older illumina runs (pre 2022)"
    for dir in $(find ${indir} -type f -name $regex | sed -r 's|/[^/]+$||' |sort |uniq);
    do
    for R in R1 R2; do cat $dir/*_${R}_*.fastq.gz > $(basename $dir)_${R}.fastq.gz;
    	echo "$(basename $dir)_${R}.fastq.gz $(ls $dir | grep ${R} | grep L001)" >> rename.sh
    done; done
    cat rename.sh | sed s/_[A-Z0-9_]*_S[0-9]*_L[0-9]*// | sed s/_001// | xargs -L 1 mv

fi