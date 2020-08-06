indir=$1
regex=$2

for dir in $(find ${indir} -type f -name $regex | sed -r 's|/[^/]+$||' |sort |uniq);
do
for R in R1 R2; do cat $dir/*_${R}_*.fastq.gz > $(basename $dir)_${R}.fastq.gz;
	echo "$(basename $dir)_${R}.fastq.gz $(ls $dir | grep ${R} | grep L001)" >> rename.sh
done; done

cat rename.sh | sed s/_[A-Z0-9_]*_S[0-9]*_L[0-9]*// | sed s/_001// | xargs -L 1 mv
