
## convert gtf to bed 

/package/UCSCtools/gtfToGenePred -genePredExt -geneNameAsName2 ${1} $PWD/gene.tmp
awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' $PWD/gene.tmp > ${2}
rm $PWD/gene.tmp
