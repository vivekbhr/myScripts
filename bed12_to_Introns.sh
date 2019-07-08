
## Convert bed12 of gene annotation to intron bed12
awk '{OFS="\t";
	split($11,a,","); 
	split($12,b,","); 
	A=""; B=""; 
	for(i=1;i<length(a)-1;i++) {
		A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";
		} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;
	}' ${1}
