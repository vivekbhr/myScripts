## Vivek 4/4/15
## This script runs samtools flagstat on the mapped BAM file and makes a plot of the output
## Input = bam file, usage "Rscript mapstat.R <input.bam> <samplename>"
library('reshape')
library('ggplot2')
args = commandArgs(trailingOnly = TRUE)
wd = system('pwd',intern = T)
setwd(wd)

print(paste('Running flagstat for',args[2],sep=" "))
system(command = paste0('/package/samtools-1.1/samtools flagstat ',args[1],' > flagstat_',args[1],'.out'),wait = TRUE)
print('Plotting results')
map.tot = as.data.frame(read.table(pipe("cut -f1 -d ' ' 'flagstat.out'")))
#map.filt = read.table(pipe("cut -f1 -d ' ' '/data/akhtar/bhardwaj/2015_Ibu_ChIP-exo/mySeqTest_MSL1/HISAT_mapping/flagstat_msl1-filtered.txt'"))
#map.tot = cbind(map.unfil,map.filt)
colnames(map.tot) = args[2]
map.tot$status = c('Total','Secondary','Suppli','Duplicates','Mapped','Paired','Read1','Read2',
                   'Properly Paired','Both Mapped','Singletons','Mate on diff Chr','Mate on diff Chr Q>5')
map.tot = map.tot[c(1,5:12),]                       
map.tot= melt(map.tot)
map.tot$status = factor(map.tot$status, levels = c('Total','Secondary','Mapped','Paired','Read1','Read2',
                                                   'Properly Paired','Both Mapped','Singletons','Mate on diff Chr'))
png(paste(args[2],"png",sep="."),width = 800,height = 500)
ggplot(map.tot, aes(status,value, group= variable)) + geom_bar(stat='identity',position='dodge') +
  theme_gray(base_size = 15) + theme(axis.text.x=element_text(angle = 60,vjust = 0.6)) +
  labs(x='Status',y='No. of Reads',title= paste('Mapping Stats: ',args[2],sep=""))

dev.off()

print("Done..!")