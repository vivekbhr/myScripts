#!/package/R-3.2.0/bin/Rscript

library(dplyr)
library(ggplot2)
library(reshape2)
args=commandArgs(trailingOnly = TRUE)
infolder=args[1]
setwd(infolder)
read.table("rDNA_counts_excludeMulti.txt",header=T,sep="\t") -> ribocounts 
sampledreads = 1000000
## make stats
ribocounts.n <- data.frame(Filename = ribocounts$Filename,
                           Total = ribocounts$Total.alignments/ribocounts$Total.alignments * 100,
                           perc.primary= ribocounts$Primary.alignments/ribocounts$Total.alignments * 100,
                           perc.primary.rDNA = ribocounts$Primary.rDNA.alignments/ribocounts$Total.alignments
)
ribocounts.num <- data.frame(Filename = ribocounts$Filename,
                             Total.reads = ribocounts$Total.reads/ribocounts$Total.reads * 100,
                             Total.aligned.reads = ribocounts$Mapped.reads/ribocounts$Total.reads * 100,
                             Primary.aligned.reads = ribocounts$Primary.Mapped.reads/ribocounts$Total.reads * 100,
                             rDNA.aligned.reads = ribocounts$rDNA.reads/sampledreads * 100,
                             rDNA.aligned.primary.reads = ribocounts$rDNA.reads.primary/sampledreads * 100
)
## plot out
pdf("ribocont_results.pdf")
ribocounts.n %>% melt() %>% 
  ggplot(.,aes(Filename,value,fill=variable)) + geom_bar(stat = "identity", position = "dodge") + coord_flip()
ribocounts.num %>% melt() %>% 
  ggplot(.,aes(Filename,value,fill=variable)) + geom_bar(stat = "identity", position = "dodge") + 
  scale_y_continuous(breaks = seq(0,100,10))  + coord_flip() + geom_hline(aes(yintercept=10))
dev.off()
