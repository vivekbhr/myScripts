
args <- commandArgs(trailingOnly = TRUE)

args[1] -> file
as.character(args[2]) -> seperator
args[3] -> column
args[4] -> name

#data <- read.table(file,sep=seperator,header=T)
data <- read.table(file,header=T)
png(paste0(name,"_",column,"histogram.png"))
hist(data[,column],col="steelblue",main=paste0("histogram: ",name),xlab="value",ylab="frequency")
dev.off()