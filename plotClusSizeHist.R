
#!/package/R-3.2.0/bin/R
## plot the histogram of cluster sizes produced by deeptools and hirarchichal clustering

library("reshape")
library("ggplot2")

args = commandArgs(trailingOnly=TRUE)

infile <- args[1] # input bedfile (eg chip sites)
outfile <- args[2] # hitogram name (png will be added)
motOut <- args[3] # motif outout folder
size <- args[4] # how many bp up and down to take
        
        

data <- read.table(infile,header=FALSE,comment.char="#",stringsAsFactors = FALSE)
data$size <- data$V3 - data$V2
# plot hist
png(paste0(outfile,".png"))
hist(data$size,col="steelblue")
dev.off()
print(paste0("Histogram plotted : ",paste0(outfile,".png")))
# write the file for motif analysis

if(!(is.null(motOut))){
        print("Splitting cluster for motif analysis")
        system(paste0("mkdir ",motOut))
        system(paste0("cp ",infile," ",motOut))
        setwd(motOut)
        mkfiles <- paste0("csplit -z --suppress-matched ",infile," /\\#Cluster_*/ {*} -f ",infile,"-split_")
        system(mkfiles,intern = TRUE)
        files <- list()
        
        tryCatch({
                for(file in list.files(pattern = paste0(infile,"-split_"))){
                data <- read.table(file,blank.lines.skip = TRUE,header = FALSE)
                data$size <- data$V3 - data$V2
                data$middle <- data$V2 + as.integer(data$size/2)
                data$begin <- data$middle - as.numeric(size)
                data$end <- data$middle + as.numeric(size)
                out <- data[,c("V1","begin","end")]
                write.table(out,paste0(file,"_forMotifs.bed"),sep="\t",quote=FALSE,row.names = F,col.names = F)
                }
        }, error = function(e) print(conditionMessage(e))
        )
        
        print(paste0("DONE..! Output in folder :",motOut))
        } else { print("DONE..! Find output ") }
