

## Load a bunch of libraries that I mostly use (NOT WORKING STILL)

## default pkg options
usual <- c("plyr","dplyr","reshape2","ggplot2")
genomic <- c("GenomicRanges","DESeq2","DEXSeq","devtools")
plotlibs <- c("UpSetR","RColorBrewer","pheatmaps","gridExtra")
efficient <- c("parallel","doParallel","readr","data.table")


loadMyLibs <- function(groupName = "usual",silent=TRUE){
  
  ## Load packages
  pkgs <- groupName
  if(silent){
    for(pkg in pkgs) suppressPackageStartupMessages(require(pkg))
  } else {
    for(pkg in pkgs) require(pkg)
  }
  
}


