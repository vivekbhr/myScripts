#!/usr/bin/env Rscript
## get input folder and output names as tsv
Args <- commandArgs(trailingOnly = TRUE)
indir <- Args[1]
outdir <- Args[2]
mynames <- read.delim(Args[3], header = F, stringsAsFactors = F)

makecmd <- function(folder, out) {
  lapply(c("R1", "R2"), function(x) {
    fs <- list.files(folder, pattern = x, full.names = TRUE)
    fs <- paste0(fs, collapse = " ")
    cmd <- paste0("cat ", fs, " > ", out, "_", x, ".fastq.gz")
    return(cmd)
  }) -> li
  return(li)
  } 

## execute cmd
f <- "link_cmd.sh"
if(file.exists(f)) file.remove(f)
apply(mynames, 1, function(r){
  cmd <- makecmd(file.path(indir, r[1]), file.path(outdir, r[2]))
  cat(cmd[[1]], file = f, append = TRUE, sep = "\n")
  cat(cmd[[2]], file = f, append = TRUE, sep = "\n")
}) -> put