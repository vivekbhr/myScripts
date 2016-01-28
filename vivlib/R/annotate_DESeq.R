
#' annotate DESeq Output file
#'
#' @param DESeqOutput tab seperated DESeq2 output file
#' @param Output Annotated output file
#' @param remote Whether use biomart to annotate file
#' @param genome When remote = TRUE, which genome to use? (available = "hg38","mm10","dm6")
#' @param map_file If remote = FALSE, provide a map file (with ENS id and Gene id in column 1 and 2) respectively
#'
#' @return annotated output file
#' @export
#'
#' @examples
#' annotate_DESeqOutput(DESeqOutput = "test.out", Output = "test_annotated.out", remote = TRUE, genome = "hg38")
#'

annotate_DESeqOutput <- function(DESeqOutput, Output, remote = TRUE, genome = "hg38", map_file ){

        seqout <- read.delim(DESeqOutput,header = TRUE)
        if(remote == TRUE) {
                genomes = data.frame(id = c("hg38","mm10","dm6"),
                                     path = c("hsapiens_gene_ensembl","mmusculus_gene_ensembl",
                                              "dmelanogaster_gene_ensembl"),
                                     stringsAsFactors = FALSE)

                ifelse(!(genome %in% genomes$id), exit("Genome not in list!"),
                       message("fetching annotation") )
                tofetch <- dplyr::filter(genomes,id == genome)$path
                ensembl = biomaRt::useMart(tofetch,mart=ensembl)
                ext.data <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","description"),
                                mart = ensembl)
                outfile <- merge(ext.data,seqout,by.x = 1,by.y = 0)
                write.table(outfile,Output,sep="\t",quote=FALSE)

        } else {
                assertthat::assert_that(assertthat::is.readable(map_file))
                bed <- read.delim(map_file,header = TRUE)
                outfile <- merge(bed,seqout,by.x = 1,by.y = 0)
                write.table(outfile,Output,sep="\t",quote=FALSE)
        }
}
