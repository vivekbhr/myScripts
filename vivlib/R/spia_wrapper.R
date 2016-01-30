

#' A wrapper over SPIA pathway analysis for mouse and human
#'
#' @param df a data frame of DESeq2 output
#' @param organism select from "mmu" (mouse) or "hsa" (human)
#' @param padjCutoff FDR cutoff for the genes to include in the list
#' @param outFile output file to write affected pathways
#' @param outPlot output SPIA two-way evidence plot + output summary (pdf)
#'
#' @return a txt file and a pdf file
#' @export
#'
#' @examples
#'
#'
spia_wrapper <- function(df, organism = "mmu", padjCutoff, outFile, outPlot){
        orgmap <- data.frame(org = c("mmu","hsa"), db = c("org.Mm.eg.db","org.Hs.eg.db"),stringsAsFactors = FALSE)
        assertthat::assert_that(organism %in% orgmap$org)
        db <- orgmap[grep(organism,orgmap$org),2]
        library(db)
        ensToEntrez <- AnnotationDbi::select(db,as.character(df$Row.names),"ENTREZID", keytype = "ENSEMBL" )
        df <- merge(df,ensToEntrez,by = 1)
        df.dg <- dplyr::filter(df, padj < padjCutoff,!(is.na(ENTREZID)),!(duplicated(ENTREZID)))
        df.map <- df.dg$log2FoldChange
        names(df.map) <- df.dg$ENTREZID
        allgenes <- na.omit(as.character(df$ENTREZID)) # can take from any df. it's the universe
        spia.degenes <- SPIA::spia(df.map,allgenes,organism = organism, nB = 2000)
        spia.degenes$Name <- substr(spia.degenes$Name,1,20)
        spia.degenes = spia.degenes[order(spia.degenes$pGFWER),]
        top <- head(spia.degenes[c(1:5,7,9,11)])
        colnames(top) <- c("TOP PATHWAYS : NAME",colnames(top[2:7]))
        pdf(outPlot)
        SPIA::plotP(spia.degenes)
        gplots::textplot(top)
        dev.off()
        write.table(spia.degenes,outFile,sep="\t",quote=FALSE,row.names = FALSE)
}
