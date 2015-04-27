library(cummeRbund)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

#############################################################
#															#
#			Command line parameters setting					#
#															#
#############################################################


path_to_cuffdiff_output_folder <- args[1]				### path the the cuffdiff output folder
output_folder <- args[2]								### path where all the output plots will be saved
gene_list_file <- args[3]								### path to the text file with the gene list of interest

#### Usage: Rscript cufflinks.R /path/to/cuffdiff/output /path/to/genelist/file/ /path/to/output/directory/

dir.create(output_folder, showWarnings = FALSE)			### create the output directory if it does not exist

print(output_folder)

#############################################################
#															#
#			READ and parse cuffdiff output					#
#															#
#############################################################

l <- strsplit(path_to_cuffdiff_output_folder ,'/')
cuffdiff_output_name <- l[[1]][length(l[[1]])]								## to get cuffdiff output name from the path
cuff<-readCufflinks(dir=path_to_cuffdiff_output_folder, rebuild=T)			### reading the cuffdiff output into R


print(paste(output_folder,cuffdiff_output_name,sep=""))
dir.create(paste(output_folder,cuffdiff_output_name,sep=""), showWarnings = FALSE)


gene_ids <- c()
AddItemNaive <- function(item)
{
    .GlobalEnv$gene_ids[[length(.GlobalEnv$gene_ids)+1]] <- item
}




############################################################
#															#
#		Global statistics and Quality Control				#
#															#
#############################################################


print(paste("\n ------- ",cuffdiff_output_name," ------- \n",sep=""))


dir.create(paste(output_folder,cuffdiff_output_name,"/","stats_and_QC",sep=""), showWarnings = FALSE)

#### Dispersion plot

pdf(file=paste(output_folder,cuffdiff_output_name,"/","stats_and_QC/","dispersionPlot_",cuffdiff_output_name,".pdf",sep=""))
 
dispersionPlot(genes(cuff))

dev.off()



#### FPKM squared coefficient of variation plot 

# The squared coefficient of variation allows visualization of cross-replicate variability between conditions 
# and can be a useful metric in determining data quality at the gene level (left) or isoform level (right).



pdf(file=paste(output_folder,cuffdiff_output_name,"/","stats_and_QC/","FPKM_SCV_plot_",cuffdiff_output_name,".pdf",sep=""))

fpkmSCVPlot(genes(cuff))
fpkmSCVPlot(isoforms(cuff))

dev.off()


#### Density plot replicate FPKM distributions

pdf(file=paste(output_folder,cuffdiff_output_name,"/","stats_and_QC/","densityPlot_",cuffdiff_output_name,".pdf",sep=""))

csDensity(genes(cuff),replicates=T)

dev.off()


#### Boxplots of FPKM distribution

pdf(file=paste(output_folder,cuffdiff_output_name,"/","stats_and_QC/","boxPlot_",cuffdiff_output_name,".pdf",sep=""))

csBoxplot(genes(cuff),replicates=T)

dev.off()



#### Scatterplot for Pairwise comparisons between samples

# Scatterplots can be useful to identify global changes and trends in gene expression between pairs of conditions.
#

pdf(file=paste(output_folder,cuffdiff_output_name,"/","stats_and_QC/","scatterPlot_",cuffdiff_output_name,".pdf",sep=""))

csScatterMatrix(genes(cuff))

dev.off()




#### Dendrogram of all samples

# Dendrogram with replicates can identify outlier replicates.
#

pdf(file=paste(output_folder,cuffdiff_output_name,"/","stats_and_QC/","dendroPlot_",cuffdiff_output_name,".pdf",sep=""))

csDendro(genes(cuff),replicates=T)

dev.off()







#### Volcano plots

# MvsA plots can be useful to determine any systematic bias that may be present between conditions.
#

pdf(file=paste(output_folder,cuffdiff_output_name,"/","stats_and_QC/","volcanoPlot_",cuffdiff_output_name,".pdf",sep=""))

csVolcanoMatrix(genes(cuff))

dev.off()





#############################################################
#															#
#		Accessing data and gene level plots					#
#															#
#############################################################


dir.create(paste(output_folder,cuffdiff_output_name,"/","gene_level_plots", sep=""), showWarnings = FALSE)


gene.featurenames<-featureNames(genes(cuff))
gene.matrix<-fpkmMatrix(genes(cuff))
gene.count.matrix<-countMatrix(genes(cuff))


singleString <- paste(readLines(gene_list_file), collapse=" ")
gene_list <- strsplit(singleString, ' ')



### handling the case where two gene names are in the cufflinks output, this will yield an error that the gene cannot be found. ###

for(g in gene_list[[1]])
{
	gene.features<-annotation(genes(cuff))		### get gene features and names as a dataframe from the cuff object
	grep_out<-grep(g, gene.features$gene_short_name)		### grep for the gene inside the table
	gene_name_inside_cufflinks <- gene.features[grep_out,]$gene_short_name		### return the gene name(s)
	gene_id_inside_cufflinks <- gene.features[grep_out,]$gene_id		### return the gene id 
	all_gene_names <- strsplit(gene_name_inside_cufflinks, ",")					### split by the separator to get all the gene names separated

	
	if(length(all_gene_names[[1]]) == 1)		## only one gene name can be found, take the gene as it is.
	{
			AddItemNaive(g)
	}
	else			### get the gene_id and retrieve the gene from cufflinks using the gene_id
	{	
		if(length(all_gene_names[[1]]) > 1) 
			{
			AddItemNaive(gene_id_inside_cufflinks)
			}
		else
		{
			stop(paste("-------------------- Are you sure that the gene ",g," is in database? --------------------",sep=""))
		}
	}
}
print(gene_ids)

myGenes<-getGenes(cuff,gene_ids)		### get the data for the genes that are listed in the gene_list_file.

### Expression Heat map for the selected set of genes that in gene_list_file

pdf(file=paste(output_folder,cuffdiff_output_name,"/","gene_level_plots/","ExpressionHeatmap_",cuffdiff_output_name,".pdf",sep=""))

csHeatmap(myGenes,cluster='both',replicates=T)

dev.off()



### Expression  Bar plot for the selected set of genes that in gene_list_file

pdf(file=paste(output_folder,cuffdiff_output_name,"/","gene_level_plots/","ExpressionBarPlot_",cuffdiff_output_name,".pdf",sep=""))

expressionBarplot(myGenes)

dev.off()


### A dendrogram of the relationship between conditions based on the expression of genes.

#
# Dendrograms can provide insight into the relationships between conditions for various genesets
# (e.g. significant genes used to draw relationships between conditions)
#

pdf(file=paste(output_folder,cuffdiff_output_name,"/","gene_level_plots/","GeneExpressionLevelDendrogram_",cuffdiff_output_name,".pdf",sep=""))

den<-csDendro(myGenes,replicates=T)

dev.off()




#############################################################
#															#
#					Individual gene plots					#
#															#
#############################################################

dir.create(paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots",sep=""), showWarnings = FALSE)

for(g in gene_list[[1]])
{
print(g)
dir.create(paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots/",g,sep=""), showWarnings = FALSE)		### create a folder for each gene

myGene<- NULL


### handling the case where two gene names are in the cufflinks output, this will yield an error that the gene cannoebe found. ###

gene.features<-annotation(genes(cuff))		### get gene features and names as a dataframe from the cuff object

grep_out<-grep(g, gene.features$gene_short_name)		### grep for the gene inside the table

gene_name_inside_cufflinks <- gene.features[grep_out,]$gene_short_name		### return the gene name(s)

gene_id_inside_cufflinks <- gene.features[grep_out,]$gene_id		### return the gene id 

all_gene_names <- strsplit(gene_name_inside_cufflinks, ",")					### split by the separator to get all the gene names separated

if(length(all_gene_names[[1]]) == 1)		## only one gene name can be found, take the gene as it is.
{
myGene<-getGene(cuff,g)
}
else			### get the gene_id and retrieve the gene from cufflinks using the gene_id
{	
		if(length(all_gene_names[[1]]) > 1) 
			{
			myGene <- getGene(cuff, gene_id_inside_cufflinks)
			}
		else
		{
			stop(paste("-------------------- Are you sure that the gene ",g," is in database? --------------------",sep=""))
		}
}


### checking if the given gene names are in the cuffdiff database
#
#check_gene_name_in_db <- try(getGene(cuff,g), silent=T)
#if(is(check_gene_name_in_db,"try-error"))
#{
#stop(paste("-------------------- Are you sure that the gene ",g," is in database? --------------------",sep=""))
#}
#else
#{
#myGene<-getGene(cuff,g)
#}


### Gene expression plot for individual genes

pdf(file=paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots/",g,"/ExpressionPlot_",cuffdiff_output_name,".pdf",sep=""))

print(expressionPlot(myGene,replicates=TRUE))

dev.off()


### Gene expression plot for all isoforms of individual genes

pdf(file=paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots/",g,"/IsoformsExpressionPlot_",cuffdiff_output_name,".pdf",sep=""))

print(expressionPlot(isoforms(myGene),replicates=T))

dev.off()



### Gene expression plot for the CDS of individual genes -- can only be used with Coding sequences, an error will be raised if applied for non-coding sequences

if(length(CDS(myGene))==0)
{
print(paste("Cannot generate CDS expression plot for gene ",g," because it is not a coding sequence",sep=""))
} else {
pdf(file=paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots/",g,"/CDSExpressionPlot_",cuffdiff_output_name,".pdf",sep=""))

print(expressionPlot(CDS(myGene),replicates=T))

dev.off()
}


### Gene expression Bar plot for individual genes

pdf(file=paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots/",g,"/BarExpressionPlot_",cuffdiff_output_name,".pdf",sep=""))

print(expressionBarplot(myGene,replicates=T))

dev.off()



### Gene expression Bar plot for all isoforms of individual genes

pdf(file=paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots/",g,"/BarExpressionPlot_",cuffdiff_output_name,".pdf",sep=""))

print(expressionBarplot(isoforms(myGene),replicates=T))

dev.off()




### Gene expression Bar plot for all CDS of individual genes

if(length(CDS(myGene))==0)
{
print(paste("Cannot generate CDS expression Bar plot for gene ",g," because it is not a coding sequence",sep=""))
} else {
pdf(file=paste(output_folder,cuffdiff_output_name,"/","Individual_genes_plots/",g,"/BarExpressionPlot_",cuffdiff_output_name,".pdf",sep=""))

print(expressionBarplot(CDS(myGene),replicates=T))

dev.off()
}

}

#############################################################
#															#
#					Distance matrix							#
#															#
#############################################################


## Heatmap of the pairwise similarities between conditions

#
# Similarities between conditions and/or replicates can provide useful insight into the relationship between various groupings of conditions and can aid in identifying outlier replicates that do not behave as expected.
#

pdf(file=paste(output_folder,cuffdiff_output_name,"/","DistHeatmapPlot_",cuffdiff_output_name,".pdf",sep=""))

csDistHeat(genes(cuff),replicates=T)

dev.off()




#############################################################
#				Dimensionality reduction					#
#			principal component analysis (PCA) 				#
#     					and 								#
#			multi-dimensional scaling (MDS)					#
#															#
#############################################################

# Dimensionality reduction is an informative approach for clustering and exploring the relationships between conditions.
# It can be useful for feature selection as well as identifying the sources of variability within your data.


### gene level features PCA plot


pdf(file=paste(output_folder,cuffdiff_output_name,"/","PCAPlot_",cuffdiff_output_name,".pdf",sep=""))

PCAplot(genes(cuff),"PC1","PC2",replicates=T)

dev.off()



### gene level features MDS plot


pdf(file=paste(output_folder,cuffdiff_output_name,"/","MDSPlot_",cuffdiff_output_name,".pdf",sep=""))

MDSplot(genes(cuff),replicates=T)

dev.off()



#############################################################
#															#
#					Significantly regulated genes			#
#															#
#############################################################

### Overview of significant features, This plot describes the number of significant genes at a 5% FDR for each pairwise interaction tested.

pdf(file=paste(output_folder,cuffdiff_output_name,"/","SignificantGenesOverview_",cuffdiff_output_name,".pdf",sep=""))

sigMatrix(cuff,level='genes',alpha=0.05)

dev.off()



mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')		### get tracking ids of significant genes that reject the null hypothesis in any condition tested
mySigGenes<-getGenes(cuff,mySigGeneIds)					### get a CuffGeneSet object of signigicantly regulated genes - this can be used for plotting



### Expression Heat map for the selected set of significantly regulated genes

pdf(file=paste(output_folder,cuffdiff_output_name,"/","TopRegulatedGenesExpressionHeatmap_",cuffdiff_output_name,".pdf",sep=""))

csHeatmap(mySigGenes,cluster='both',replicates=T,labRow=F)

dev.off()





