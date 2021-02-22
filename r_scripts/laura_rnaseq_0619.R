##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/laura_rnaseq_0619";
setwd(workingDir);

###all data
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('laura_rnaseq_0619.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('laura_rnaseq_0619.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)


##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~batch + dx)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="laura_rnaseq_0619.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="laura_rnaseq_0619.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$dx ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='laura_rnaseq_0619.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("batch"))
ggsave('laura_rnaseq_0619.batch_pca.pdf', width=6, height = 6)
##principal components analysis 
plotPCA(rld, intgroup = c("dx"))
ggsave('laura_rnaseq_0619.dx_pca.pdf', width=6, height = 6)


##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_name", "batch", "dx")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='laura_rnaseq_0619.25_var_gene_clustering.pdf', width = 9, height = 6)
