##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/cunn_rnaseq_0220";
setwd(workingDir);

##all samples
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('strain_rnaseq_0220_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('strain_rnaseq_0220_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
#colData1$Age_GW <- as.factor(colData1$Age_GW)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~sample_only + strain)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="strain_rnaseq_0220_all.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="kim_rnaseq_0720_b6_bulk.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$strain, rld$sample_only ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='strain_rnaseq_0220_all.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("gene"))
ggsave('strain_rnaseq_0220_all.gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("strain"))
ggsave('strain_rnaseq_0220_all.strain_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sample_only"))
ggsave('strain_rnaseq_0220_all.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sex"))
ggsave('strain_rnaseq_0220_all.sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("batch"))
ggsave('strain_rnaseq_0220_all.batch_pca.pdf', width=6, height = 6)

