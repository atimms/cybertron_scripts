##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/sophie_splotch_1220";
setwd(workingDir);

###all data
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('sophie_splotch_1220.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sophie_splotch_1220.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$genotype <- as.factor(colData1$genotype)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~genotype)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="sophie_splotch_1220.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sophie_splotch_1220.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='sophie_splotch_1220.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('sophie_splotch_1220.genotype_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_name", "genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='sophie_splotch_1220.25_var_gene_clustering.pdf', width = 9, height = 6)




###removed 2 samples and added sex to teh metadata
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('sophie_splotch_2_1220.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sophie_splotch_2_1220.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)


##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~sex + genotype)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="sophie_splotch_2_1220.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sophie_splotch_2_1220.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='sophie_splotch_2_1220.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('sophie_splotch_2_1220.genotype_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_name", "genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='sophie_splotch_2_1220.25_var_gene_clustering.pdf', width = 9, height = 6)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "mut", "wt"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18788,]
write.csv(resOrdered2DF, file="sophie_splotch_2_1220.mut_vs_wt.csv")






