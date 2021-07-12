##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )

##using R 4.03

##load libarys 
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/vishal_rnaseq_0621";
setwd(workingDir);

###vishal_rnaseq_0621
##read in count and metadata
countData1 <- read.table('vishal_rnaseq_0621.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_rnaseq_0621.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$condition <- as.factor(colData1$condition)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="vishal_rnaseq_0621.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="vishal_rnaseq_0621.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$condition ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='vishal_rnaseq_0621.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('vishal_rnaseq_0621.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("condition"))
ggsave('vishal_rnaseq_0621.condition_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("condition", "sample_name")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='vishal_rnaseq_0621.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "group1", "control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19587,]
write.csv(resOrdered2DF, file="vishal_rnaseq_0621.group1_vs_control.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("condition", "group2", "control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19587,]
write.csv(resOrdered2DF, file="vishal_rnaseq_0621.group2_vs_control.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("condition", "group3", "control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19587,]
write.csv(resOrdered2DF, file="vishal_rnaseq_0621.group3_vs_control.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("condition", "group2", "group3"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19587,]
write.csv(resOrdered2DF, file="vishal_rnaseq_0621.group2_vs_group3.csv")

###vishal_rnaseq_0621_ctl_g1
##read in count and metadata
countData1 <- read.table('vishal_rnaseq_0621_ctl_g1.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_rnaseq_0621_ctl_g1.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$condition <- as.factor(colData1$condition)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~condition)
dds
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="vishal_rnaseq_0621_ctl_g1.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="vishal_rnaseq_0621_ctl_g1.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$condition ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='vishal_rnaseq_0621_ctl_g1.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('vishal_rnaseq_0621_ctl_g1.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("condition"))
ggsave('vishal_rnaseq_0621_ctl_g1.condition_pca.pdf', width=6, height = 6)
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("condition", "sample_name")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='vishal_rnaseq_0621_ctl_g1.25_var_gene_clustering.pdf', width = 10, height = 10)
##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "group1", "control"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18275,]
write.csv(resOrdered2DF, file="vishal_rnaseq_0621_ctl_g1.group1_vs_control.csv")



###vishal_rnaseq_0621_g2_g3
##read in count and metadata
countData1 <- read.table('vishal_rnaseq_0621_g2_g3.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('vishal_rnaseq_0621_g2_g3.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$condition <- as.factor(colData1$condition)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="vishal_rnaseq_0621_g2_g3.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="vishal_rnaseq_0621_g2_g3.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$condition ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='vishal_rnaseq_0621_g2_g3.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('vishal_rnaseq_0621_g2_g3.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("condition"))
ggsave('vishal_rnaseq_0621_g2_g3.condition_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("condition", "sample_name")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='vishal_rnaseq_0621_g2_g3.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "group2", "group3"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19587,]
write.csv(resOrdered2DF, file="vishal_rnaseq_0621_g2_g3.group2_vs_group3.csv")
