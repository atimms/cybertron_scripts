##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321";
setwd(workingDir);

###cunn_flexcell_rnaseq_0321
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Gene, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0321.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.25_var_gene_clustering.pdf', width = 10, height = 10)


###cunn_flexcell_rnaseq_0321.CTRL
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_CTRL.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_CTRL.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.CTRL.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.CTRL.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.CTRL.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.CTRL.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19911,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.CTRL.Strained_vs_Unstrained.csv")




###cunn_flexcell_rnaseq_0321.FLNA
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNA.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNA.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.FLNA.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.FLNA.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNA.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNA.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19885,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.FLNA.Strained_vs_Unstrained.csv")


###cunn_flexcell_rnaseq_0321.FLNB
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNB.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNB.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.FLNB.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.FLNB.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNB.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNB.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19444,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.FLNB.Strained_vs_Unstrained.csv")


###cunn_flexcell_rnaseq_0321.FLNC
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNC.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNC.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.FLNC.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.FLNC.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNC.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNC.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19586,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.FLNC.Strained_vs_Unstrained.csv")


###cunn_flexcell_rnaseq_0321.PIEZO1
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_PIEZO1.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_PIEZO1.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.PIEZO1.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.PIEZO1.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.PIEZO1.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.PIEZO1.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20253,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.PIEZO1.Strained_vs_Unstrained.csv")

