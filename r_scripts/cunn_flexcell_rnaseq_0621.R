##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
33using r 3.6

##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0621";
setwd(workingDir);

###cunn_flexcell_rnaseq_0321 --- all samples
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0621.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0621.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~sample_only + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 12, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0621.norm_counts.csv")
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
dev.copy2pdf(file='cunn_flexcell_rnaseq_0621.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0621.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0621.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0621.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0621.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0621.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0621.25_var_gene_clustering.pdf', width = 10, height = 10)






###cunn_flexcell_rnaseq_0321 --- just samples with duplicates
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_dups_0621.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_dups_0621.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Cohort + sample_condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 12, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_dups_0621.norm_counts.csv")
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
dev.copy2pdf(file='cunn_flexcell_rnaseq_dups_0621.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_dups_0621.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_dups_0621.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_dups_0621.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_dups_0621.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_dups_0621.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_dups_0621.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("sample_condition", "C1653_Strained", "C1653_Unstrained"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21019,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0621.C1653.Strained_vs_Unstrained.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("sample_condition", "C2012_Strained", "C2012_Unstrained"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21019,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0621.C2012.Strained_vs_Unstrained.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("sample_condition", "C2013_Strained", "C2013_Unstrained"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21019,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0621.C2013.Strained_vs_Unstrained.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("sample_condition", "C3049_Strained", "C3049_Unstrained"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21019,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0621.C3049.Strained_vs_Unstrained.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("sample_condition", "C3066_Strained", "C3066_Unstrained"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21019,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0621.C3066.Strained_vs_Unstrained.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("sample_condition", "C4019_Strained", "C4019_Unstrained"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21019,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0621.C4019.Strained_vs_Unstrained.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("sample_condition", "CTL_Strained", "CTL_Unstrained"))
##get summary
summary(res2) #lots of significant genes
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:21019,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0621.CTL.Strained_vs_Unstrained.csv")

