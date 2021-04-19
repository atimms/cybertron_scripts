library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")

workingDir = "/home/atimms/ngs_data/rnaseq/kim_rnaseq_0218";
setwd(workingDir);


##rl combined set
##read in count and metadata
countData1 <- read.table('kim_rnaseq_0218_rl_combined.star_fc.counts_adj.txt', header=T, row.names=1)
colData1 <- read.table('kim_rnaseq_0218_rl_combined.star_fc.metadata_adj.txt', header=T, row.names=1)
#countData1 <- read.table('kim_rnaseq_0218_rl_combined.star_fc.counts.txt', header=T, row.names=1)
#colData1 <- read.table('kim_rnaseq_0218_rl_combined.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
colData1$rnaaccess_batch <- as.factor(colData1$rnaaccess_batch)
colData1$age_pcw <- as.factor(colData1$age_pcw)
colData1$donor <- as.factor(colData1$donor)
#dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~rnaaccess_batch + sex + age_pcw + tissue)
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~rnaaccess_batch + age_pcw + tissue)
dds
##remove rows of the DESeqDataSet that have no counts, or less than 10 counts across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 15, ]
nrow(dds)
##rlog transform and check
rld <- vst(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="kim_rnaseq_0218_rl_combined.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$tissue, rld$donor, sep="-" )
#rownames(sampleDistMatrix) <- rownames(colData1)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='kim_rnaseq_0218_rl_combined.sample_heatmap.pdf', width = 7, height = 5)
#principal components analysis 
plotPCA(rld, intgroup = "tissue")
ggsave('kim_rnaseq_0218_rl_combined.pca.tissue.pdf', width=6, height = 6) 
plotPCA(rld, intgroup = "age_pcw")
ggsave('kim_rnaseq_0218_rl_combined.pca.age.pdf', width=6, height = 6) 
plotPCA(rld, intgroup = "donor")
ggsave('kim_rnaseq_0218_rl_combined.pca.donor.pdf', width=6, height = 6) 
plotPCA(rld, intgroup = "rnaaccess_batch")
ggsave('kim_rnaseq_0218_rl_combined.pca.batch.pdf', width=6, height = 6) 
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
##remove first 2 genes as giving strange signal in one sample
#topVarGenes = topVarGenes[3:25]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("tissue","sex","rnaaccess_batch", "age_pcw", "donor")])
pheatmap(mat, annotation_col=df, fontsize_row=5, fontsize_col=5)
#pheatmap(mat, annotation_col=df)
dev.copy2pdf(file='kim_rnaseq_0218_rl_combined.deseq.gene_clustering.pdf', width = 7, height = 5)
##differential expression
##do the test, need to add filtering to remove rows with low counts
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 10
dds <- dds[filter,]
#dds <- DESeq(dds)
##didn't work so use.. https://support.bioconductor.org/p/65091/
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=2000)

##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("tissue", "EGL", "PCL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.egl_vs_pcl.csv")
res2 <- results(dds, contrast=c("tissue", "RL", "PCL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.rl_vs_pcl.csv")
res2 <- results(dds, contrast=c("tissue", "RL", "EGL"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.rl_vs_egl.csv")

##new analysis 0819
res2 <- results(dds, contrast=c("tissue", "RL", "Bulk"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.rl_vs_bulk.csv")
res2 <- results(dds, contrast=c("tissue", "PCL", "Bulk"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.pcl_vs_bulk.csv")
res2 <- results(dds, contrast=c("tissue", "EGL", "Bulk"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:15836,]
write.csv(resOrdered2DF, file="kim_rnaseq_0218_rl_combined.de.egl_vs_bulk.csv")


