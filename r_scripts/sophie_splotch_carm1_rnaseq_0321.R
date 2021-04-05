##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##install sva for combat-seq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("sva")
##need to use R 4 plus for combat-seq

##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("sva") ##for combat-seq

workingDir = "/home/atimms/ngs_data/rnaseq/sophie_carm1_splotch_rnaseq_0321";
setwd(workingDir);

###all data
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('sophie_carm1_splotch_rnaseq_0321.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('sophie_carm1_splotch_rnaseq_0321.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
#colData1$genotype <- as.factor(colData1$genotype)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, sex~genotype)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="sophie_carm1_splotch_rnaseq_0321.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sophie_carm1_splotch_rnaseq_0321.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$sex, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='sophie_carm1_splotch_rnaseq_0321.no_correction.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('sophie_carm1_splotch_rnaseq_0321.no_correction.genotype_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sex"))
ggsave('sophie_carm1_splotch_rnaseq_0321.no_correction.sex_pca.pdf', width=6, height = 6)


##lets run combat-seq
#https://github.com/zhangyuqing/ComBat-seq
#example here https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/
batches = c(1,1,1,1,1,1,2,2,2,2,2,2)
genotype_group = c(1,1,1,2,2,2,1,2,1,1,2,2)
sex_group = c(1,2,2,1,1,2,2,2,2,2,2,1)
covariate_matrix = cbind(genotype_group, sex_group)
##? remove 1st column, use covar_mod
#corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data[,sample_names]), batch = batches, group = groups)
corrected_data = ComBat_seq(counts = as.matrix(countData1), batch = batches, group = genotype_group)
#corrected_data = ComBat_seq(counts = as.matrix(countData1), batch = batches, covar_mod = covariate_matrix)
head(corrected_data)
head(countData1)

###using genotype i.e. keep mice seperate

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = corrected_data, colData = colData1, sex~genotype)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="sophie_carm1_splotch_rnaseq_0321.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sophie_carm1_splotch_rnaseq_0321.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$sex, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='sophie_carm1_splotch_rnaseq_0321.ComBat.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('sophie_carm1_splotch_rnaseq_0321.ComBat.genotype_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sex"))
ggsave('sophie_carm1_splotch_rnaseq_0321.ComBat.sex_pca.pdf', width=6, height = 6)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "het", "splotch"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19309,]
write.csv(resOrdered2DF, file="sophie_splotch_0121.ComBat.het_vs_splotch.csv")
res2 <- results(dds, contrast=c("genotype", "wt", "carm1"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19309,]
write.csv(resOrdered2DF, file="sophie_splotch_0121.ComBat.wt_vs_carm1.csv")


###using genotype_simple i.e. just mut and control

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = corrected_data, colData = colData1,  ~ genotypesimple)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="sophie_carm1_splotch_rnaseq_0321.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="sophie_carm1_splotch_rnaseq_0321.deseq.vst_counts.csv")

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))

##to get a specific test:
res2 <- results(dds, contrast=c("genotypesimple", "ctl", "mut"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19309,]
write.csv(resOrdered2DF, file="sophie_splotch_0121.ComBat.ctl_vs_combined_mut.csv")


