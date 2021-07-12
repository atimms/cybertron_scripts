##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/gumbiner_rnaseq_0421";
setwd(workingDir);

###all samples
##read in count and metadata
countData1 <- read.table('gumbiner_rnaseq_0421_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('gumbiner_rnaseq_0421_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
#colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~time_point + condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="gumbiner_rnaseq_0421_all.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="gumbiner_rnaseq_0421_all.deseq.vst_counts.csv")
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
dev.copy2pdf(file='gumbiner_rnaseq_0421_all.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('gumbiner_rnaseq_0421_all.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("condition"))
ggsave('gumbiner_rnaseq_0421_all.condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("time_point"))
ggsave('gumbiner_rnaseq_0421_all.time_point_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_all.25_var_gene_clustering.pdf', width = 10, height = 10)

##make heatmap of genes
genes = c("EPHB3", "TM4SF1", "ACTG1", "TXNIP", "NDUFB9", "SLC25A25", "PROX1", "FNIP1", "SLC28A2", "NCR3LG1")
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_all.10_genes_clustering.pdf', width = 9, height = 6)




###4 hours
##read in count and metadata
countData1 <- read.table('gumbiner_rnaseq_0421_4h.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('gumbiner_rnaseq_0421_4h.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
#colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="gumbiner_rnaseq_0421_4h.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="gumbiner_rnaseq_0421_4h.deseq.vst_counts.csv")
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
dev.copy2pdf(file='gumbiner_rnaseq_0421_4h.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('gumbiner_rnaseq_0421_4h.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("condition"))
ggsave('gumbiner_rnaseq_0421_4h.condition_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_4h.25_var_gene_clustering.pdf', width = 10, height = 10)

##make heatmap of genes
genes = c("EPHB3", "TM4SF1", "ACTG1", "TXNIP", "NDUFB9", "SLC25A25", "PROX1", "FNIP1", "SLC28A2", "NCR3LG1")
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_4h.10_genes_clustering.pdf', width = 9, height = 6)

##make heatmap of DE genes
genes = c("EPHB3", "TXNIP", "PROX1", "TM4SF1", "ACTG1", "RP11-395P17.3", "SLC38A2", "IGF1R", "NCR3LG1", "FNIP1", "MDM4", "RAPH1", "ZMYND19", "TNRC6B", "POMP", "RPS28", "KLF2", "SKI", "RPL24", "NPC2", "INADL", "MAN2C1", "RPS16", "BDP1", "RPL18", "TMEM131", "RBM33", "EFHD2", "RPL13", "RSRC1")
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("condition")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=6, fontsize_col=8, cluster_cols = F)
dev.copy2pdf(file='gumbiner_rnaseq_0421_4h.30_DE_genes.clustering.pdf', width = 9, height = 6)



##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "treated", "untreated"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18872,]
write.csv(resOrdered2DF, file="gumbiner_rnaseq_0421_4h.treated_vs_untreated.csv")




###4 hours
##read in count and metadata
countData1 <- read.table('gumbiner_rnaseq_0421_24h.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('gumbiner_rnaseq_0421_24h.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
#colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="gumbiner_rnaseq_0421_24h.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="gumbiner_rnaseq_0421_24h.deseq.vst_counts.csv")
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
dev.copy2pdf(file='gumbiner_rnaseq_0421_24h.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('gumbiner_rnaseq_0421_24h.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("condition"))
ggsave('gumbiner_rnaseq_0421_24h.condition_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_24h.25_var_gene_clustering.pdf', width = 10, height = 10)

##make heatmap of genes
genes = c("EPHB3", "TM4SF1", "ACTG1", "TXNIP", "NDUFB9", "SLC25A25", "PROX1", "FNIP1", "SLC28A2", "NCR3LG1")
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_24h.10_genes_clustering.pdf', width = 9, height = 6)

##make heatmap of DE genes
genes = c("ACSS2", "CD24", "MYO1A", "SPCS3", "SLC25A1", "CLDN15", "ZDHHC20", "ITGA3", "SLC38A1", "HSPH1", "ARPP19", "SPTAN1", "EPHB3", "PROM2", "DNAJC3", "FASN", "CD68", "PRKCSH", "COPB1", "MUC20", "CDT1", "MSI2", "MVD", "PROX1", "MAP2K2", "CS", "TRIB3", "IQGAP3", "XPO1", "SF3B1")
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("condition")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=6, fontsize_col=8, cluster_cols = F)
dev.copy2pdf(file='gumbiner_rnaseq_0421_24h.30_DE_genes.clustering.pdf', width = 9, height = 6)



##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "treated", "untreated"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18023,]
write.csv(resOrdered2DF, file="gumbiner_rnaseq_0421_24h.treated_vs_untreated.csv")




###30 minutes
##read in count and metadata
countData1 <- read.table('gumbiner_rnaseq_0421_30m.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('gumbiner_rnaseq_0421_30m.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
#colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="gumbiner_rnaseq_0421_30m.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="gumbiner_rnaseq_0421_30m.deseq.vst_counts.csv")
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
dev.copy2pdf(file='gumbiner_rnaseq_0421_30m.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('gumbiner_rnaseq_0421_30m.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("condition"))
ggsave('gumbiner_rnaseq_0421_30m.condition_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_30m.25_var_gene_clustering.pdf', width = 10, height = 10)

##make heatmap of genes
genes = c("EPHB3", "TM4SF1", "ACTG1", "TXNIP", "NDUFB9", "SLC25A25", "PROX1", "FNIP1", "SLC28A2", "NCR3LG1")
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("condition", "time_point")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='gumbiner_rnaseq_0421_30m.10_genes_clustering.pdf', width = 9, height = 6)

##make heatmap of DE genes
genes = c("NDRG1", "E2F8", "FABP1", "SUV420H2", "SRXN1", "TFF3", "WASH7P", "MIR6859-1", "LOC101927589", "LOC729737", "LOC100996442", "LOC102723897", "LOC102723917", "RP4-669L17.10", "LOC100134822", "MIR6723", "LOC100133331", "LOC101928706", "LOC101060494", "LOC100288069", "LOC100287934", "FAM87B", "LINC00115", "LINC01128", "FAM41C", "LOC284600", "LOC100130417", "LOC101928801", "SAMD11", "NOC2L")
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("condition")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=6, fontsize_col=8, cluster_cols = F)
dev.copy2pdf(file='gumbiner_rnaseq_0421_30m.30_DE_genes.clustering.pdf', width = 9, height = 6)



##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("condition", "treated", "untreated"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18678,]
write.csv(resOrdered2DF, file="gumbiner_rnaseq_0421_30m.treated_vs_untreated.csv")




