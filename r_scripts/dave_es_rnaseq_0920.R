##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/dave_es_cell_rnaseq_0920";
setwd(workingDir);

###es_cells_0920
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('dave_esc_0920.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('dave_esc_0920.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$genotype <- as.factor(colData1$genotype)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~time_point + genotype)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="es_cells_0920.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="es_cells_0920.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$time_point, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='es_cells_0920.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("time_point"))
ggsave('es_cells_0920.timepoint_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("genotype"))
ggsave('es_cells_0920.genotype_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("time_point", "genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='es_cells_0920.25_var_gene_clustering.pdf', width = 9, height = 6)

########################

###es_cells_0920_day0
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('dave_esc_d0_0920.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('dave_esc_d0_0920.star_fc.metadata.txt', header=T, row.names=1)
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
#write.csv(count_data, file="es_cells_0920_day0.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="es_cells_0920_day0.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$time_point, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='es_cells_0920_day0.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('es_cells_0920_day0.genotype_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='es_cells_0920_day0.25_var_gene_clustering.pdf', width = 9, height = 6)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "C7", "F1"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:17193,]
write.csv(resOrdered2DF, file="es_cells_0920_day0.C7_vs_F1.csv")


###es_cells_0920_day4
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('dave_esc_d4_0920.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('dave_esc_d4_0920.star_fc.metadata.txt', header=T, row.names=1)
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
#write.csv(count_data, file="es_cells_0920_day4.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="es_cells_0920_day4.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$time_point, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='es_cells_0920_day4.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('es_cells_0920_day4.genotype_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='es_cells_0920_day4.25_var_gene_clustering.pdf', width = 9, height = 6)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "C7", "F1"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:17998,]
write.csv(resOrdered2DF, file="es_cells_0920_day4.C7_vs_F1.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "419", "F1"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:17998,]
write.csv(resOrdered2DF, file="es_cells_0920_day4.419_vs_F1.csv")



###es_cells_0920_day8
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('dave_esc_d8_0920.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('dave_esc_d8_0920.star_fc.metadata.txt', header=T, row.names=1)
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
#write.csv(count_data, file="es_cells_0920_day8.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="es_cells_0920_day8.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$time_point, rld$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='es_cells_0920_day8.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("genotype"))
ggsave('es_cells_0920_day8.genotype_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("genotype")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='es_cells_0920_day8.25_var_gene_clustering.pdf', width = 9, height = 6)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "C7", "F1"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:18254,]
write.csv(resOrdered2DF, file="es_cells_0920_day8.C7_vs_F1.csv")


