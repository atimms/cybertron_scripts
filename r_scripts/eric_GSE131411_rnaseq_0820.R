##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/eric_rnaseq_0720/GSE131411";
setwd(workingDir);

###GSE131411_ss_hom_all -- so all homozygote SepticShock patients
##read in count and metadata
countData1 <- read.table('GSE131411_ss_hom_all.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('GSE131411_ss_hom_all.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$age <- as.factor(colData1$age)

##add to deseq, give countdata and metadata and then design information with variable we're using listed last
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ timepoint + sex + age + rs3184504)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="GSE131411_ss_hom_all.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="GSE131411_ss_hom_all.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$patient ,rld$timepoint, rld$rs3184504,   sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='GSE131411_ss_hom_all.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("patient"))
ggsave('GSE131411_ss_hom_all.patient_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("timepoint"))
ggsave('GSE131411_ss_hom_all.timepoint_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("rs3184504"))
ggsave('GSE131411_ss_hom_all.rs3184504_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sex"))
ggsave('GSE131411_ss_hom_all.sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("age"))
ggsave('GSE131411_ss_hom_all.age_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("mortality"))
ggsave('GSE131411_ss_hom_all.mortality_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("patient", "timepoint", "rs3184504", "sex", "age", "mortality")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='GSE131411_ss_hom_all.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("rs3184504", "CC", "TT"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:26503,]
write.csv(resOrdered2DF, file="GSE131411_ss_hom_all.dwm_vs_ctl.csv")

###GSE131411_ss_hom_t1
##read in count and metadata
countData1 <- read.table('GSE131411_ss_hom_t1.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('GSE131411_ss_hom_t1.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$age <- as.factor(colData1$age)

##add to deseq, give countdata and metadata and then design information with variable we're using listed last
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sex + rs3184504)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="GSE131411_ss_hom_t1.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="GSE131411_ss_hom_t1.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$patient ,rld$timepoint, rld$rs3184504,   sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='GSE131411_ss_hom_t1.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("patient"))
ggsave('GSE131411_ss_hom_t1.patient_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("rs3184504"))
ggsave('GSE131411_ss_hom_t1.rs3184504_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sex"))
ggsave('GSE131411_ss_hom_t1.sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("age"))
ggsave('GSE131411_ss_hom_t1.age_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("mortality"))
ggsave('GSE131411_ss_hom_t1.mortality_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("patient", "rs3184504", "sex", "age", "mortality")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='GSE131411_ss_hom_t1.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("rs3184504", "CC", "TT"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22636,]
write.csv(resOrdered2DF, file="GSE131411_ss_hom_t1.dwm_vs_ctl.csv")

###GSE131411_ss_hom_t2
##read in count and metadata
countData1 <- read.table('GSE131411_ss_hom_t2.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('GSE131411_ss_hom_t2.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$age <- as.factor(colData1$age)

##add to deseq, give countdata and metadata and then design information with variable we're using listed last
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sex + rs3184504)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="GSE131411_ss_hom_t2.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="GSE131411_ss_hom_t2.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$patient ,rld$timepoint, rld$rs3184504,   sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='GSE131411_ss_hom_t2.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("patient"))
ggsave('GSE131411_ss_hom_t2.patient_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("rs3184504"))
ggsave('GSE131411_ss_hom_t2.rs3184504_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sex"))
ggsave('GSE131411_ss_hom_t2.sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("age"))
ggsave('GSE131411_ss_hom_t2.age_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("mortality"))
ggsave('GSE131411_ss_hom_t2.mortality_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("patient", "rs3184504", "sex", "age", "mortality")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='GSE131411_ss_hom_t2.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("rs3184504", "CC", "TT"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22732,]
write.csv(resOrdered2DF, file="GSE131411_ss_hom_t2.dwm_vs_ctl.csv")


###GSE131411_ss_hom_t3
##read in count and metadata
countData1 <- read.table('GSE131411_ss_hom_t3.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('GSE131411_ss_hom_t3.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$age <- as.factor(colData1$age)

##add to deseq, give countdata and metadata and then design information with variable we're using listed last
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sex + rs3184504)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="GSE131411_ss_hom_t3.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="GSE131411_ss_hom_t3.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$patient ,rld$timepoint, rld$rs3184504,   sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='GSE131411_ss_hom_t3.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("patient"))
ggsave('GSE131411_ss_hom_t3.patient_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("rs3184504"))
ggsave('GSE131411_ss_hom_t3.rs3184504_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("sex"))
ggsave('GSE131411_ss_hom_t3.sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("age"))
ggsave('GSE131411_ss_hom_t3.age_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("mortality"))
ggsave('GSE131411_ss_hom_t3.mortality_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("patient", "rs3184504", "sex", "age", "mortality")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='GSE131411_ss_hom_t3.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("rs3184504", "CC", "TT"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:24472,]
write.csv(resOrdered2DF, file="GSE131411_ss_hom_t3.dwm_vs_ctl.csv")
