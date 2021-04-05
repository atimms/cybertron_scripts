##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("sva") ##for combat-seq

workingDir = "/home/atimms/ngs_data/rnaseq/laura_rnaseq_0619";
setwd(workingDir);

###all data together
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('laura_rnaseq_0619.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('laura_rnaseq_0619.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)


##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~batch + dx)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="laura_rnaseq_0619.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="laura_rnaseq_0619.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$sample_name, rld$dx ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='laura_rnaseq_0619.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("batch"))
ggsave('laura_rnaseq_0619.batch_pca.pdf', width=6, height = 6)
##principal components analysis 
plotPCA(rld, intgroup = c("dx"))
ggsave('laura_rnaseq_0619.dx_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_name", "batch", "dx")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='laura_rnaseq_0619.25_var_gene_clustering.pdf', width = 9, height = 6)

###Dysplasia datase only
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('laura_rnaseq_dys_0619.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('laura_rnaseq_dys_0619.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~dx)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="laura_rnaseq_dys_0619.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="laura_rnaseq_dys_0619.deseq.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('laura_rnaseq_dys_0619.sample_name_pca.pdf', width=6, height = 6)
##principal components analysis 
plotPCA(rld, intgroup = c("dx"))
ggsave('laura_rnaseq_dys_0619.dx_pca.pdf', width=6, height = 6)
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_name", "dx")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='laura_rnaseq_dys_0619.25_var_gene_clustering.pdf', width = 9, height = 6)

##clustering/heatmap from list
genes = c('Vps37a', 'Gm9573', 'Adcyap1r1', 'Zfyve26', 'Krt4', 'Fat1', 'Neat1', 'Vprbp', 'Gsto1', 'Muc6', 'S100a11', 'Rian', 'Zfp36l1', 'Auts2', 'Bdp1', 'Ccnd2', 'Etl4', 'Pofut1', 'Slc23a3', 'Trap1', 'Siah3', 'Dlgap1', 'Dnajc6', 'Tagln3', 'Akt1s1', 'Pglyrp4', 'Ephb6', 'Slc15a3', 'Ppm1f', 'Slc45a4', 'Cul1', 'Ccnd3', 'Mdm2', 'Cdkn1a', 'Tnfrsf1a', 'Muc4', 'Cd44', 'Dsp', 'Tmem173', 'Kif5b', 'Dusp9', 'Runx1', 'Col27a1', 'Dsc2', 'Tnrc6b', 'Chd3', 'Lama3', 'Cdk16', 'Ly6g6c')
#genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("sample_name","dx")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='laura_rnaseq_dys_0619.dys_gene_clustering.pdf', width = 9, height = 6)

##differential expression
##do the test
dds <- DESeq(dds)
##get a specific test:
res2 <- results(dds, contrast=c("dx", "control", "dysplasia"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:12022,]
write.csv(resOrdered2DF, file="laura_rnaseq_dys_0619.ctl_vs_dys_mut.csv")

###SCC datase only
##read in count and metadata -- using counts using subset of genes
countData1 <- read.table('laura_rnaseq_scc_0619.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('laura_rnaseq_scc_0619.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~dx)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="laura_rnaseq_scc_0619.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="laura_rnaseq_scc_0619.deseq.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('laura_rnaseq_scc_0619.sample_name_pca.pdf', width=6, height = 6)
##principal components analysis 
plotPCA(rld, intgroup = c("dx"))
ggsave('laura_rnaseq_scc_0619.dx_pca.pdf', width=6, height = 6)
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_name", "dx")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='laura_rnaseq_scc_0619.25_var_gene_clustering.pdf', width = 9, height = 6)

##clustering/heatmap from list
genes = c('Clca3b', 'Dsg3', 'Krt6a', 'Krt13', 'Serpinb2', 'Claml3', 'Col17a1', 'Tmprss11b', 'Dsc3', 'Krt17', 'Psca', 'Krt14', 'Krt5', 'Dsc2', 'Pkp1', 'Arg1', 'Krt15', 'Ltf', 'Gm9573', 'Gsto1', 'Snord34', 'Bpifa1', 'Fn1', 'Sfi1', 'Snord12', 'Snord32a', 'Snord58b', 'Snord99', 'Anxa1', 'Eef2', 'Fat1', 'Kcnq1ot1', 'Krt4', 'Ly6g6c', 'Muc4', 'Pisd-ps1', 'Snora74a', 'Snord104', 'Snord49a', 'Dsp', 'Itgb1', 'S100a11', 'Cd44', 'Anxa8')
#genes <- genes$genes #vectorize
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[,c("sample_name","dx")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='laura_rnaseq_scc_0619.scc_gene_clustering.pdf', width = 9, height = 6)

##differential expression
##do the test
dds <- DESeq(dds)
##get a specific test:
res2 <- results(dds, contrast=c("dx", "control", "scc"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22562,]
write.csv(resOrdered2DF, file="laura_rnaseq_scc_0619.ctl_vs_scc_mut.csv")



##apply combat...
countData1 <- read.table('laura_rnaseq_0619.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('laura_rnaseq_0619.star_fc.metadata.txt', header=T, row.names=1)
##lets run combat-seq
#https://github.com/zhangyuqing/ComBat-seq
#example here https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/
batches = c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2)
groups = c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
corrected_data = ComBat_seq(counts = as.matrix(countData1), batch = batches, group = groups)
head(corrected_data)
head(countData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = corrected_data, colData = colData1, ~dx)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
#write.csv(count_data, file="laura_rnaseq_combat_0619.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="laura_rnaseq_combat_0619.deseq.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("sample_name"))
ggsave('laura_rnaseq_combat_0619.sample_name_pca.pdf', width=6, height = 6)
##principal components analysis 
plotPCA(rld, intgroup = c("dx"))
ggsave('laura_rnaseq_combat_0619.dx_pca.pdf', width=6, height = 6)
##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("sample_name", "dx")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='laura_rnaseq_combat_0619.25_var_gene_clustering.pdf', width = 9, height = 6)





