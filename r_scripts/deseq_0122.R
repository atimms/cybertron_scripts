##load libraries
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/tmp/rnaseq_0122/deseq";
setwd(workingDir);

##read in count data
countData1 <- read.table('splotch.star_fc.counts.txt', header=T, row.names=1)
head(countData1)
##read in count and metadata -- using counts using subset of genes
colData1 <- read.table('splotch.star_fc.metadata.txt', header=T, row.names=1)
head(colData1)


##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~sex + genotype)
dds

##remove rows of the DESeqDataSet that have less than 10 reads total
nrow(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds)

##vst transform data
vst <- vst(dds, blind=TRUE)
##rlog transform data
rld <- rlog(dds, blind=TRUE)

##check
head(assay(dds), 3)
head(assay(vst), 3)
head(assay(rld), 3)

##calculate sample distances from rlog 
sampleDists <- dist(t( assay(vst)))
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )      
rownames(sampleDistMatrix) <- paste(vst$sample_name, vst$genotype ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='splotch.sample_heatmap.pdf', width = 7, height = 5)

##principal components analysis 
plotPCA(vst, intgroup = c("genotype"))
ggsave('splotch.vst.genotype_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("genotype"))
ggsave('splotch.rld.genotype_pca.pdf', width=6, height = 6)
plotPCA(vst, intgroup = c("sex"))
ggsave('splotch.vst.sex_pca.pdf', width=6, height = 6)
plotPCA(vst, intgroup = c("sex", "genotype"))
ggsave('splotch.vst.sex_genotype_pca.pdf', width=6, height = 6)


##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(vst)),decreasing=TRUE),25)
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(vst)[c("sample_name", "genotype", "sex")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='splotch.25_var_gene_clustering.pdf', width = 9, height = 6)


##size factors and normalized counts
ddsSF <- estimateSizeFactors(dds)
sizeFactors(ddsSF)
head(counts(ddsSF, normalized = FALSE))
head(counts(ddsSF, normalized = TRUE))

##estimate variance
ddsLocal <- estimateDispersions(ddsSF, fitType = "local") 
ddsMean <- estimateDispersions(ddsSF, fitType = "mean") 
ddsParam <- estimateDispersions(ddsSF, fitType = "parametric")
op <- par(mfrow = c(1, 3))
plotDispEsts(ddsMean) 
plotDispEsts(ddsParam) 
plotDispEsts(ddsLocal)
par(op)

 
##Performing gene-wise comparisons
ddsParam <- nbinomWaldTest(ddsParam) 
res <- results(ddsParam)
res

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary, this just gets the last test
res <- results(dds)
summary(res)
##to get a specific test:
res2 <- results(dds, contrast=c("genotype", "ko", "wt"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:17620,]
write.csv(resOrdered2DF, file="splotch.mut_vs_het.csv")
