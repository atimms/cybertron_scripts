# Load the ggplot2 package
library(ggplot2)
library(plyr)
library("DESeq2")
library(pheatmap)
##set working directory
setwd('/home/atimms/ngs_data/misc/eric_graph_SH2B_0920')


##Human Bone marrow, from noramalized data
dat <- read.csv('human_bm.SH2B_formatted.txt', sep = '\t')

##get the value for length, mean, and sd data, doesn't work in nans
##from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
cdata <- ddply(dat, c("tissue", "gene"), summarise,
               N    = length(value), mean = mean(value),
               sd   = sd(value), se   = sd / sqrt(N))
cdata

##graph with all 3 genes in one graph
ggplot(data=cdata, aes(x=tissue, y=mean, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3, width=.2, position=position_dodge(.9)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7))
ggsave('human_bm.SH2B_combined.pdf', width = 10, height = 6)
##graph with all 3 genes in seperate graphs i.e. facetted as bargraph
ggplot(data=cdata, aes(x=tissue, y=mean, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3, width=.2, position=position_dodge(.9)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('human_bm.SH2B.bar_graph.pdf', width = 10, height = 8)
##facetted boxplot
ggplot(dat, aes(x=tissue, y=value, fill=gene)) + geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('human_bm.SH2B.boxplot.pdf', width = 10, height = 8)



##mouse bone marrow -- GSE151682
##manually change raw data, remove " and .bam and columns about gene

##read in count and metadata
countData1 <- read.table('GSE151682_genecounts_20200601_formatted.txt', header=T, row.names=1)
colData1 <- read.table('GSE151682_metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ group)
dds
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="GSE151682_norm_counts.csv")

##graph reformatted norm counts for Sh1b genes
dat <- read.csv('GSE151682_for_graphing.txt', sep = '\t')

##get the value for length, mean, and sd data, doesn't work in nans
##from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
cdata <- ddply(dat, c("group", "gene"), summarise,
               N    = length(value), mean = mean(value),
               sd   = sd(value), se   = sd / sqrt(N))
cdata

##graph with all 3 genes in seperate graphs i.e. facetted
ggplot(data=cdata, aes(x=group, y=mean, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3, width=.2, position=position_dodge(.9)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =12)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('mouse_bm.GSE151682.SH2B.bar_graph.pdf', width = 10, height = 8)
##facetted boxplot
ggplot(dat, aes(x=group, y=value, fill=gene)) + geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('mouse_bm.GSE151682.SH2B.boxplot.pdf', width = 10, height = 8)

##DE for different groups
##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 16, ]
nrow(dds)
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
write.csv(assay(rld), file="GSE151682_vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("group"))
ggsave('mouse_bm.GSE151682.group_pca.pdf', width=6, height = 6)

##make heatmap for candidate genes
#genes <- read.csv('tim_marker_genes_1019.txt', header = FALSE, col.names = 'genes', as.is = TRUE)
#genes <- genes$genes #vectorize
#genes = c('Sh2b3', 'Ly86', 'Csf1r', 'Elane')
#foo <- assay(rld) #save assay as matrix
#assaygenes <- row.names(foo) #extract rownames
#idx <- assaygenes %in% genes #find target genes
#newmat <- assay(rld)[idx,] #subset to target genes
#newmat <- newmat - rowMeans(newmat)
#newdf <- as.data.frame(colData(rld)["group"])
#pheatmap(newmat, annotation_col=newdf, fontsize_row=5, fontsize_col=5)
#dev.copy2pdf(file='mouse_bm_GSE151682.heatmap.pdf', width = 7, height = 10)

##new way... get vst values and get pertinat row and change to gene names manually
heatmap.data <- read.csv('GSE151682_vst_counts_genes2.csv', header=T, row.names=1)
mat <- heatmap.data - rowMeans(heatmap.data)
pheatmap(heatmap.data, fontsize_row=10, fontsize_col=10)
dev.copy2pdf(file='mouse_bm_GSE151682.heatmap.pdf', width = 7, height = 10)
pheatmap(mat, fontsize_row=10, fontsize_col=10)
dev.copy2pdf(file='mouse_bm_GSE151682.heatmap_averaged.pdf', width = 7, height = 10)

##differential expression
dds <- DESeq(dds)
##to get a specific test:
res2 <- results(dds, contrast=c("group", "d3_mid_grade_sepsis_ly6c_neg_gmp", "ly6c_neg_gmp"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22734,]
write.csv(resOrdered2DF, file="mouse_bm_GSE151682.ly6c_neg_gmp.sepsis_vs_ctl.csv")
##to get a specific test:
res2 <- results(dds, contrast=c("group", "d3_mid_grade_sepsis_preneu", "preneu"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22734,]
write.csv(resOrdered2DF, file="mouse_bm_GSE151682.preneu.sepsis_vs_ctl.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("group", "d3_mid_grade_sepsis_proneu1", "proneu1"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22734,]
write.csv(resOrdered2DF, file="mouse_bm_GSE151682.proneu1.sepsis_vs_ctl.csv")

##to get a specific test:
res2 <- results(dds, contrast=c("group", "d3_mid_grade_sepsis_proneu2", "proneu2"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:22734,]
write.csv(resOrdered2DF, file="mouse_bm_GSE151682.proneu2.sepsis_vs_ctl.csv")




##Human blood, from TPM noramalized data
dat <- read.csv('human_blood_GSE107011.TPM.txt', sep = '\t')

##get the value for length, mean, and sd data, doesn't work with nans
##from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
cdata <- ddply(dat, c("tissue", "gene"), summarise,
               N    = length(TPM), mean = mean(TPM),
               sd   = sd(TPM), se   = sd / sqrt(N))
cdata

##graph with all 3 genes in seperate graphs i.e. facetted as bargraph
ggplot(data=cdata, aes(x=tissue, y=mean, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3, width=.2, position=position_dodge(.9)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('human_blood_GSE107011.SH2B.bar_graph.pdf', width = 10, height = 8)
##facetted boxplot
ggplot(dat, aes(x=tissue, y=TPM, fill=gene)) + geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('human_blood_GSE107011.SH2B.boxplot.pdf', width = 10, height = 8)


##Mouse Peripheral Blood -- GSE122597
##read in count and metadata
countData1 <- read.table('GSE122597_Gene_count_table.txt', header=T, row.names=1)
colData1 <- read.table('GSE122597_metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sample_name)
dds
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="GSE122597_norm_counts.csv")

##graph reformatted norm counts for Sh1b genes
dat <- read.csv('GSE122597_for_graphing.txt', sep = '\t')

##get the value for length, mean, and sd data, doesn't work in nans
##from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
cdata <- ddply(dat, c("group", "gene"), summarise,
               N    = length(value), mean = mean(value),
               sd   = sd(value), se   = sd / sqrt(N))
cdata

##graph with all 3 genes in seperate graphs i.e. facetted
ggplot(data=cdata, aes(x=group, y=mean, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3, width=.2, position=position_dodge(.9)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =12)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('mouse_blood.GSE122597.SH2B.bar_graph.pdf', width = 10, height = 8)
##facetted boxplot
ggplot(dat, aes(x=group, y=value, fill=gene)) + geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =12)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('mouse_blood.GSE122597.SH2B.boxplot.pdf', width = 10, height = 8)


##Mouse Peripheral Blood -- GSE109125
##read in count and metadata
countData1 <- read.csv('GSE109125_Gene_count_table_GENCODE_vM25.csv', header=T, row.names=1)
colData1 <- read.table('GSE109125_metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~ sample_name)
dds
##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="GSE122597_norm_counts.csv")

##graph reformatted norm counts for Sh1b genes
dat <- read.csv('GSE122597_for_graphing.txt', sep = '\t')

##get the value for length, mean, and sd data, doesn't work in nans
##from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
cdata <- ddply(dat, c("group", "gene"), summarise,
               N    = length(value), mean = mean(value),
               sd   = sd(value), se   = sd / sqrt(N))
cdata

##graph with all 3 genes in seperate graphs i.e. facetted
ggplot(data=cdata, aes(x=group, y=mean, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3, width=.2, position=position_dodge(.9)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =12)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('mouse_blood.GSE122597.SH2B.bar_graph.pdf', width = 10, height = 8)
##facetted boxplot
ggplot(dat, aes(x=group, y=value, fill=gene)) + geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =12)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('mouse_blood.GSE122597.SH2B.boxplot.pdf', width = 10, height = 8)


