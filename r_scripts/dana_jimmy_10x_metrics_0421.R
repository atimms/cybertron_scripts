
# source libraries needed for this. the only non-generic one you might need here is reshape2 
# but if something fails because a command can't be found, find that package it's in and install / load it
# for example:
#install.packages("reshape2")
library(reshape2)
library(Seurat)
library("ggplot2")

workingDir = "/archive/bennett_j/Dana10xpractice/andrew_analysis";
setwd(workingDir);

#load("path/to/seurat/object/for/scRNA/or/Visium/object.Rdata")
lm1 <- readRDS(file = "/archive/bennett_j/Dana10xpractice/LM1_harmonysubset.RDS")
lm2 <- readRDS(file = "/archive/bennett_j/Dana10xpractice/LM2_harmonysubset.RDS")
lm38 <- readRDS(file = "/archive/bennett_j/Dana10xpractice/LM38_harmonysubset.RDS")

gene_list = c("PIK3CA","KRAS","DNMT3A", "CALR", "JAK2", "SF3B1") #put whatever genes you want, even just one
#gene_list = c("PIK3CA") #put whatever genes you want, even just one


##lm1
df = as.data.frame(lm1@assays$RNA@counts[gene_list,]) #this is for RNA; the assay will change to SCT instead of RNA for Visium. Also change "SeuratObjectName" to your object name
df$gene = rownames(df)
df = melt(df)

# for plots of genes with very low expression, think about using a y axis log scale or pseudolog scale. example:
# scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 2, 5, 10, 25, 50, 75))

#violin + box plot
ggplot(df, aes(x = reorder(gene, value), y = value, fill = gene, color = gene)) + 
  geom_violin(width = 1)+
  geom_boxplot(width = 0.05, outlier.size = 0.5, outlier.colour = "black", color = "black")+
  labs(x =  "Gene", y = "UMI Count")
ggsave('lm1.violin.pdf', width=6, height = 4) 

#boxplot only
ggplot(df, aes(x = reorder(gene, value), y = value, fill = gene)) + 
  geom_boxplot(width = 0.5, outlier.size = 0.1)+
  labs(x =  "Gene", y = "UMI Count")
ggsave('lm1.boxplot.pdf', width=6, height = 4) 

#histogram
ggplot(df, aes(x = value,  fill = gene))+
  geom_histogram(binwidth = 1)+
  scale_y_continuous(trans ="log10")+
  # facet_wrap(~reorder(gene, value), ncol = 6)+ #only needed for multiple genes if you want to order them differently
  labs(x = "UMI Count", y = "Number of Cells")
ggsave('lm1.histogram.pdf', width=6, height = 4)

#calculate mean of expression for each gene
df_by_gene = split(df, df$gene)
means = unlist(lapply(df_by_gene, function(x)mean(x$value)))
head(means)
means = as.data.frame(means)
means$gene = rownames(means)
means = means[order(means$means),]
means
write.csv(means, file="lm1.mean_expression.csv")

# plot of non-zero expression, aka what % of cells/spots have detectable expression of the gene
temp = as.data.frame(lm1@assays$RNA@counts[gene_list,]) > 0
nonzero = data.frame(gene = rownames(temp), percent_nonzero = rowSums(temp)/ncol(temp)*100)

ggplot(nonzero, aes(x = reorder(gene, percent_nonzero), y = percent_nonzero, fill = gene, color = gene))+
  geom_bar(stat = "identity", color = "black")+
  labs( x = "Gene", y = "% of Cells with >0 UMIs")
ggsave('lm1.non_zero_expression.pdf', width=6, height = 4)

##lm2
df = as.data.frame(lm2@assays$RNA@counts[gene_list,]) #this is for RNA; the assay will change to SCT instead of RNA for Visium. Also change "SeuratObjectName" to your object name
df$gene = rownames(df)
df = melt(df)

# for plots of genes with very low expression, think about using a y axis log scale or pseudolog scale. example:
# scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 2, 5, 10, 25, 50, 75))

#violin + box plot
ggplot(df, aes(x = reorder(gene, value), y = value, fill = gene, color = gene)) + 
  geom_violin(width = 1)+
  geom_boxplot(width = 0.05, outlier.size = 0.5, outlier.colour = "black", color = "black")+
  labs(x =  "Gene", y = "UMI Count")
ggsave('lm2.violin.pdf', width=6, height = 4) 

#boxplot only
ggplot(df, aes(x = reorder(gene, value), y = value, fill = gene)) + 
  geom_boxplot(width = 0.5, outlier.size = 0.1)+
  labs(x =  "Gene", y = "UMI Count")
ggsave('lm2.boxplot.pdf', width=6, height = 4) 

#histogram
ggplot(df, aes(x = value,  fill = gene))+
  geom_histogram(binwidth = 1)+
  scale_y_continuous(trans ="log10")+
  # facet_wrap(~reorder(gene, value), ncol = 6)+ #only needed for multiple genes if you want to order them differently
  labs(x = "UMI Count", y = "Number of Cells")
ggsave('lm2.histogram.pdf', width=6, height = 4)

#calculate mean of expression for each gene
df_by_gene = split(df, df$gene)
means = unlist(lapply(df_by_gene, function(x)mean(x$value)))
head(means)
means = as.data.frame(means)
means$gene = rownames(means)
means = means[order(means$means),]
means
write.csv(means, file="lm2.mean_expression.csv")

# plot of non-zero expression, aka what % of cells/spots have detectable expression of the gene
temp = as.data.frame(lm2@assays$RNA@counts[gene_list,]) > 0
nonzero = data.frame(gene = rownames(temp), percent_nonzero = rowSums(temp)/ncol(temp)*100)

ggplot(nonzero, aes(x = reorder(gene, percent_nonzero), y = percent_nonzero, fill = gene, color = gene))+
  geom_bar(stat = "identity", color = "black")+
  labs( x = "Gene", y = "% of Cells with >0 UMIs")
ggsave('lm2.non_zero_expression.pdf', width=6, height = 4)

##lm38
df = as.data.frame(lm38@assays$RNA@counts[gene_list,]) #this is for RNA; the assay will change to SCT instead of RNA for Visium. Also change "SeuratObjectName" to your object name
df$gene = rownames(df)
df = melt(df)

# for plots of genes with very low expression, think about using a y axis log scale or pseudolog scale. example:
# scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 2, 5, 10, 25, 50, 75))

#violin + box plot
ggplot(df, aes(x = reorder(gene, value), y = value, fill = gene, color = gene)) + 
  geom_violin(width = 1)+
  geom_boxplot(width = 0.05, outlier.size = 0.5, outlier.colour = "black", color = "black")+
  labs(x =  "Gene", y = "UMI Count")
ggsave('lm38.violin.pdf', width=6, height = 4) 

#boxplot only
ggplot(df, aes(x = reorder(gene, value), y = value, fill = gene)) + 
  geom_boxplot(width = 0.5, outlier.size = 0.1)+
  labs(x =  "Gene", y = "UMI Count")
ggsave('lm38.boxplot.pdf', width=6, height = 4) 

#histogram
ggplot(df, aes(x = value,  fill = gene))+
  geom_histogram(binwidth = 1)+
  scale_y_continuous(trans ="log10")+
  # facet_wrap(~reorder(gene, value), ncol = 6)+ #only needed for multiple genes if you want to order them differently
  labs(x = "UMI Count", y = "Number of Cells")
ggsave('lm38.histogram.pdf', width=6, height = 4)

#calculate mean of expression for each gene
df_by_gene = split(df, df$gene)
means = unlist(lapply(df_by_gene, function(x)mean(x$value)))
head(means)
means = as.data.frame(means)
means$gene = rownames(means)
means = means[order(means$means),]
means
write.csv(means, file="lm38.mean_expression.csv")

# plot of non-zero expression, aka what % of cells/spots have detectable expression of the gene
temp = as.data.frame(lm38@assays$RNA@counts[gene_list,]) > 0
nonzero = data.frame(gene = rownames(temp), percent_nonzero = rowSums(temp)/ncol(temp)*100)

ggplot(nonzero, aes(x = reorder(gene, percent_nonzero), y = percent_nonzero, fill = gene, color = gene))+
  geom_bar(stat = "identity", color = "black")+
  labs( x = "Gene", y = "% of Cells with >0 UMIs")
ggsave('lm38.non_zero_expression.pdf', width=6, height = 4)



