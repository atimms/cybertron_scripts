##load libraries
library("pheatmap")


workingDir = "/archive/mirzaa_g/exomes/result_files_0321/graphing/wei_pathways_1221";
setwd(workingDir);

##heatmap from dxgroup/pathway

##read in data
phenotype_data <- read.table('dx1_pathway1.txt', header=T, row.names=1, sep='\t')

##graph and make pdf
pheatmap(phenotype_data, cluster_rows = F, cluster_cols = F)
dev.copy2pdf(file='dx1_pathway1.raw_counts.heatmap.pdf', height = 5, width = 10 )

##by percentages
phenotype_data_percentages_dx = phenotype_data / rowSums(phenotype_data)[row(phenotype_data)] * 100
pheatmap(phenotype_data_percentages_dx, cluster_rows = F, cluster_cols = F)
dev.copy2pdf(file='dx1_pathway1.percentage_dx1.heatmap.pdf', height = 5, width = 10 )
phenotype_data_percentages_pathway = phenotype_data / colSums(phenotype_data)[col(phenotype_data)] * 100
pheatmap(phenotype_data_percentages_pathway, cluster_rows = F, cluster_cols = F)
dev.copy2pdf(file='dx1_pathway1.percentage_pathway1.heatmap.pdf', height = 5, width = 10 )
