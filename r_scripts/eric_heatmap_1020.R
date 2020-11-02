##load libraries
library("pheatmap")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/misc/eric_heatmap_1020";
setwd(workingDir);


##read in file
hm_data <- read.table('heatmap_data1.txt', header=T, row.names=1, sep='\t')
head(hm_data)

##scale on columns
pheatmap(hm_data, scale="row", cluster_cols = F)
dev.copy2pdf(file='eric_hm_scaled_1020.pdf', height = 5 )

pheatmap(hm_data, scale="row", cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50))
dev.copy2pdf(file='eric_hm_scaled_br_1020.pdf', height = 5 )


##get %s accross each row
hm_perc <- hm_data/rowSums(hm_data) * 100
##make heatmap
pheatmap(hm_perc, cluster_cols = F)
dev.copy2pdf(file='eric_hm1_1020.pdf', height = 5)


