##load libraries
library("pheatmap")
library("ggplot2")
library(ggbiplot)
library(data.table)

workingDir = "/home/atimms/ngs_data/misc/cunn_geomx_0421";
setwd(workingDir);


##take1

##using slide1_calrite

##heatmap/clustering
##sample/location info
s1c_sample_info = read.table('slide1_calrite_samples.txt', header=T, row.names=1, sep='\t')
##all the genes
s1c_all_data <- read.table('slide1_calrite_all.txt', header=T, row.names=1, sep='\t')
s1c_all_data_log <- log(s1c_all_data +1, 2)
pheatmap(s1c_all_data_log, annotation = s1c_sample_info, show_rownames=F)
dev.copy2pdf(file='slide1_calrite_all.heatmap.pdf', width = 8, height = 6)
##only genes with a value >30
s1c_30_data <- read.table('slide1_calrite_30.txt', header=T, row.names=1, sep='\t')
s1c_30_data_log <- log(s1c_all_data +1, 2)
pheatmap(s1c_30_data_log, annotation = s1c_sample_info, show_rownames=F)
dev.copy2pdf(file='slide1_calrite_30.heatmap.pdf', width = 8, height = 6)

##pca
#https://www.datacamp.com/community/tutorials/pca-analysis-r
s1c_all_data <- read.table('slide1_calrite_all.txt', header=T, row.names=1, sep='\t')
t_s1c_all_data <- transpose(s1c_all_data)
colnames(t_s1c_all_data) <- rownames(s1c_all_data)
rownames(t_s1c_all_data) <- colnames(s1c_all_data)
s1c_all.pca <- prcomp(t_s1c_all_data, center = TRUE,scale. = TRUE)
s1c_all.pca <- prcomp(t_s1c_all_data)
summary(s1c_all.pca)
s1c_sample_info = read.table('slide1_calrite_samples.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1c_all.pca, labels=rownames(t_s1c_all_data),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_all.pca.pdf', width = 8, height = 6)
s1c_all.pca <- prcomp(t_s1c_all_data)
ggbiplot(s1c_all.pca, labels=rownames(t_s1c_all_data),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_all.pca_no_scale.pdf', width = 8, height = 6)

s1c_30_data <- read.table('slide1_calrite_30.txt', header=T, row.names=1, sep='\t')
t_s1c_30_data <- transpose(s1c_30_data)
colnames(t_s1c_30_data) <- rownames(s1c_30_data)
rownames(t_s1c_30_data) <- colnames(s1c_30_data)
s1c_30.pca <- prcomp(t_s1c_30_data, center = TRUE,scale. = TRUE)
summary(s1c_30.pca)
s1c_sample_info = read.table('slide1_calrite_samples.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1c_all.pca, labels=rownames(t_s1c_30_data),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_30.pca.pdf', width = 8, height = 6)
s1c_30.pca <- prcomp(t_s1c_30_data)
ggbiplot(s1c_all.pca, labels=rownames(t_s1c_30_data),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_30.pca_no_scale.pdf', width = 8, height = 6)

##using slide1_edta
##heatmap/clustering
##sample/location info
s1e_sample_info = read.table('slide1_edta_samples.txt', header=T, row.names=1, sep='\t')
##all the genes
s1e_all_data <- read.table('slide1_edta_all.txt', header=T, row.names=1, sep='\t')
s1e_all_data_log <- log(s1e_all_data +1, 2)
pheatmap(s1e_all_data_log, annotation = s1e_sample_info, show_rownames=F)
dev.copy2pdf(file='slide1_edta_all.heatmap.pdf', width = 8, height = 6)

##pca
s1e_all_data <- read.table('slide1_edta_all.txt', header=T, row.names=1, sep='\t')
t_s1e_all_data <- transpose(s1e_all_data)
colnames(t_s1e_all_data) <- rownames(s1e_all_data)
rownames(t_s1e_all_data) <- colnames(s1e_all_data)
s1e_all.pca <- prcomp(t_s1e_all_data, center = TRUE,scale. = TRUE)
summary(s1e_all.pca)
s1e_sample_info = read.table('slide1_edta_samples.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1e_all.pca, labels=rownames(t_s1e_all_data),var.axes=FALSE, groups=s1e_sample_info$location)
dev.copy2pdf(file='slide1_edta_all.pca.pdf', width = 8, height = 6)

##remove sample 36

##heatmap/clustering
##sample/location info
s1e_sample_info = read.table('slide1_edta_samples_m36.txt', header=T, row.names=1, sep='\t')
##all the genes
s1e_all_data <- read.table('slide1_edta_m36.txt', header=T, row.names=1, sep='\t')
s1e_all_data_log <- log(s1e_all_data +1, 2)
pheatmap(s1e_all_data_log, annotation = s1e_sample_info, show_rownames=F)
dev.copy2pdf(file='slide1_edta_minus36.heatmap.pdf', width = 8, height = 6)

##pca
s1e_all_data <- read.table('slide1_edta_m36.txt', header=T, row.names=1, sep='\t')
t_s1e_all_data <- transpose(s1e_all_data)
colnames(t_s1e_all_data) <- rownames(s1e_all_data)
rownames(t_s1e_all_data) <- colnames(s1e_all_data)
s1e_all.pca <- prcomp(t_s1e_all_data, center = TRUE,scale. = TRUE)
summary(s1e_all.pca)
s1e_sample_info = read.table('slide1_edta_samples_m36.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1e_all.pca, labels=rownames(t_s1e_all_data),var.axes=FALSE, groups=s1e_sample_info$location)
dev.copy2pdf(file='slide1_edta_minus36.pca.pdf', width = 8, height = 6)
s1e_all.pca <- prcomp(t_s1e_all_data)
ggbiplot(s1e_all.pca, labels=rownames(t_s1e_all_data),var.axes=FALSE, groups=s1e_sample_info$location)
dev.copy2pdf(file='slide1_edta_minus36.pca_noscale.pdf', width = 8, height = 6)




##take2 0514 using slide1_calrite

##heatmap/clustering
##sample/location info
s1c_sample_info = read.table('slide1_calrite_samples_0514.txt', header=T, row.names=1, sep='\t')
##all the genes
s1c_all_data <- read.table('slide1_calrite_0514.txt', header=T, row.names=1, sep='\t')
s1c_all_data_log <- log(s1c_all_data +1, 2)
pheatmap(s1c_all_data_log, annotation = s1c_sample_info, show_rownames=F, border_color = NA)
dev.copy2pdf(file='slide1_calrite_all.heatmap_0514.pdf', width = 8, height = 6)
##genes max >30
s1c_all_data <- read.table('slide1_calrite30_0514.txt', header=T, row.names=1, sep='\t')
s1c_all_data_log <- log(s1c_all_data +1, 2)
pheatmap(s1c_all_data_log, annotation = s1c_sample_info, show_rownames=F, border_color = NA)
dev.copy2pdf(file='slide1_calrite_30.heatmap_0514.pdf', width = 8, height = 6)
##on individual gene sets
s1c_sample_info = read.table('slide1_calrite_samples_0514.txt', header=T, row.names=1, sep='\t')
s1c_gs_data <- read.table('slide1_calrite_PI3K_0514.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
pheatmap(s1c_gs_data_log, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_pi3k.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_pi3k.heatmap_mean_0514.pdf', width = 8, height = 6)
#pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4,color=rainbow(10))
#dev.copy2pdf(file='slide1_calrite_pi3k.heatmap_mean2_0514.pdf', width = 8, height = 6)
#pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4,color=colorRampPalette(c("green", "white", "red"))(10))
#dev.copy2pdf(file='slide1_calrite_pi3k.heatmap_mean3_0514.pdf', width = 8, height = 6)
s1c_gs_data <- read.table('slide1_calrite_ecm_0514.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
pheatmap(s1c_gs_data_log, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_ecm.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_ecm.heatmap_mean_0514.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_focal_adhesion_0514.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
pheatmap(s1c_gs_data_log, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_focal_adhesion.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_focal_adhesion.heatmap_mean_0514.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_pluri_0514.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
pheatmap(s1c_gs_data_log, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_pluri.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_pluri.heatmap_mean_0514.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_tgfb_0514.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
pheatmap(s1c_gs_data_log, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_tgfb.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_tgfb.heatmap_mean_0514.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_wnt_0514.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
pheatmap(s1c_gs_data_log, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_wnt.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_wnt.heatmap_mean_0514.pdf', width = 8, height = 6)

##pca
#https://www.datacamp.com/community/tutorials/pca-analysis-r
s1c_all_data <- read.table('slide1_calrite_0514.txt', header=T, row.names=1, sep='\t')
t_s1c_all_data <- transpose(s1c_all_data)
colnames(t_s1c_all_data) <- rownames(s1c_all_data)
rownames(t_s1c_all_data) <- colnames(s1c_all_data)
s1c_all.pca <- prcomp(t_s1c_all_data, center = TRUE,scale. = TRUE)
#s1c_all.pca <- prcomp(t_s1c_all_data)
summary(s1c_all.pca)
s1c_sample_info = read.table('slide1_calrite_samples_0514.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1c_all.pca, labels=rownames(t_s1c_all_data),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_all.pca12.0514.pdf', width = 8, height = 6)
ggbiplot(s1c_all.pca, labels=rownames(t_s1c_all_data), choices=c(3,4),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_all.pca34.0514.pdf', width = 8, height = 6)
##>30 data
s1c_30_data <- read.table('slide1_calrite30_0514.txt', header=T, row.names=1, sep='\t')
t_s1c_30_data <- transpose(s1c_30_data)
colnames(t_s1c_30_data) <- rownames(s1c_30_data)
rownames(t_s1c_30_data) <- colnames(s1c_30_data)
s1c_30.pca <- prcomp(t_s1c_30_data, center = TRUE,scale. = TRUE)
summary(s1c_30.pca)
s1c_sample_info = read.table('slide1_calrite_samples_0514.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1c_30.pca, labels=rownames(t_s1c_30_data),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_30.pca12.0514.pdf', width = 8, height = 6)
ggbiplot(s1c_30.pca, labels=rownames(t_s1c_30_data), choices=c(3,4),var.axes=FALSE, groups=s1c_sample_info$location)
dev.copy2pdf(file='slide1_calrite_30.pca34.0514.pdf', width = 8, height = 6)


##take2 0514 using slide1_edta

##heatmap/clustering
##sample/location info
s1e_sample_info = read.table('slide1_edta_samples_0514.txt', header=T, row.names=1, sep='\t')
##all the genes
s1e_all_data <- read.table('slide1_edta_0514.txt', header=T, row.names=1, sep='\t')
s1e_all_data_log <- log(s1e_all_data +1, 2)
pheatmap(s1e_all_data_log, annotation = s1e_sample_info, show_rownames=F, border_color = NA)
dev.copy2pdf(file='slide1_edta_all.heatmap_0514.pdf', width = 8, height = 6)
##genes max >30
s1e_all_data <- read.table('slide1_edta30_0514.txt', header=T, row.names=1, sep='\t')
s1e_all_data_log <- log(s1e_all_data +1, 2)
pheatmap(s1e_all_data_log, annotation = s1e_sample_info, show_rownames=F, border_color = NA)
dev.copy2pdf(file='slide1_edta_30.heatmap_0514.pdf', width = 8, height = 6)

##pca
#https://www.datacamp.com/community/tutorials/pca-analysis-r
s1e_all_data <- read.table('slide1_edta_0514.txt', header=T, row.names=1, sep='\t')
t_s1e_all_data <- transpose(s1e_all_data)
colnames(t_s1e_all_data) <- rownames(s1e_all_data)
rownames(t_s1e_all_data) <- colnames(s1e_all_data)
s1e_all.pca <- prcomp(t_s1e_all_data, center = TRUE,scale. = TRUE)
#s1e_all.pca <- prcomp(t_s1e_all_data)
summary(s1e_all.pca)
s1e_sample_info = read.table('slide1_edta_samples_0514.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1e_all.pca, labels=rownames(t_s1e_all_data),var.axes=FALSE, groups=s1e_sample_info$location)
dev.copy2pdf(file='slide1_edta_all.pca12.0514.pdf', width = 8, height = 6)
ggbiplot(s1e_all.pca, labels=rownames(t_s1e_all_data), choices=c(3,4),var.axes=FALSE, groups=s1e_sample_info$location)
dev.copy2pdf(file='slide1_edta_all.pca34.0514.pdf', width = 8, height = 6)
##>30 data
s1e_30_data <- read.table('slide1_edta30_0514.txt', header=T, row.names=1, sep='\t')
t_s1e_30_data <- transpose(s1e_30_data)
colnames(t_s1e_30_data) <- rownames(s1e_30_data)
rownames(t_s1e_30_data) <- colnames(s1e_30_data)
s1e_30.pca <- prcomp(t_s1e_30_data, center = TRUE,scale. = TRUE)
summary(s1e_30.pca)
s1e_sample_info = read.table('slide1_edta_samples_0514.txt', header=T, row.names=1, sep='\t')
ggbiplot(s1e_30.pca, labels=rownames(t_s1e_30_data),var.axes=FALSE, groups=s1e_sample_info$location)
dev.copy2pdf(file='slide1_edta_30.pca12.0514.pdf', width = 8, height = 6)
ggbiplot(s1e_30.pca, labels=rownames(t_s1e_30_data), choices=c(3,4),var.axes=FALSE, groups=s1e_sample_info$location)
dev.copy2pdf(file='slide1_edta_30.pca34.0514.pdf', width = 8, height = 6)

##on individual gene sets
s1e_sample_info = read.table('slide1_edta_samples_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data <- read.table('slide1_edta_PI3K_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
pheatmap(s1e_gs_data_log, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_pi3k.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_pi3k.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_ecm_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
pheatmap(s1e_gs_data_log, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_ecm.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_ecm.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_focal_adhesion_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
pheatmap(s1e_gs_data_log, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_focal_adhesion.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_focal_adhesion.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_pluri_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
pheatmap(s1e_gs_data_log, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_pluri.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_pluri.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_tgfb_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
pheatmap(s1e_gs_data_log, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_tgfb.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_tgfb.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_wnt_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
pheatmap(s1e_gs_data_log, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_wnt.heatmap_0514.pdf', width = 8, height = 6)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_wnt.heatmap_mean_0514.pdf', width = 8, height = 6)




##repeat individual gene sets on edta without sample 36
s1e_sample_info = read.table('slide1_edta_m36_samples_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data <- read.table('slide1_edta_m36_PI3K_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_m36_pi3k.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_m36_ecm_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_m36_ecm.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_m36_focal_adhesion_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_m36_focal_adhesion.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_m36_pluri_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_m36_pluri.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_m36_tgfb_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_m36_tgfb.heatmap_mean_0514.pdf', width = 8, height = 6)

s1e_gs_data <- read.table('slide1_edta_m36_wnt_0514.txt', header=T, row.names=1, sep='\t')
s1e_gs_data_log <- log(s1e_gs_data +1, 2)
mat <- s1e_gs_data_log - rowMeans(s1e_gs_data_log)
pheatmap(mat, annotation = s1e_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_edta_m36_wnt.heatmap_mean_0514.pdf', width = 8, height = 6)



##take3 0515 using slide1_calrite but a subset of sample

##heatmap/clustering
##sample/location info
s1c_sample_info = read.table('slide1_calrite_samples_0515.txt', header=T, row.names=1, sep='\t')
##all the genes
s1c_all_data <- read.table('slide1_calrite_0515.txt', header=T, row.names=1, sep='\t')
s1c_all_data_log <- log(s1c_all_data +1, 2)
mat <- s1c_all_data_log - rowMeans(s1c_all_data_log)
pheatmap(mat, annotation = s1c_sample_info, show_rownames=F, border_color = NA)
dev.copy2pdf(file='slide1_calrite_all.heatmap_0515.pdf', width = 8, height = 6)

##genesets
s1c_gs_data <- read.table('slide1_calrite_wnt_0515.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_wnt.heatmap_mean_0515.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_ecm_0515.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_ecm.heatmap_mean_0515.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_focal_adhesion_0515.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_focal_adhesion.heatmap_mean_0515.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_PI3K_0515.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_PI3K.heatmap_mean_0515.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_pluri_0515.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_pluri.heatmap_mean_0515.pdf', width = 8, height = 6)

s1c_gs_data <- read.table('slide1_calrite_tgfb_0515.txt', header=T, row.names=1, sep='\t')
s1c_gs_data_log <- log(s1c_gs_data +1, 2)
mat <- s1c_gs_data_log - rowMeans(s1c_gs_data_log)
pheatmap(mat, annotation = s1c_sample_info, fontsize_row=4, border_color = NA)
dev.copy2pdf(file='slide1_calrite_tgfb.heatmap_mean_0515.pdf', width = 8, height = 6)
