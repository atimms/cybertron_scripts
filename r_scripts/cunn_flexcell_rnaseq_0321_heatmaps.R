##load libraries
library("pheatmap")
library("ggplot2")

#workingDir = "/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321";
#workingDir = "/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0621";
workingDir = "/home/atimms/ngs_data/rnaseq/cunn_rnaseq_0220";
setwd(workingDir);


##test different methods...
##norm counts ratio
norm_ratio_data <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.ratios.focal_adhesion.heatmap.txt', header=T, row.names=1, sep='\t')
pheatmap(norm_ratio_data, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='ctl_flna.norm_counts.ratios.focal_adhesion.heatmap.pdf', width = 9, height = 6)
##norm log2fc
norm_logfc_data <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.log2fc.focal_adhesion.heatmap.txt', header=T, row.names=1, sep='\t')
pheatmap(norm_logfc_data, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='ctl_flna.norm_counts.log2fc.focal_adhesion.heatmap.pdf', width = 9, height = 6)
##norm counts
norm_counts_data <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.counts.focal_adhesion.heatmap.txt', header=T, row.names=1, sep='\t')
pheatmap(norm_counts_data, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='ctl_flna.norm_counts.counts.focal_adhesion.heatmap.pdf', width = 9, height = 6)
newmat <- norm_counts_data - rowMeans(norm_counts_data)
pheatmap(newmat, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='ctl_flna.norm_counts.counts_mean.focal_adhesion.heatmap.pdf', width = 9, height = 6)

##norm counts ratio
norm_ratio_data <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.vst_counts.ratios.focal_adhesion.heatmap.txt', header=T, row.names=1, sep='\t')
pheatmap(norm_ratio_data, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='ctl_flna.vst_counts.ratios.focal_adhesion.heatmap.pdf', width = 9, height = 6)
##norm counts
norm_counts_data <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.vst_counts.counts.focal_adhesion.heatmap.txt', header=T, row.names=1, sep='\t')
pheatmap(norm_counts_data, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='ctl_flna.vst_counts.counts.focal_adhesion.heatmap.pdf', width = 9, height = 6)
newmat <- norm_counts_data - rowMeans(norm_counts_data)
pheatmap(newmat, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='ctl_flna.vst_counts.counts_mean.focal_adhesion.heatmap.pdf', width = 9, height = 6)

##loop through file -- log2fc
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321', pattern ="log2fc.*.heatmap.txt")
for(i in 1:length(filename)){
  print(filename[i])
  norm_logfc_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  pheatmap(norm_logfc_data, fontsize_row=6, fontsize_col=8)
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  
}

##loop through file -- ratios
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321', pattern ="ratios.*.heatmap.txt")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  ##not sure pmax is working correctly
  norm_ratio_data = pmax(norm_ratio_data, -10)
  norm_ratio_data = pmin(norm_ratio_data, 10)
  pheatmap(norm_ratio_data, fontsize_row=6, fontsize_col=8)
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  
}

##function for making pngs
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

##loop through file -- ratios another take, using +/- 3 as the max/min
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321', pattern ="ratios.*.heatmap.txt")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data < -3,-3)
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data > 3,3)
  pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  #gene_pdf <- gsub(".txt","_3.pdf",filename[i])
  #dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  a = pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  gene_png <- gsub(".txt","_3.png",filename[i])
  save_pheatmap_png(a, gene_png)
}

##loop through file -- all samples combined, using ratios and +/- 3 as the max/min
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321', pattern ="all_samples.norm_counts.ratios")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data < -3,-3)
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data > 3,3)
  pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  #gene_pdf <- gsub(".txt","_3.pdf",filename[i])
  #dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  a = pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  gene_png <- gsub(".txt","_3.png",filename[i])
  save_pheatmap_png(a, gene_png)
}



##loop through file -- all samples combined, using ratios and +/- 3 as the max/min but actually clustering the samples
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321', pattern ="all_samples.norm_counts.ratios")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data < -3,-3)
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data > 3,3)
  pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8)
  #gene_pdf <- gsub(".txt","_3.pdf",filename[i])
  #dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  a = pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8)
  gene_png <- gsub(".txt","_cluster_samples.png",filename[i])
  save_pheatmap_png(a, gene_png)
}


##loop through file -- all samples combined, use the ratios to calculate correlation
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321', pattern ="all_samples.norm_counts.ratios")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  corr_res <- round(cor(norm_ratio_data, method="spearman"),3)
  pheatmap(corr_res, fontsize_row=8, fontsize_col=8)
  #gene_pdf <- gsub(".txt","_3.pdf",filename[i])
  #dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  a = pheatmap(corr_res, fontsize_row=8, fontsize_col=8)
  gene_png <- gsub(".txt","_correlation.png",filename[i])
  save_pheatmap_png(a, gene_png)
}


##for 0621 data -- flexcell
##loop through file -- all samples combined, using ratios and +/- 3 as the max/min
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0621', pattern ="heatmap.txt")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data < -3,-3)
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data > 3,3)
  pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  #gene_pdf <- gsub(".txt","_3.pdf",filename[i])
  #dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  a = pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  gene_png <- gsub(".txt","_3.png",filename[i])
  save_pheatmap_png(a, gene_png)
}

##for 0621 data -- matrix data
##loop through file -- all samples combined, using ratios and +/- 3 as the max/min
filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_rnaseq_0220', pattern ="heatmap.txt")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data < -3,-3)
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data > 3,3)
  pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  #gene_pdf <- gsub(".txt","_3.pdf",filename[i])
  #dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  a = pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  gene_png <- gsub(".txt","_3.png",filename[i])
  save_pheatmap_png(a, gene_png)
}

##redo some graphs for deema 1021
workingDir = "/home/atimms/ngs_data/rnaseq/deema_heatmaps_1021";
setwd(workingDir);

##function for making pngs
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
##loop through file -- all samples combined, using ratios and +/- 3 as the max/min
filename <- dir('/home/atimms/ngs_data/rnaseq/deema_heatmaps_1021', pattern ="heatmap.txt")
for(i in 1:length(filename)){
  print(filename[i])
  norm_ratio_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data < -3,-3)
  norm_ratio_data = replace(norm_ratio_data, norm_ratio_data > 3,3)
  pheatmap(norm_ratio_data, fontsize_row=4, fontsize_col=8, cluster_cols = F)
  gene_pdf <- gsub(".txt","_3.pdf",filename[i])
  dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
}








##testting stuff
norm_ratio_data <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.ratios.focal_adhesion.heatmap.txt', header=T, row.names=1, sep='\t')
to_remove<-c("ARHGAP6")
norm_ratio_data = norm_ratio_data[!(row.names(norm_ratio_data) %in% to_remove),]
pheatmap(norm_ratio_data, fontsize_row=6, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.ratios.focal_adhesion.heatmap.pdf', width = 9, height = 6)

filename <- dir('/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321', pattern ="ratios.*.heatmap.txt")
for(i in 1:length(filename)){
  print(filename[i])
  norm_logfc_data <- read.csv(filename[i],  header=T, row.names=1, sep='\t')
  pheatmap(norm_logfc_data, fontsize_row=6, fontsize_col=8)
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  #dev.copy2pdf(file=gene_pdf, width = 9, height = 6)
  
}

norm_ratio_data <- read.table('all_samples.norm_counts.ratios.adherens.heatmap.txt', header=T, row.names=1, sep='\t')
##calculate correlation
corr_res <- round(cor(norm_ratio_data, method="spearman"),3)
##graph
pheatmap(corr_res, fontsize = 3)

norm_ratio_data = pmax(norm_ratio_data, -3)
norm_ratio_data = pmin(norm_ratio_data, 3)
head(norm_ratio_data)
norm_ratio_data = replace(norm_ratio_data, norm_ratio_data < -3,-3)
norm_ratio_data = replace(norm_ratio_data, norm_ratio_data > 3,3)
head(norm_ratio_data)


save_pheatmap_png(a, "my_heatmap.png")

a = pheatmap(norm_ratio_data, fontsize_row=6, fontsize_col=8, cluster_cols = F)
png(filename="test.png")
plot(pheatmap(norm_ratio_data, fontsize_row=6, fontsize_col=8, cluster_cols = F))
dev.off()

dev.print(a, 'test.png')

dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.ratios.focal_adhesion.heatmap_3.pdf', width = 9, height = 6)
##get row names
rownames(norm_ratio_data[a$tree_row[["order"]],])
