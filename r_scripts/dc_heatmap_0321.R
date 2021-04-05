library("ggplot2")

workingDir = "/home/atimms/ngs_data/misc/dc_heatmap";
setwd(workingDir);

Data1 <- read.table('qpcr_data_nogapdh.txt', header=T, row.names=1)
mat = Data1
mat <- mat - rowMeans(mat)
pheatmap(log10(Data1), fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='qpcr_nogapdh.log10.pdf', width = 10, height = 7)
pheatmap(Data1, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='qpcr_nogapdh.std.pdf', width = 10, height = 7)

pheatmap(mat, fontsize_row=8, fontsize_col=8)
Data2 <- read.table('qpcr_data.txt', header=T, row.names=1)
pheatmap(log10(Data2), fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='qpcr.log10.pdf', width = 10, height = 7)
pheatmap(Data2, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='qpcr.std.pdf', width = 10, height = 7)
