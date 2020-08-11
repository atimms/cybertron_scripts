##load libraries
library("pheatmap")
library("ggplot2")

#workingDir = "/home/atimms/ngs_data/misc/cherry_human_scatac_call_peaks_0720";
#setwd(workingDir);
workingDir = "/home/atimms/ngs_data/misc/cherry_scatac_call_peaks_0820";
setwd(workingDir);

##test on one file
corr_data <- read.table('correlation.100d.macs2.bampe_p1e-10_keepdups_summits.txt', header=T, row.names=1)
head(corr_data)
##get corr to 3 decimal points
corr_res <- round(cor(corr_data, method="spearman"),3)
head(corr_res)
pheatmap(corr_res)
dev.copy2pdf(file='correlation.100d.macs2.bampe_p1e-10_keepdups_summits.pdf', width = 7, height = 5)


##loop through files -0720 data
filename <- dir('/home/atimms/ngs_data/misc/cherry_human_scatac_call_peaks_0720', pattern=glob2rx("correlation*txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##graph
  pheatmap(corr_res)
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 7, height = 5)
}

##loop through files - 0820 data
filename <- dir('/home/atimms/ngs_data/misc/cherry_scatac_call_peaks_0820', pattern=glob2rx("correlation*txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##graph
  pheatmap(corr_res)
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 7, height = 5)
}
