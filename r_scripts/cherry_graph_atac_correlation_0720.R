##load libraries
library("pheatmap")
library("ggplot2")

#workingDir = "/home/atimms/ngs_data/misc/cherry_human_scatac_call_peaks_0720";
#setwd(workingDir);
#workingDir = "/home/atimms/ngs_data/misc/cherry_scatac_call_peaks_0820";
#setwd(workingDir);
##re doing this with wk12 samples 1021
workingDir = "/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_corr_rpts_1021";
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
  ##write files
  corr_csv <- gsub(".txt",".csv",filename[i])
  write.csv(corr_res, file=corr_csv)
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 7, height = 5)
}


##loop through files - 1021 data
filename <- dir('/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_corr_rpts_1021', pattern=glob2rx("correlation*keepdups_summits.txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##graph
  pheatmap(corr_res)
  ##write files
  corr_csv <- gsub(".txt",".csv",filename[i])
  write.csv(corr_res, file=corr_csv)
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 7, height = 5)
}

##make heatmaps that have specific order
new_order <- c('IPSC_c4_bulk', 'IPSC_c5_1_bulk', 'X5wk_bulk', 'X5wk_c5_1_bulk', 'd53_bulk', 'd59_bulk', 'd74_bulk', 'd78_bulk', 'X12wk1_bulk', 'X12wk2_bulk', 'X12wk3_bulk', 'd113_bulk', 'd132_bulk', 'X20wk_bulk', 'X20wk_c5_1_bulk', 'X28.1_bulk', 'X28.2_bulk', 'hu5_bulk', 'hu7_bulk', 'hu8_bulk')
filename <- dir('/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_corr_rpts_1021', pattern=glob2rx("correlation*keepdups_summits.txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##reoder data
  corr_res_ordered <- corr_res[new_order, new_order]
  ##graph
  pheatmap(corr_res_ordered, cluster_cols = F, cluster_rows = F)
  ##write files
  corr_csv <- gsub(".txt","_ordered.csv",filename[i])
  write.csv(corr_res, file=corr_csv)
  corr_pdf <- gsub(".txt","_ordered.pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 7, height = 5)
}



##scatterplots comparing human/organoid data

##one file
corr_data <- read.table('annotated.collapsed_cc.100d.2000size.macs2_q0.000001_summits.average_for_scatterplot.txt', header=T, row.names=1)
head(corr_data)
ggplot(corr_data, aes(x=human, y=organoid)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE) 

##all files
filename <- dir('/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_corr_rpts_1021', pattern=glob2rx("annotated*average_for_scatterplot.txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##graph
  ggplot(corr_data, aes(x=human, y=organoid)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE) 
  ##write files
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 7, height = 5)
}
