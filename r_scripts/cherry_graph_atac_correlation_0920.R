##load libraries
library("pheatmap")
library("ggplot2")

##test on one file
corr_data <- read.table('correlation.all_samples_cell_classes.combined_beds.100d.2000size.macs2_q0.01_summits.txt', header=T, row.names=1)
head(corr_data)
##get corr to 3 decimal points
corr_res <- round(cor(corr_data, method="spearman"),3)
head(corr_res)
pheatmap(corr_res, fontsize = 3)
dev.copy2pdf(file='test.pdf', width = 15, height = 10)


##data set1 - combined_peaks a bunch of different ways
workingDir = "/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/all_samples_cell_classes_092920/combined_peaks";
setwd(workingDir);
##loop through files
filename <- dir('/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/all_samples_cell_classes_092920/combined_peaks', pattern=glob2rx("correlation*txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##graph
  pheatmap(corr_res, fontsize = 3)
  ##write files
  corr_csv <- gsub(".txt",".csv",filename[i])
  write.csv(corr_res, file=corr_csv)
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 15, height = 10)
}

##data set2 - individual_peaks a bunch of different ways
workingDir = "/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/all_samples_cell_classes_092920/individual_peaks";
setwd(workingDir);
##loop through files
filename <- dir('/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/all_samples_cell_classes_092920/individual_peaks', pattern=glob2rx("correlation*txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##graph
  pheatmap(corr_res, fontsize = 3)
  ##write files
  corr_csv <- gsub(".txt",".csv",filename[i])
  write.csv(corr_res, file=corr_csv)
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 15, height = 10)
}


##data set3 - data collapsed by cell class
workingDir = "/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/collapsed_cell_classes_093020";
setwd(workingDir);
##loop through files
filename <- dir('/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/collapsed_cell_classes_093020', pattern=glob2rx("correlation*txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##graph
  pheatmap(corr_res, fontsize = 3)
  ##write files
  corr_csv <- gsub(".txt",".csv",filename[i])
  write.csv(corr_res, file=corr_csv)
  corr_pdf <- gsub(".txt",".pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 15, height = 10)
}


##data set4 - data collapsed by cell class properly 10/20
workingDir = "/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/collapsed_cell_classes_100520";
setwd(workingDir);
##loop through files
filename <- dir('/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_calls_0920/collapsed_cell_classes_100520', pattern=glob2rx("correlation*txt"))
for(i in 1:length(filename)){
  ##read in data
  corr_data <- read.table(filename[i], header=T, row.names=1)
  ##calculate correlation
  corr_res <- round(cor(corr_data, method="spearman"),3)
  ##graph
  pheatmap(corr_res, fontsize = 7)
  ##write files
  corr_csv <- gsub(".txt",".csv",filename[i])
  write.csv(corr_res, file=corr_csv)
  corr_pdf <- gsub(".txt",".counts.pdf",filename[i])
  dev.copy2pdf(file=corr_pdf, width = 15, height = 10)
  ##difference from mean accross the rows
  mat <- corr_data - rowMeans(corr_data)
  corr_res2 <- round(cor(mat, method="spearman"),3)
  pheatmap(corr_res2, fontsize = 7)
  corr_pdf2 <- gsub(".txt",".dev_from_mean.pdf",filename[i])
  dev.copy2pdf(file=corr_pdf2, width = 15, height = 10)
  ##calculate % accross rows
  mat2 <- corr_data/rowSums(corr_data) * 100
  mat3 = na.omit(mat2)
  corr_res3 <- round(cor(mat3, method="spearman"),3)
  pheatmap(corr_res3, fontsize = 7)
  corr_pdf2 <- gsub(".txt",".percetage.pdf",filename[i])
  dev.copy2pdf(file=corr_pdf2, width = 15, height = 10) 
}

corr_data <- read.table('correlation.collapsed_cc_only.100d.2000size.macs2_q0.01_summits.txt', header=T, row.names=1)
mat <- corr_data - rowMeans(corr_data)
mat2 <- corr_data - rowMeans(corr_data)
mat2 <- corr_data/rowSums(corr_data) * 100
write.csv(mat, file="from_mean_temp.csv")
write.csv(mat2, file="perc_temp.csv")
corr_res3 <- round(cor(mat, method="spearman"),3)

mat3 = na.omit(mat2)
