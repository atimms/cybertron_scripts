##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)

##move to working directory
setwd('/home/atimms/ngs_data/exomes/working/ghayda_laser_exomes_0421')

##single file
k10_r5_hgdp <- read.csv('LP97-105.k10_r5.hgdp_all.txt', sep ='\t')
ggplot(aes(x = pc1, y = pc2, color = continent), data = k10_r5_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))

##loop through file -- k6_r5
filename <- dir('/home/atimms/ngs_data/exomes/working/ghayda_laser_exomes_0421', pattern =".k6_r5.hgdp_all.txt")
for(i in 1:length(filename)){
  k6_r5_hgdp <- read.csv(filename[i], sep ='\t')
  ggplot(aes(x = pc1, y = pc2, color = continent), data = k6_r5_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
  gene_pdf <- gsub(".txt",".pc1_pc2.pdf",filename[i])
  ggsave(gene_pdf, height = 4, width = 12)
  ggplot(aes(x = pc3, y = pc4, color = continent), data = k6_r5_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
  gene_pdf <- gsub(".txt",".pc3_pc4.pdf",filename[i])
  ggsave(gene_pdf, height = 4, width = 12)
}

##loop through file -- k10_r10
filename <- dir('/home/atimms/ngs_data/exomes/working/ghayda_laser_exomes_0421', pattern =".k10_r10.hgdp_all.txt")
for(i in 1:length(filename)){
  k10_r10_hgdp <- read.csv(filename[i], sep ='\t')
  ggplot(aes(x = pc1, y = pc2, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
  gene_pdf <- gsub(".txt",".pc1_pc2.pdf",filename[i])
  ggsave(gene_pdf, height = 4, width = 12)
  ggplot(aes(x = pc3, y = pc4, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
  gene_pdf <- gsub(".txt",".pc3_pc4.pdf",filename[i])
  ggsave(gene_pdf, height = 4, width = 12)
  ggplot(aes(x = pc5, y = pc6, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
  gene_pdf <- gsub(".txt",".pc5_pc6.pdf",filename[i])
  ggsave(gene_pdf, height = 4, width = 12)
  ggplot(aes(x = pc7, y = pc8, color = continent), data = k10_r10_hgdp) + geom_point() + scale_color_manual("Population", values=c("orange", "purple", "red", "pink", "blue", "yellow", "green", "#000000"))
  gene_pdf <- gsub(".txt",".pc7_pc8.pdf",filename[i])
  ggsave(gene_pdf, height = 4, width = 12)
}
