##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)
library(dplyr)
library(scales)

setwd('/home/atimms/ngs_data/enu_mapping/dave_zf_sex_1220')

##testing on one
hom_data <- read.csv('dave_zf_sex_1220.mom_het_dad_wt.q50_cov10.no_alt.danRer11_10000kb_1000kb.combined_for_r.txt', sep = '\t')
goodChrOrder <- paste("chr",c(1:25),sep="")
hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
hom_data
#ggplot(data=hom_data,aes(x=start, y=average_aaf, fill=test),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid( ~ chr, scales="free", space="free_x") + 
#  theme_bw() + theme(axis.text.x = element_blank())
ggplot(data=hom_data,aes(x=start, y=value, group=test, color=test)) +
  geom_line() + facet_grid(test ~ chr, scales="free", space="free_x") + theme(axis.text.x = element_blank())

hom_data = read.csv('dave_zf_sex_1220.mom_het_dad_wt.q50_cov10.no_alt.danRer11_10000kb_1000kb.aaf_for_all.txt', sep = '\t')
goodChrOrder <- paste("chr",c(1:25),sep="")
hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
ggplot(data=hom_data,aes(x=start, y=average_aaf, group=test, color=test)) +
  geom_line() + facet_grid( ~ chr, scales="free", space="free_x") + theme(axis.text.x = element_blank()) +
  scale_color_manual(values=c("pink", "red", "skyblue", 'navyblue'))

##loop through files with combined aaf values 011921
filename <- dir('/home/atimms/ngs_data/enu_mapping/dave_zf_sex_1220', pattern ="aaf_for_all.txt")
for(i in 1:length(filename)){
  ##read in hom data
  hom_data <- read.csv(filename[i], sep ='\t')
  ##reorder chromosome, if using 'chr'
  goodChrOrder <- paste("chr",c(1:25),sep="")
  hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  ##line graph
  ggplot(data=hom_data,aes(x=start, y=average_aaf, group=test, color=test)) +
    geom_line() + facet_grid( ~ chr, scales="free", space="free_x") + theme(axis.text.x = element_blank()) +
    scale_color_manual(values=c("pink", "red", "skyblue", 'navyblue'))
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  ggsave(gene_pdf, width = 20, height = 5)
}


##loop through files with aaf values
filename <- dir('/home/atimms/ngs_data/enu_mapping/dave_zf_sex_1220', pattern ="aaf_for_r.txt")
for(i in 1:length(filename)){
  ##read in hom data
  hom_data <- read.csv(filename[i], sep ='\t')
  ##reorder chromosome, if using 'chr'
  goodChrOrder <- paste("chr",c(1:25),sep="")
  hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  ##line graph
  ggplot(data=hom_data,aes(x=start, y=average_aaf, group=test, color=test)) +
    geom_line() + facet_grid( ~ chr, scales="free", space="free_x") + theme(axis.text.x = element_blank())
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  ggsave(gene_pdf, width = 20, height = 5)
}

##loop through files with all the values
filename <- dir('/home/atimms/ngs_data/enu_mapping/dave_zf_sex_1220', pattern ="combined_for_r.txt")
for(i in 1:length(filename)){
  ##read in hom data
  hom_data <- read.csv(filename[i], sep ='\t')
  ##reorder chromosome, if using 'chr'
  goodChrOrder <- paste("chr",c(1:25),sep="")
  hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  ##line graph
  ggplot(data=hom_data,aes(x=start, y=value, group=test, color=test)) +
    geom_line() + facet_grid(test ~ chr, scales="free", space="free_x") + theme(axis.text.x = element_blank())
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  ggsave(gene_pdf, width = 30, height = 10)
}
