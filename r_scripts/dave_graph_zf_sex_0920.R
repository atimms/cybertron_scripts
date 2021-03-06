##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)
library(dplyr)
library(scales)

setwd('/home/atimms/ngs_data/enu_mapping/dave_zf_sex_0920')

##testing on one
hom_data <- read.csv('dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.combined_for_r.txt', sep = '\t')
goodChrOrder <- paste("chr",c(1:25),sep="")
hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
hom_data
#ggplot(data=hom_data,aes(x=start, y=average_aaf, fill=test),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid( ~ chr, scales="free", space="free_x") + 
#  theme_bw() + theme(axis.text.x = element_blank())
ggplot(data=hom_data,aes(x=start, y=value, group=test, color=test)) +
  geom_line() + facet_grid(test ~ chr, scales="free", space="free_x") + theme(axis.text.x = element_blank())

##loop through files with aaf values
filename <- dir('/home/atimms/ngs_data/enu_mapping/dave_zf_sex_0920', pattern ="aaf_for_r.txt")
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
filename <- dir('/home/atimms/ngs_data/enu_mapping/dave_zf_sex_0920', pattern ="combined_for_r.txt")
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
  ggsave(gene_pdf, width = 20, height = 5)
}

##just graph chr20 and chr21
##combined values
hom_data <- read.csv('dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.combined_for_r.txt', sep = '\t')
hom_data21 = hom_data %>% filter(chr == "chr21")
hom_data20 = hom_data %>% filter(chr == "chr20")
ggplot(data=hom_data20,aes(x=start, y=value, group=test, color=test)) +
  geom_line() + facet_grid(test ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.combined_for_r.chr20.pdf', width = 20, height = 5)
ggplot(data=hom_data21,aes(x=start, y=value, group=test, color=test)) +
  geom_line() + facet_grid(test ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.combined_for_r.chr21.pdf', width = 20, height = 5)
hom_data <- read.csv('dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.combined_for_r.txt', sep = '\t')
hom_data21 = hom_data %>% filter(chr == "chr21")
hom_data20 = hom_data %>% filter(chr == "chr20")
#ggplot(data=hom_data,aes(x=start, y=average_aaf, fill=test),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid( ~ chr, scales="free", space="free_x") + 
#  theme_bw() + theme(axis.text.x = element_blank())
ggplot(data=hom_data20,aes(x=start, y=value, group=test, color=test)) +
  geom_line() + facet_grid(test ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.combined_for_r.chr20.pdf', width = 20, height = 5)
ggplot(data=hom_data21,aes(x=start, y=value, group=test, color=test)) +
  geom_line() + facet_grid(test ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.combined_for_r.chr21.pdf', width = 20, height = 5)
##just aaf
hom_data <- read.csv('dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.aaf_for_r.txt', sep = '\t')
hom_data21 = hom_data %>% filter(chr == "chr21")
hom_data20 = hom_data %>% filter(chr == "chr20")
ggplot(data=hom_data20,aes(x=start, y=average_aaf, group=test, color=test)) +
  geom_line() + facet_grid( ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.aaf_for_r.chr20.pdf', width = 20, height = 5)
ggplot(data=hom_data21,aes(x=start, y=average_aaf, group=test, color=test)) +
  geom_line() + facet_grid( ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.aaf_for_r.chr21.pdf', width = 20, height = 5)
hom_data <- read.csv('dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.aaf_for_r.txt', sep = '\t')
hom_data21 = hom_data %>% filter(chr == "chr21")
hom_data20 = hom_data %>% filter(chr == "chr20")
ggplot(data=hom_data20,aes(x=start, y=average_aaf, group=test, color=test)) +
  geom_line() + facet_grid( ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.aaf_for_r.chr20.pdf', width = 20, height = 5)
ggplot(data=hom_data21,aes(x=start, y=average_aaf, group=test, color=test)) +
  geom_line() + facet_grid( ~ chr, scales="free", space="free_x") + scale_x_continuous(labels = comma)
ggsave('dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.danRer11_2000kb_500kb.aaf_for_r.chr21.pdf', width = 20, height = 5)


