##load libraries needed
library(ggplot2)
library(gridBase)
library(gridExtra)   

##dave/chunyue_abcb11b_recessive_0321
##homozygosity mapping
setwd('/home/atimms/ngs_data/enu_mapping/chunyue_abcb11b_recessive_0321')
##loop through files
filename <- dir('/home/atimms/ngs_data/enu_mapping/chunyue_abcb11b_recessive_0321', pattern ="combined_hom_mapping.txt")
for(i in 1:length(filename)){
  ##read in hom data
  hom_data <- read.csv(filename[i], sep ='\t')
  ##reorder chromosome, if using 'chr'
  #goodChrOrder <- paste("chr",c(1:25,"X","Y"),sep="")
  #hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  #goodChrOrder <- paste(c(1:25,"X","Y"),sep="")
  #hom_data$chr <- factor(hom_data$chr,levels=goodChrOrder)
  ##bar graph
  ggplot(data=hom_data,aes(x=chromosome, y=value),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(analysis ~ chr, scales="free", space="free_x") + 
    theme_bw() + theme(axis.text.x = element_blank())
  gene_pdf <- gsub(".txt",".pdf",filename[i])
  ggsave(gene_pdf, height = 4, width = 12)
  #ggsave(paste(filename[i], 'pdf',sep="."))
}
##test just graphing one at a time
hom_data <- read.csv('TP1FUnr_danRer11_1000kb_100kb_combined_hom_mapping.txt', sep ='\t')
ggplot(data=hom_data,aes(x=chromosome, y=value),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(analysis ~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('test.pdf', height = 4, width = 12)

##graph the count data
hom_data <- read.csv('TP1FRes.resc_hom.danRer11_1000kb_100kb.bed', sep ='\t')
ggplot(data=hom_data,aes(x=start, y=snp_number),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('TP1FRes.resc_hom.danRer11_1000kb_100kb.pdf', height = 4, width = 12)
hom_data <- read.csv('TP1FRes.resc_hom.danRer11_100kb_100kb.bed', sep ='\t')
ggplot(data=hom_data,aes(x=start, y=snp_number),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('TP1FRes.resc_hom.danRer11_100kb_100kb.pdf', height = 4, width = 12)
hom_data <- read.csv('TP1FRes.resc_hom.danRer11_500kb_100kb.bed', sep ='\t')
ggplot(data=hom_data,aes(x=start, y=snp_number),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('TP1FRes.resc_hom.danRer11_500kb_100kb.pdf', height = 4, width = 12)

##graph the pseudo snptrack data
hom_data <- read.csv('TP1FRes_danRer11_1000kb_100kb.snptrack_hom.bed', sep ='\t')
ggplot(data=hom_data,aes(x=start, y=hom_score_rescued),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('TP1FRes_danRer11_1000kb_100kb_snptrack_hom.pdf', height = 4, width = 12)
ggplot(data=hom_data,aes(x=start, y=hom_score_all),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('TP1FRes_danRer11_1000kb_100kb_snptrack_hom_alt.pdf', height = 4, width = 12)
hom_data <- read.csv('TP1FRes_danRer11_100kb_10kb.snptrack_hom.bed', sep ='\t')
ggplot(data=hom_data,aes(x=start, y=hom_score_rescued),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('TP1FRes_danRer11_100kb_10kb_snptrack_hom.pdf', height = 4, width = 12)
ggplot(data=hom_data,aes(x=start, y=hom_score_all),binwidth=1e5 ) + geom_bar(stat="identity", colour="black") + facet_grid(~ chr, scales="free", space="free_x") + 
  theme_bw() + theme(axis.text.x = element_blank())
ggsave('TP1FRes_danRer11_100kb_10kb_snptrack_hom_alt.pdf', height = 4, width = 12)

