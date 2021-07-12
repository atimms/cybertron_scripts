#install.packages("UpSetR")
#install.packages('ComplexUpset')
# Library
##using r4.0.3
library(UpSetR)
#library(ComplexUpset) ##used for stacked upset
library(ggplot2)

setwd('/archive/mirzaa_g/exomes/result_files_0321/graphing_0321')

##upset graphs for phenotypes

##read in data -- adding summary meg/mic to phenotypes 0/1 data
links <- read.table('upset_0621_1.txt', header=T, row.names = 1, sep='\t')
##graph interactions
upset(links, sets = c('MEG', 'MIC', 'WM', 'CTX', 'Ventricle', 'CBL', 'BS', 'BG.TH', 'CC'), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.copy2pdf(file='phenotypes_1.upset.060921.pdf', width = 9, height = 6)

##read in data -- adding summary meg/mic to phenotypes using new split data
links <- read.table('upset_0621_2.txt', header=T, row.names = 1, sep='\t')
##graph interactions
upset(links, sets = c('MEG', 'MIC', 'LIS', 'PMG', 'SIMP', 'CTXD', 'HET', 'DWM', 'non.DWM', 'MTM', 'CBTE', 'ACC', 'pACC', 'thinCC', 'MegaCC'), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.copy2pdf(file='phenotypes_2.upset.060921.pdf', width = 9, height = 6)


##solved graphs

##solved/brain size
dat = read.table('solved.brain_size.txt', header=T, row.names = 1, sep='\t')
dat$absolute_SD = abs(dat$SD) ##get absolute values
##comapare brainsize/summary boxplots
ggplot(dat, aes(x=summary, y=absolute_SD, fill=summary)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust=1))  + guides(fill=FALSE)
ggsave('brain_size_summary.boxplot.060921.pdf', width=10, height = 6)
ggplot(dat, aes(x=solved, y=absolute_SD, fill=solved)) + geom_boxplot() + guides(fill=FALSE)
ggsave('brain_size_solved.boxplot.060921.pdf', width=10, height = 6) 
ggplot(dat, aes(x=solved2, y=absolute_SD, fill=solved2)) + geom_boxplot() + guides(fill=FALSE)
ggsave('brain_size_solved2.boxplot.060921.pdf', width=10, height = 6)
##density/histograms
ggplot(dat, aes(x=SD, fill=solved2)) + geom_density(alpha=.3)
ggsave('brain_size_solved2.density.060921.pdf', width=10, height = 6)

##solved/dx1
dat = read.table('solved.dx1.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=DxGroup1, fill=summary)) +
  geom_bar(stat="count", position="fill")
ggsave('DxGroup1_summary.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=DxGroup1, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('DxGroup1_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=DxGroup1, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('DxGroup1_solved2.stacked_bar.060921.pdf', width=10, height = 6)

##solved/meg or mic
dat = read.table('solved.meg_mic.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=MEG_MIC, fill=summary)) +
  geom_bar(stat="count", position="fill")
ggsave('meg_mic_summary.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=MEG_MIC, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('meg_mic_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=MEG_MIC, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('meg_mic_solved2.stacked_bar.060921.pdf', width=10, height = 6)

##solved/meg or mic
dat = read.table('solved.severe_meg_mic.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=MEG_MIC, fill=summary)) +
  geom_bar(stat="count", position="fill")
ggsave('severe_meg_mic_summary.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=MEG_MIC, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('severe_meg_mic_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=MEG_MIC, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('severe_meg_mic_solved2.stacked_bar.060921.pdf', width=10, height = 6)

##solved/ctx
dat = read.table('solved.ctx.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=CTX, fill=summary)) +
  geom_bar(stat="count", position="fill")
ggsave('ctx_summary.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=CTX, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('ctx_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=CTX, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('ctx_solved2.stacked_bar.060921.pdf', width=10, height = 6)

##solved/het
dat = read.table('solved.het.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=HET, fill=summary)) +
  geom_bar(stat="count", position="fill")
ggsave('het_summary.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=HET, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('het_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=HET, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('het_solved2.stacked_bar.060921.pdf', width=10, height = 6)

##solved/cbl
dat = read.table('solved.cbl.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=CBL, fill=summary)) +
  geom_bar(stat="count", position="fill")
ggsave('cbl_summary.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=CBL, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('cbl_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=CBL, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('cbl_solved2.stacked_bar.060921.pdf', width=10, height = 6)

##solved/cc
dat = read.table('solved.cc.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=CC, fill=summary)) +
  geom_bar(stat="count", position="fill")
ggsave('cc_summary.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=CC, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('cc_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=CC, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('cc_solved2.stacked_bar.060921.pdf', width=10, height = 6)

##solved/ head size group
dat = read.table('solved.hs_group.txt', header=T, row.names = 1, sep='\t')
ggplot(data=dat, aes(x=SD.group, fill=solved)) +
  geom_bar(stat="count", position="fill")
ggsave('head_size_solved.stacked_bar.060921.pdf', width=10, height = 6)
ggplot(data=dat, aes(x=SD.group, fill=solved2)) +
  geom_bar(stat="count", position="fill")
ggsave('head_size_solved2.stacked_bar.060921.pdf', width=10, height = 6)
