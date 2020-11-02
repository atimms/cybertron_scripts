# Load the ggplot2 package
library(ggplot2)
##set working directory
setwd('/home/atimms/ngs_data/misc/eric_graph_SH2B_0920')

##Human Bone marrow, from noramalized data
dat <- read.csv('human_bm.SH2B_formatted.txt', sep = '\t')

##facetted boxplot
ggplot(dat, aes(x=tissue, y=value, fill=gene)) + geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('human_bm.SH2B.boxplot.pdf', width = 10, height = 8)

##Human blood, from TPM noramalized data
dat <- read.csv('human_blood_GSE107011.TPM.txt', sep = '\t')
##facetted boxplot
ggplot(dat, aes(x=tissue, y=TPM, fill=gene)) + geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size =7)) +
  facet_grid(gene ~ .) + guides(fill=FALSE)
ggsave('human_blood_GSE107011.SH2B.boxplot.pdf', width = 10, height = 8)
