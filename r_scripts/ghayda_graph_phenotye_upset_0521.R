#install.packages("UpSetR")
#install.packages('ComplexUpset')
# Library
##using r4.0.3
library(UpSetR)
library(ComplexUpset) ##if loaded can't use upset in upsetr
library(ggplot2)
##https://github.com/hms-dbmi/UpSetR
##https://jku-vds-lab.at/tools/upset/#:~:text=UpSet%20concept,the%20figure%20on%20the%20right.&text=The%20first%20row%20in%20the,B%20or%20C)%2C%20etc.
##https://github.com/krassowski/complex-upset

setwd('/archive/mirzaa_g/exomes/result_files_0321/graphing_0321')

##read in data -- adding summary meg/mic to phenotypes
links <- read.table('upset_test_0521.txt', header=T, row.names = 1, sep='\t')

##graph interactions
upset(links, sets = c("MEG", "MIC", "WM", "CTX", "Ventricle", "CBL", "BS", "BG_TH", "CC"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.copy2pdf(file='just_phenotypes.upset.052721.pdf', width = 9, height = 6)



##fill in with solved rates etc
##to get stacked bar plots need ComplexUpset and maybe ggplot
#https://stackoverflow.com/questions/54770795/stacked-barplot-in-upsetr/56704255#56704255

##read in data -- adding summary meg/mic to phenotypes
links <- read.table('upset_test_0521.txt', header=T, row.names = 1, sep='\t')
##get names of columns to get intersections
genres = colnames(links)[3:10]
##make stacked graph and copy
upset(links,genres,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=summary)
    )
  ),
  width_ratio=0.1)
dev.copy2pdf(file='phenotypes_summary.upset.052721.pdf', width = 15, height = 9)
##make stacked graph with 
upset(links,genres,
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=summary_combined)
        )
      ),
      width_ratio=0.1)
dev.copy2pdf(file='phenotypes_summary_combined.upset.052721.pdf', width = 15, height = 9)

##make graph and copy
upset(links,genres,
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=gene)
        )
      ),
      width_ratio=0.1)
dev.copy2pdf(file='phenotypes_genes.upset.052721.pdf', width = 15, height = 9)
