#install.packages("UpSetR")
#install.packages('ComplexUpset')
# Library
##using r4.0.3
library(UpSetR)
#library(ComplexUpset) ##used for stacked upset
library(ggplot2)
#library("viridis") ##for color palette, not used
library(RColorBrewer)
library(stringr) #for string manipulation
library(tidyr)

setwd('/archive/mirzaa_g/exomes/result_files_0321/graphing/graphing_0721')


##for changing colors on upset
#https://krassowski.github.io/complex-upset/articles/Examples_R.html#5-adjusting-other-aesthetics

##color palettes
#https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
##using RColorBrewer (which i must load)
##look at palettes for color blind
display.brewer.all(colorblindFriendly = TRUE)


##read in all factors
f2g <- read.table('factors_to_graph_070921.txt', header=T, sep='\t')
f2g


####classification i.e. ped types

##all ped types
##random order
ggplot(f2g,aes(y = classification, fill=classification)) + geom_bar() + theme_minimal()
##order by frequency
ggplot(f2g,aes(y = reorder(classification, table(classification)[classification]), fill=classification)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("classification")
ggsave("pedigree_types_all.barchart.071221.pdf", width = 9, height = 6)
##stacked by dxgroup
ggplot(f2g,aes(y = reorder(classification, table(classification)[classification]), fill=DxGroup1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("classification")
ggsave("pedigree_types_all.stacked_barchart_dx.071221.pdf", width = 9, height = 6)

##just the mosaic peds
f2g_mosaic = subset(f2g, classification %in% c("Trio*", "Singleton*", "Duo*"))
ggplot(f2g_mosaic,aes(y = reorder(classification, table(classification)[classification]), fill=classification)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("classification")
ggsave("pedigree_types_mosaic.barchart.071221.pdf", width = 9, height = 6)
ggplot(f2g_mosaic,aes(y = reorder(classification, table(classification)[classification]), fill=DxGroup1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("classification")
ggsave("pedigree_types_mosaic.stacked_barchart_dx.071221.pdf", width = 9, height = 6)

##combine the mosaic peds to regular peds
##remove the asterix
f2g$classification_no_mosaic = str_replace_all(f2g$classification,"\\*","")
##order by frequency
ggplot(f2g,aes(y = reorder(classification_no_mosaic, table(classification_no_mosaic)[classification_no_mosaic]), fill=classification_no_mosaic)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("classification")
ggsave("pedigree_types_no_mosaic.barchart.071221.pdf", width = 9, height = 6)
##stacked by dxgroup
ggplot(f2g,aes(y = reorder(classification_no_mosaic, table(classification_no_mosaic)[classification_no_mosaic]), fill=DxGroup1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("classification")
ggsave("pedigree_types_no_mosaic.stacked_barchart_dx.071221.pdf", width = 9, height = 6)


#Ethnicity
##known... split out column with commas
f2g$eth1 = str_replace_all(f2g$NIH.Ethnicity, ",.*", "")
##order by frequency
ggplot(f2g,aes(y = reorder(eth1, table(eth1)[eth1]), fill=eth1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("NIH Ethnicity")
ggsave("NIH_ethnicity.barchart.071221.pdf", width = 9, height = 6)
##stacked by dxgroup
ggplot(f2g,aes(y = reorder(eth1, table(eth1)[eth1]), fill=DxGroup1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("NIH Ethnicity")
ggsave("NIH_ethnicity.stacked_barchart_dx.071221.pdf", width = 9, height = 6)
##with laser
f2g$eth2 = str_replace_all(f2g$laser.assisted.ethnicity, ",.*", "")##order by frequency
ggplot(f2g,aes(y = reorder(eth2, table(eth2)[eth2]), fill=eth2)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("laser assisted ethnicity")
ggsave("laser_ethnicity.barchart.071221.pdf", width = 9, height = 6)
##stacked by dxgroup
ggplot(f2g,aes(y = reorder(eth2, table(eth2)[eth2]), fill=DxGroup1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("laser assisted ethnicity")
ggsave("laser_ethnicity.stacked_barchart_dx.071221.pdf", width = 9, height = 6)

#Gender
##known... split out column with commas
f2g_sex = separate_rows(f2g,"Gender.",sep = ", ")
f2g_sex$gender = str_replace_all(f2g_sex$Gender.," ","")
##order by frequency
ggplot(f2g_sex,aes(y = reorder(gender, table(gender)[gender]), fill=gender)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("Gender")
ggsave("gender.barchart.071221.pdf", width = 9, height = 6)
##stacked by dxgroup
ggplot(f2g_sex,aes(y = reorder(gender, table(gender)[gender]), fill=DxGroup1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("Gender")
ggsave("gender.stacked_barchart_dx.071221.pdf", width = 9, height = 6)

##sample types
f2g_st = separate_rows(f2g,"proband_sample_type",sep = ", ")
##remove parentheses
f2g_st$sample_type = str_replace_all(f2g_st$proband_sample_type, " \\s*\\([^\\)]+\\)", "")
##order by frequency
ggplot(f2g_st,aes(y = reorder(sample_type, table(sample_type)[sample_type]), fill=sample_type)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("sample type")
ggsave("sample_type.barchart.071221.pdf", width = 9, height = 6)
##stacked by dxgroup
ggplot(f2g_st,aes(y = reorder(sample_type, table(sample_type)[sample_type]), fill=DxGroup1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("sample type")
ggsave("sample_type.stacked_barchart_dx.071221.pdf", width = 9, height = 6)

##sample number i.e. multiple samples/proband and stacked by dxgroup
ggplot(f2g,aes(reorder(DxGroup1, -table(DxGroup1)[DxGroup1]), fill=factor(sample_number, levels=c( "2 or more samples", "1 sample")))) + geom_bar() + scale_fill_brewer(palette = "Set2") + theme_minimal() + xlab("DxGroup1") + guides(fill=guide_legend(title="sample counts"))
ggsave("sample_counts.stacked_barchart_dx.071221.pdf", width = 9, height = 6)


###coverage data

##read in data
cov_all <- read.table('coverage_all_peds.txt', header=T, sep='\t')
cov_all
cov_brain_peds <- read.table('coverage_brain_peds.txt', header=T, sep='\t')
cov_brain_peds
##add extra cols for ped_types
cov_all$ped_type_comb = str_replace_all(cov_all$ped_type,"\\*","")

##make box plots
ggplot(cov_all, aes(x=dx, y=coverage, fill=dx)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal() 
ggsave("coverage.all.dx1.071221.pdf", width = 9, height = 6)
ggplot(cov_all, aes(x=ped_type, y=coverage, fill=ped_type)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal() 
ggsave("coverage.all.all_ped_types.071221.pdf", width = 9, height = 6)
ggplot(cov_all, aes(x=ped_type_comb, y=coverage, fill=ped_type_comb)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal() 
ggsave("coverage.all.combined_ped_types.071221.pdf", width = 9, height = 6)
ggplot(cov_all, aes(x=sample_type, y=coverage, fill=sample_type)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal() 
ggsave("coverage.all.sample_type.071221.pdf", width = 9, height = 6)
ggplot(cov_brain_peds, aes(x=sample_type_alt, y=coverage, fill=sample_type_alt)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal() 
ggsave("coverage.brain_samples.sample_type.071221.pdf", width = 9, height = 6)

##upset graphs for phenotypes

##read in data -- adding summary meg/mic to phenotypes 0/1 data
links <- read.table('upset1.txt', header=T, row.names = 1, sep='\t')
##graph interactions
upset(links, sets = c('MEG', 'MIC', 'CTX', 'CBL', 'BS', 'BG.TH', 'CC'), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.copy2pdf(file='phenotypes_1.upset.071221.pdf', width = 9, height = 6)

##read in data -- adding summary meg/mic to phenotypes using new split data
links <- read.table('upset2.txt', header=T, row.names = 1, sep='\t')
##graph interactions
upset(links, sets = c('MEG', 'MIC', 'CTXD', 'LIS', 'PMG', 'CTX_other', 'SIMP', 'HET', 'DWM', 'non.DWM', 'CBTE', 'ACC', 'megaCC'), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.copy2pdf(file='phenotypes_2.upset.071221.pdf', width = 9, height = 6)

##assymetry data

##read in data
asy_data <- read.table('asymmetry_072221.txt', header=T, sep='\t')
asy_data

ggplot(asy_data,aes(x = asymmetry, fill=summary)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("summary")
ggsave("asymmetry_summary.stacked_barchart.072221.pdf", width = 9, height = 6)
ggplot(asy_data,aes(x = asymmetry, fill=solved1)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved")
ggsave("asymmetry_solved1.stacked_barchart.072221.pdf", width = 9, height = 6)
ggplot(asy_data,aes(x = asymmetry, fill=solved2)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved")
ggsave("asymmetry_solved2.stacked_barchart.072221.pdf", width = 9, height = 6)
ggplot(asy_data,aes(x = asymmetry, fill=DxGroup1)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal()
ggsave("asymmetry_DxGroup1.stacked_barchart.072221.pdf", width = 9, height = 6)
ggplot(asy_data,aes(x = mosaic, fill=asymmetry)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal()
ggsave("asymmetry_mosaic.stacked_barchart.072221.pdf", width = 9, height = 6)


##solved graphs

##read in data
asy_dat <- read.table('asymmetry_072221.txt', header=T, sep='\t')
asy_dat

##all ped types
##random order
ggplot(solved_dat,aes(y = SUMMARY..interpretation., fill=SUMMARY..interpretation.)) + geom_bar() + theme_minimal()
##order by frequency
ggplot(solved_dat,aes(y = reorder(SUMMARY..interpretation., table(SUMMARY..interpretation.)[SUMMARY..interpretation.]), fill=SUMMARY..interpretation.)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("summary")
ggsave("summary.barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat,aes(y = reorder(solved1, table(solved1)[solved1]), fill=solved1)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved")
ggsave("solved1.barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat,aes(y = reorder(solved2, table(solved2)[solved2]), fill=solved2)) + geom_bar() + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved")
ggsave("solved2.barchart.071321.pdf", width = 9, height = 6)

##stacked by dxgroup
ggplot(solved_dat,aes(x = DxGroup1, fill=SUMMARY..interpretation.)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("summary")
ggsave("summary_dx1.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat,aes(x = DxGroup1, fill=solved1)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved")
ggsave("solved1_dx1.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat,aes(x = DxGroup1, fill=solved2)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved")
ggsave("solved2_dx1.stacked_barchart.071321.pdf", width = 9, height = 6)

##ctx
##remove ND from ctx column
solved_dat_ctx= subset(solved_dat, CTX..LIS.PMG.SIMP.CTXD.NL. != "ND")
##plot
ggplot(solved_dat_ctx,aes(x = CTX..LIS.PMG.SIMP.CTXD.NL., fill=SUMMARY..interpretation.)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("summary") + xlab("CTX")
ggsave("summary_ctx.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_ctx,aes(x = CTX..LIS.PMG.SIMP.CTXD.NL., fill=solved1)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("CTX")
ggsave("solved1_ctx.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_ctx,aes(x = CTX..LIS.PMG.SIMP.CTXD.NL., fill=solved2)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("CTX")
ggsave("solved2_ctx.stacked_barchart.071321.pdf", width = 9, height = 6)

##cbl
##remove ND from cbl column
solved_dat_1= subset(solved_dat, CBL..DWW.non.DWM.NL.CBTE. != "ND")
solved_dat_cbl= subset(solved_dat_1, CBL..DWW.non.DWM.NL.CBTE. != "no data")
##plot
ggplot(solved_dat_cbl,aes(x = CBL..DWW.non.DWM.NL.CBTE., fill=SUMMARY..interpretation.)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("summary") + xlab("CBL")
ggsave("summary_cbl.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_cbl,aes(x = CBL..DWW.non.DWM.NL.CBTE., fill=solved1)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("CBL")
ggsave("solved1_cbl.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_cbl,aes(x = CBL..DWW.non.DWM.NL.CBTE., fill=solved2)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("CBL")
ggsave("solved2_cbl.stacked_barchart.071321.pdf", width = 9, height = 6)

##cc
##remove ND from cc column
solved_dat_cc= subset(solved_dat, CC..ACC.pACC.thin.NL.MegaCC. != "ND")
##plot
ggplot(solved_dat_cc,aes(x = CC..ACC.pACC.thin.NL.MegaCC., fill=SUMMARY..interpretation.)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("summary") + xlab("CC")
ggsave("summary_cc.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_cc,aes(x = CC..ACC.pACC.thin.NL.MegaCC., fill=solved1)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("CC")
ggsave("solved1_cc.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_cc,aes(x = CC..ACC.pACC.thin.NL.MegaCC., fill=solved2)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("CC")
ggsave("solved2_cc.stacked_barchart.071321.pdf", width = 9, height = 6)

##brain_size
##remove ND from brain_size column
solved_dat_bs= subset(solved_dat, brain_size != "ND")
##plot
ggplot(solved_dat_bs,aes(x = brain_size, fill=SUMMARY..interpretation.)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("summary") + xlab("brain size")
ggsave("summary_brain_size.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_bs,aes(x = brain_size, fill=solved1)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("brain size")
ggsave("solved1_brain_size.stacked_barchart.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_bs,aes(x = brain_size, fill=solved2)) + geom_bar(position="fill") + scale_fill_brewer(palette = "RdYlBu") + theme_minimal() + ylab("solved") + xlab("brain size")
ggsave("solved2_brain_size.stacked_barchart.071321.pdf", width = 9, height = 6)

##brain_size using size sds
##remove ND/NL from brain_size column and get abs values
solved_dat_1= subset(solved_dat, brain_size_sd != "ND")
solved_dat_sd= subset(solved_dat_1, brain_size_sd != "NL")
solved_dat_sd$abs_sd = abs(as.numeric(solved_dat_sd$brain_size_sd))
#boxplots
ggplot(solved_dat_sd, aes(x=SUMMARY..interpretation., y=abs_sd, fill=SUMMARY..interpretation.)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal()  + ylab("absolute SD")  + xlab("summary") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("summary_brain_sd.boxplot.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_sd, aes(x=solved1, y=abs_sd, fill=solved1)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal()  + ylab("absolute SD")  + xlab("solved") 
ggsave("solved1_brain_sd.boxplot.071321.pdf", width = 9, height = 6)
ggplot(solved_dat_sd, aes(x=solved2, y=abs_sd, fill=solved2)) + geom_boxplot()+ scale_fill_brewer(palette = "RdYlBu") + theme_minimal()  + ylab("absolute SD")  + xlab("solved") 
ggsave("solved2_brain_sd.boxplot.071321.pdf", width = 9, height = 6)



