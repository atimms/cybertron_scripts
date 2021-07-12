# install.packages("dendsort")

##load libraries
library("pheatmap")
library("ggplot2")
library('dendsort')

workingDir = "/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_gwas_atac_peaks_0621";
setwd(workingDir);

##graph snp data
human_hm <- read.table('combined_snps.cc_peaks.counts_for_heatmap.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm)
dev.copy2pdf(file='combined_snps.cc_peaks.heatmap_070821.pdf', width = 9, height = 6)

##graph region data
human_hm <- read.table('combined_gwas_regions.cc_peaks.for_heatmap.count.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm, fontsize_row=4, fontsize_col = 5)
dev.copy2pdf(file='combined_gwas_regions.cc_peaks.for_heatmap.count.070821.pdf', width = 9, height = 6)
human_hm <- read.table('combined_gwas_regions.cc_peaks.for_heatmap.count_snp.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm, fontsize_row=4, fontsize_col = 5)
dev.copy2pdf(file='combined_gwas_regions.cc_peaks.for_heatmap.count_snp.070821.pdf', width = 9, height = 6)
human_hm <- read.table('combined_gwas_regions.cc_peaks.for_heatmap.count_size.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm, fontsize_row=4, fontsize_col = 5)
dev.copy2pdf(file='combined_gwas_regions.cc_peaks.for_heatmap.count_size.070821.pdf', width = 9, height = 6)


##graph lead snp data
human_hm <- read.table('gwas_lead_snp_per_region_0721.cc_peaks.counts_for_heatmap.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm, fontsize_row=3, fontsize_col = 5)
dev.copy2pdf(file='gwas_lead_snp_per_region_0721.cc_peaks.counts_for_heatmap.pdf', width = 9, height = 6)
human_hm <- read.table('gwas_lead_snp_per_set_0721.cc_peaks.counts_for_heatmap.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm, fontsize_row=2, fontsize_col = 5)
dev.copy2pdf(file='gwas_lead_snp_per_set_0721.cc_peaks.counts_for_heatmap.pdf', width = 9, height = 6)

