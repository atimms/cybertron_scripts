# Bioconductor version (stable)
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MAGeCKFlute")
library(MAGeCKFlute)
library(ggplot2)

##Set working directory and specific dir in active RSS where the data is sitting
setwd('/home/atimms/ngs_data/misc/vishal_crispr_screen_1020')

##from mageck mle analysis
file3 = "./vishal_crispr_1020.mageck_mle.gene_summary.txt"
# Read and visualize the file format
gdata = ReadBeta(file3)
head(gdata)
##run pipeline
FluteMLE(gdata, treatname="LH2_top10", ctrlname="LH3_bottom25", proj="LH2_mle", organism="hsa")
FluteMLE(gdata, treatname="LH3_bottom25", ctrlname="LH2_top10", proj="LH3_mle", organism="hsa")

##magic test on LH2
file1 = file.path("./vishal_crispr_1020.mageck_test_lh2.gene_summary.txt")
# Read and visualize the file format
gdata = ReadRRA(file1)
head(gdata)
file2 = file.path("./vishal_crispr_1020.mageck_test_lh2.sgrna_summary.txt")
sdata = ReadsgRRA(file2)
head(sdata)
FluteRRA(gdata, sdata, proj="LH2_RRA_default", organism="hsa",
         scale_cutoff = 1, outdir = "./")
FluteRRA(gdata, proj="LH2_RRA_depmap_omit", organism="hsa", incorporateDepmap = TRUE,
         omitEssential = TRUE, outdir = "./")

##magic test on LH3
file1 = file.path("./vishal_crispr_1020.mageck_test_lh3.gene_summary.txt")
# Read and visualize the file format
gdata = ReadRRA(file1)
head(gdata)
file2 = file.path("./vishal_crispr_1020.mageck_test_lh3.sgrna_summary.txt")
sdata = ReadsgRRA(file2)
head(sdata)
FluteRRA(gdata, sdata, proj="LH3_RRA_default", organism="hsa",
         scale_cutoff = 1, outdir = "./")
FluteRRA(gdata, proj="LH3_RRA_depmap_omit", organism="hsa", incorporateDepmap = TRUE,
         omitEssential = TRUE, outdir = "./")







##step by step
lh2 = 'LH2_top10'
lh3 = 'LH3_bottom25'
##Normalization of beta scores
gdata_cc = NormalizeBeta(gdata, samples=c(lh2, lh3), method="cell_cycle")
head(gdata_cc)
##Distribution of all gene beta scores
DensityView(gdata_cc, samples=c(lh2, lh3))
dev.copy2pdf(file='vishal_crispr_1020.mageck_mle.beta_dist.pdf', width = 7, height = 5)
ConsistencyView(gdata_cc, lh2, lh3)
dev.copy2pdf(file='vishal_crispr_1020.mageck_mle.beta_consistancy.pdf', width = 7, height = 5)
MAView(gdata_cc, lh2, lh3)
dev.copy2pdf(file='vishal_crispr_1020.mageck_mle.beta_MAView.pdf', width = 7, height = 5)
##positive and negative selection
p1 = ScatterView(gdata_cc, "LH2_top10", "LH3_bottom25", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE)
print(p1)
dev.copy2pdf(file='vishal_crispr_1020.mageck_mle.beta_top_bottom.pdf', width = 7, height = 5)
##rankplot
rankdata = gdata_cc$LH2_top10 - gdata_cc$LH3_bottom25
names(rankdata) = gdata_cc$Gene
RankView(rankdata)
dev.copy2pdf(file='vishal_crispr_1020.mageck_mle.beta_rank.pdf', width = 7, height = 5)
gdata_cc$Diff = gdata_cc$LH2_top10 - gdata_cc$LH3_bottom25
gdata_cc$Rank = rank(gdata_cc$Diff)
p1 = ScatterView(gdata_cc, x = "Diff", y = "Rank", label = "Gene", 
                 top = 5, model = "rank")
print(p1)
dev.copy2pdf(file='vishal_crispr_1020.mageck_mle.beta_rank2.pdf', width = 7, height = 5)
##square plot
p1 = ScatterView(gdata_cc, x = "LH2_top10", y = "LH3_bottom25", label = "Gene", 
                 model = "ninesquare", top = 10, display_cut = TRUE, force = 2)
print(p1)
dev.copy2pdf(file='vishal_crispr_1020.mageck_mle.beta_square.pdf', width = 7, height = 5)
# 9-square groups
Square9 = p1$data
idx=Square9$group=="midleft"
geneList = Square9$Diff
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene
# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 6, bottom = 0)

