##install archr and related packages
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
#library(ArchR)
#ArchR::installExtraPackages()
#devtools::install_github("GreenleafLab/chromVARmotifs")
##install GenomicRanges to analyze the co-accessability data
#if (!require("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#devtools::install_github("immunogenomics/harmony")
#devtools::install_github("immunogenomics/presto")
#install.packages("magick")

##get archr and set up env
#library(ggplot2)
#library(GenomicRanges)
library(chromVARmotifs)
library(pheatmap)
library(ArchR)
library(ComplexHeatmap)
library(circlize)
set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")
##Set working directory and specific dir in active RSS where the data is sitting
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all_third')


##load archr project
proj1 = loadArchRProject(path = 'human_all')

##call peaks and then retrieve
##make pseudo-bulk replicates - change minCells to 20 (for mature ganglions)
proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "Clusters2", minCells = 20, force = TRUE)
##call peaks with macs2
pathToMacs2 <- findMacs2()
proj1 <- addReproduciblePeakSet(ArchRProj = proj1, groupBy = "Clusters2", 
                                pathToMacs2 = pathToMacs2, cutOff = 0.000001)
##get GRanges object
getPeakSet(proj1)
##add peak matrix
proj1 <- addPeakMatrix(proj1)
getAvailableMatrices(proj1)
##Identifying Marker peaks between cell classes
#get peaks specific to cell class
markersPeaks <- getMarkerFeatures(ArchRProj = proj1, 
                                  useMatrix = "PeakMatrix", groupBy = "Clusters2",
                                  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")

##graph marker peaks
heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, 
                                  cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmap_cell_class_marker_peaks.macs2_q0.000001.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)
##instead of drawing plot make a matrix
heatmap.matrix <- plotMarkerHeatmap(seMarker = markersPeaks, 
                                    cutOff = "FDR <= 0.01 & Log2FC >= 1", returnMat = TRUE)
##get csv file from matrix
write.csv(heatmap.matrix, file='heatmap_cell_class_marker_peaks.macs2_q0.000001.matrix.csv')
##get all peaks 
heatmap.matrix <- plotMarkerHeatmap(seMarker = markersPeaks, 
                                    cutOff = "FDR <= 1000 & Log2FC >= -1000", returnMat = TRUE)
##get csv file from matrix
write.csv(heatmap.matrix, file='heatmap_cell_class_all_peaks.macs2_q0.000001.matrix.csv')

##get peak to gene data and atac/rna index info - repeat and save in current dir
proj1 <- addPeak2GeneLinks(ArchRProj = proj1,reducedDims = "Harmony", maxDist = 300000)
p2g <- getPeak2GeneLinks(ArchRProj = proj1, corCutOff = 0.4,
                         resolution = 10000, returnLoops = F)
write.table(p2g, 'human_all.min_cell20.p2g.m2_q0.000001.cor4.txt', row.names = F, sep="\t", quote = FALSE)
rnaidx.info = metadata(p2g)[[2]]
rnaidx.df = data.frame(rnaidx.info)
write.table(rnaidx.df, 'human_all.min_cell20.m2_q0.000001.rnaseq_info.txt', row.names = T, sep="\t", quote = FALSE)
peaks.gr <- getPeakSet(proj1)
df.peaks.gr = data.frame(peaks.gr)
write.table(df.peaks.gr, 'human_all.min_cell20.peak_info.m2_q0.000001.txt', row.names = T, sep="\t", quote = FALSE)

