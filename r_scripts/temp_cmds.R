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

##get archr and set up env
library(pheatmap)
library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")
##Set working directory and specific dir in active RSS where the data is sitting
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all')

##load project
proj1 = loadArchRProject(path = 'human')

##make pseudo-bulk replicates - changed min cells
proj2 <- addGroupCoverages(ArchRProj = proj1, groupBy = "Clusters2", minCells = 20, force = TRUE)
saveArchRProject(ArchRProj = proj2, outputDirectory = "0920_redo_macs2_q0.01")

##call peaks with macs2
pathToMacs2 <- findMacs2()

##default..
##call peaks
proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "Clusters2", 
          pathToMacs2 = pathToMacs2, force = T)
##add peak matrix
proj2 <- addPeakMatrix(proj2)
getAvailableMatrices(proj2)
##get peak info
peaks.gr <- getPeakSet(proj2)
df_gr2 = data.frame(peaks.gr)
write.table(df_gr2, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.peak_info.m2_q0.01.txt', row.names = T, sep="\t", quote = FALSE)

##Identifying Marker peaks between cell classes
#Our scRNA labels
table(proj2$Clusters2)
#get peaks specific to cell class
markersPeaks <- getMarkerFeatures(ArchRProj = proj2, useMatrix = "PeakMatrix", groupBy = "Clusters2",
              bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
##all markers
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
write.csv(markerList, file='marker_list.peaks_per_cell_class.csv')
##graph marker peaks
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, 
                              cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmap_cell_class_marker_peaks.pdf", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)
##get tracks for marker genes and then print
p <- plotBrowserTrack(ArchRProj = proj2, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000)
plotPDF(plotList = p, name = "plot_tracks_for_marker_genes_cell_class.pdf", 
        ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)


##Motif and Feature Enrichment
##add motif
proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")
##look for motif enrichment in all cell classes marker peaks
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = proj2,
                                   peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "heatmap_enriched_motifs_markers.pdf", width = 12, height = 9, ArchRProj = proj2, addDOC = FALSE)


##ChromVAR Deviatons Enrichment -- different way of looking at motif enrichment
##load peak annotation if not there
if("Motif" %ni% names(proj2@peakAnnotation)){
  proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")
}
##make background peaks
proj2 <- addBgdPeaks(proj2)
##compute per-cell deviations
proj2 <- addDeviationsMatrix(ArchRProj = proj2, 
                             peakAnnotation = "Motif", force = TRUE)
##get deviations and plot
plotVarDev <- getVarDeviations(proj2, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "variable_motif_deviation_scores.pdf", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)
## extract a subset of motifs for downstream analysis (most variable)
motifs <- c("TAL", "CRX", "PITX", "GSC", "NEUROD1")
markerMotifs <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
##get just z scores also remove motifs that shouldn't be included
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs
##get impute weights from before to smooth signal
p <- plotGroups(ArchRProj = proj2, groupBy = "Clusters2", colorBy = "MotifMatrix", 
                name = markerMotifs, imputeWeights = getImputeWeights(proj2))
##plot and then save
p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
plotPDF(p, name = "groups_deviations_imputation.pdf", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)
##overlay zscores onto umap
p <- plotEmbedding(ArchRProj = proj2, colorBy = "MotifMatrix", 
                   name = sort(markerMotifs), embedding = "UMAP", imputeWeights = getImputeWeights(proj2))
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(p, name = "motif_umaps_zscores.pdf", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)


##compare with inferred gene expression umap
markerRNA <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA
p <- plotEmbedding(ArchRProj = proj2, colorBy = "GeneScoreMatrix", name = sort(markerRNA), 
                   embedding = "UMAP",imputeWeights = getImputeWeights(proj2))
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(p, name = "motif_umaps_inferred_gene_expression.pdf", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)

##compare with gene expression from rnaseq umap
markerRNA <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA
p <- plotEmbedding(ArchRProj = proj2, colorBy = "GeneIntegrationMatrix", name = sort(markerRNA), 
                   embedding = "UMAP",imputeWeights = getImputeWeights(proj2))
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(p, name = "motif_umaps_rnaseq_gene_expression.pdf", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)


##motif footprinting
##get motif postion
motifPositions <- getPositions(proj2)
motifPositions
##get motifs we're interested in
motifs <- c("TAL", "CRX", "PITX", "GSC", "NEUROD1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
##pseudo bulk if needed
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters2")
##get footprints
seFoot <- getFootprints(ArchRProj = proj2, positions = motifPositions[markerMotifs], 
                        groupBy = "Clusters2")
##Subtracting the Tn5 Bias
plotFootprints(seFoot = seFoot, ArchRProj = proj2, normMethod = "Subtract",
               plotName = "footprints_subtract_bias.pdf", addDOC = FALSE, smoothWindow = 5)
##Dividing by the Tn5 Bias
plotFootprints(seFoot = seFoot, ArchRProj = proj2, normMethod = "Divide",
               plotName = "footprints_divide_bias.pdf", addDOC = FALSE, smoothWindow = 5)
##Footprinting Without Normalization for Tn5 Bias
plotFootprints(seFoot = seFoot, ArchRProj = proj2, normMethod = "None", 
               plotName = "footprints_no_norm.pdf", addDOC = FALSE, smoothWindow = 5)
##TSS insertion profile
seTSS <- getFootprints(ArchRProj = proj2, positions = GRangesList(TSS = getTSS(proj2)), 
                       groupBy = "Clusters2", flank = 2000)
plotFootprints(seFoot = seTSS, ArchRProj = proj2, normMethod = "None",
               plotName = "TSS_insertion_no_norm.pdf", addDOC = FALSE, flank = 2000, flankNorm = 100)


##CoAccessibility
proj2 <- addCoAccessibility(ArchRProj = proj2, reducedDims = "Harmony")
##can change resolution to get more/less interactions
cA <- getCoAccessibility(ArchRProj = proj2, corCutOff = 0.5,
                         resolution = 10000, returnLoops = TRUE)
cA[[1]]
##plot for marker genes
markerGenes  <- c(
  "RHO", "SAG",  #Rods
  "ARR3", "GNAT2", #Cones
  "LHX1", "ONECUT2", #Horizontals
  "VSX2", "LHX4", #Bipolars
  "TFAP2A", "TFAP2B", #Amacrines
  "RBPMS", "ELAVL4", "GAP43", "POU4F1", "POU4F2", "THY1", #Ganglions
  "SOX9", "SLC1A3", #Mullers
  "GFAP", "S100A1" #Astrocytes
)
p <- plotBrowserTrack(ArchRProj = proj2, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000,
                      downstream = 50000,loops = getCoAccessibility(proj2))
plotPDF(plotList = p, name = "plot_tracks_marker_genes_with_CoAccessibility.pdf", 
        ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

##peak2gene

##add p2g links and retreive, hange resolution and returnloops depending on purpose
proj2 <- addPeak2GeneLinks(ArchRProj = proj2,reducedDims = "Harmony")
p2g <- getPeak2GeneLinks(ArchRProj = proj2, corCutOff = 0.45,
                         resolution = 10000, returnLoops = TRUE)
p2g
##plot tracks (need markerGenes in env)
p <- plotBrowserTrack(ArchRProj = proj2, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000,
                      loops = getPeak2GeneLinks(proj2))
plotPDF(plotList = p, name = "plot_tracks_marker_genes_Peak2Gene.pdf", 
        ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##plot heatmaps (plotting isn't working)
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Clusters2")
p
plotPDF(plotList = p, name = "plot_heatmap_marker_genes_Peak2Gene.pdf", 
        ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)


##save/load archr project
saveArchRProject(ArchRProj = proj2)
proj2 = loadArchRProject(path = 'human')













####################



##get co-accessabilty dataget peak info
CoAc <- getCoAccessibility(ArchRProj = proj2, corCutOff = 0.4,
                             resolution = 10000, returnLoops = FALSE)
write.table(CoAc, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.co_acc.m2_q0.01.cor4.txt', row.names = F, sep="\t", quote = FALSE)
##get peak info
peaks.gr <- getPeakSet(proj2)
df_gr2 = data.frame(peaks.gr)
write.table(df_gr2, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.peak_info.m2_q0.01.txt', row.names = T, sep="\t", quote = FALSE)


##get peak to gene data and rna index info
p2g <- getPeak2GeneLinks(ArchRProj = proj2, corCutOff = 0.4,
                         resolution = 10000, returnLoops = F)
write.table(p2g, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.p2g.m2_q0.01.cor4.txt', row.names = F, sep="\t", quote = FALSE)
rnaidx.info = metadata(p2g)[[2]]
rnaidx.df = data.frame(rnaidx.info)
write.table(rnaidx.df, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.m2_q0.01.rnaseq_info.txt', row.names = T, sep="\t", quote = FALSE)


saveArchRProject(ArchRProj = proj2)









