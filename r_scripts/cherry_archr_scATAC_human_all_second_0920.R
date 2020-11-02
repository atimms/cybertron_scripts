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


##get archr and set up env
#library(ggplot2)
#library(GenomicRanges)
library(chromVARmotifs)
library(pheatmap)
library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")
##Set working directory and specific dir in active RSS where the data is sitting
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all')


##input frag files
inputFiles = c('Hu5' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu5/outs/fragments.tsv.gz',
               'Hu7' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu7/outs/fragments.tsv.gz',
               'Hu8' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu8/outs/fragments.tsv.gz',
               'd53' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/d53/outs/fragments.tsv.gz',
               'd59' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/d59/outs/fragments.tsv.gz',
               'd74' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/d74/outs/fragments.tsv.gz',
               'd78' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/d78/outs/fragments.tsv.gz',
               'd113' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/d113/outs/fragments.tsv.gz',
               'd132' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/d132/outs/fragments.tsv.gz')
inputFiles
##make arrow file
ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = names(inputFiles),
                               filterTSS = 4, #Dont set this too high because you can always increase later
                               filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE)
ArrowFiles
##doublet inferrance
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                               knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                               LSIMethod = 1)
##Creating An ArchRProject
proj1 <- ArchRProject(ArrowFiles = ArrowFiles, 
  outputDirectory = "cluster_assignments",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
##info about project
proj1
paste0("Memory Size = ", round(object.size(proj1) / 10^6, 3), " MB") #memory used
getAvailableMatrices(proj1)
##plot sample stats comparing samples
p1 <- plotGroups(ArchRProj = proj1, groupBy = "Sample", colorBy = "cellColData", 
                 name = "TSSEnrichment", plotAs = "ridges")
p2 <- plotGroups(ArchRProj = proj1, groupBy = "Sample", colorBy = "cellColData", 
                 name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p3 <- plotGroups(ArchRProj = proj1, groupBy = "Sample", colorBy = "cellColData", 
                 name = "log10(nFrags)", plotAs = "ridges")
p4 <- plotGroups(ArchRProj = proj1, groupBy = "Sample", colorBy = "cellColData", 
                 name = "log10(nFrags)", plotAs = "violin",alpha = 0.4, addBoxPlot = TRUE)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj1, addDOC = FALSE, width = 4, height = 4)
p1 <- plotFragmentSizes(ArchRProj = proj1)
p2 <- plotTSSEnrichment(ArchRProj = proj1)
plotPDF(p1,p2, name = "QC_sample_fragSizes_TSSProfile.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##save/load archr project
saveArchRProject(ArchRProj = proj1)
#proj1 = loadArchRProject(path = 'cluster_assignments')

##filter doublets, can adjust to filter more is needed
proj1 <- filterDoublets(ArchRProj = proj1)
##Dimensionality Reduction and Clustering
##lsi if samples are pretty similiar
##The most common parameters to tweak are iterations, varFeatures, and resolution
##change sample cell from 10k to 50k?
proj1 <- addIterativeLSI(ArchRProj = proj1, useMatrix = "TileMatrix", name = "IterativeLSI", 
                         iterations = 4, clusterParams = list( #See Seurat::FindClusters
                           resolution = c(0.1,0.2,0.4), sampleCells = 50000, n.start = 10), 
                         varFeatures = 25000, dimsToUse = 1:30)
##use harmony if sample need more adjustment
proj1 <- addHarmony(ArchRProj = proj1, reducedDims = "IterativeLSI",
                    name = "Harmony", groupBy = "Sample")
##find cluster using seurat (can use LSI or harmony) -- force is needed is repeating
proj1 <- addClusters(input = proj1, reducedDims = "Harmony",
                     method = "Seurat", name = "Clusters", resolution = 0.2, force = TRUE)
counts_cluster_sample = table(proj1$Clusters, proj1$Sample)
write.csv(counts_cluster_sample, file='counts_cluster_sample.harmony.res0.2.csv')
##graph
cM <- confusionMatrix(paste0(proj1$Clusters), paste0(proj1$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
                        border_color = "black")
p
plotPDF(plotList = p, name = "cluster_sample_heatmap.harmony.res0.2.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##umap using harmony
proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "Harmony", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "UMAP_sample_clusters.harmony.res0.2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##finding marker genes
markersGS <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "GeneScoreMatrix", 
                               groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6
write.csv(markerList, file='marker_list.harmony.res0.2.csv')

##view markers with imputation and plot on umap
proj1 <- addImputeWeights(proj1)
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

p <- plotEmbedding(
  ArchRProj = proj1, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj1)
)
#Rearrange for grid plotting
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
plotPDF(plotList = p, name = "plot_marker_genes_with_imputation.harmony.res0.2.pdf", 
        ArchRProj = proj1,  addDOC = FALSE, width = 5, height = 5)
##get tracks for marker genes and then print
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000)
plotPDF(plotList = p, name = "plot_tracks_for_marker_genes.harmony.res0.2.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##integrate scRNASeq
##get seurat object
seRNA =  readRDS('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/seurat_analysis/all/human_harmony_clusters_defined.rds')
seRNA
##Unconstrained Integration
proj1 <- addGeneIntegrationMatrix(ArchRProj = proj1, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix",
                                  reducedDims = "Harmony", seRNA = seRNA, addToArrow = FALSE, groupRNA = "celltype",
                                  nameCell = "predictedCell", nameGroup = "predictedGroup", nameScore = "predictedScore")
##make color palatte and then plot
pal <- paletteDiscrete(values = seRNA$celltype)
p1 <- plotEmbedding(proj1, colorBy = "cellColData", name = "predictedGroup", pal = pal, rastr = FALSE)
plotPDF(p1, name = "UMAP_clusters_rnaseq.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##looks good so repeat but change addToArrow = T and force = T
addArchRThreads(threads = 1) ##only seems to work with single thread 
proj1 <- addGeneIntegrationMatrix(ArchRProj = proj1, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix",
                                  reducedDims = "Harmony", seRNA = seRNA, addToArrow = TRUE, force = TRUE, groupRNA = "celltype",
                                  nameCell = "predictedCell", nameGroup = "predictedGroup", nameScore = "predictedScore")
addArchRThreads(threads = 10) ##switch back to 10 threads
getAvailableMatrices(proj1) ##check we have GeneIntegrationMatrix

##add cluster2 to names the cluster
cM <- confusionMatrix(proj1$Clusters, proj1$cellNames)
labelOld <- rownames(cM)
##check the weird order and replicate
labelOld
labelNew <- c('Early progenitors', 'Ganglion precursors', 'Developing ganglions', 'AC/HC/GC precursors', 'Late progenitors', 'Developing bipolars', 'Mature rods', 'Bipolar/photoreceptor precursors', 'Rod precursors', 'Developing cones', 'Mature horizontals', 'Mature amacrines', 'Mature bipolars', 'Mature bipolars', 'Mature cones', 'Mature mullers', 'Developing amacrines', 'Horizontal precursors', 'Mature ganglions', 'Developing ganglions')
proj1$Clusters2 <- mapLabels(proj1$Clusters, newLabels = labelNew, oldLabels = labelOld)
p1 <- plotEmbedding(proj1, colorBy = "cellColData", name = "Clusters2", rastr = FALSE)
plotPDF(p1, name = "UMAP_cluster_names.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj1)
#proj1 = loadArchRProject(path = 'cluster_assignments')

##call peaks using min cells = 20 and q=0.01 for macs2
##change p2g and coac to 300k

##make pseudo-bulk replicates - change minCells to 20 (for mature ganglions)
proj2 <- addGroupCoverages(ArchRProj = proj1, groupBy = "Clusters2", minCells = 20)
#proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters2", minCells = 20, force = T)
saveArchRProject(ArchRProj = proj2, outputDirectory = "macs2_q0.01")
#proj2 = loadArchRProject(path = 'macs2_q0.01')
##call peaks with macs2
pathToMacs2 <- findMacs2()
proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "Clusters2", 
        pathToMacs2 = pathToMacs2, cutOff = 0.01)
##get GRanges object
getPeakSet(proj2)
##add peak matrix
proj2 <- addPeakMatrix(proj2)
getAvailableMatrices(proj2)
##Identifying Marker peaks between cell classes
#Our scRNA labels
table(proj2$Clusters2)
#get peaks specific to cell class
markersPeaks <- getMarkerFeatures(ArchRProj = proj2, 
                                  useMatrix = "PeakMatrix", groupBy = "Clusters2",
                                  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
##all markers
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
write.csv(markerList, file='marker_list.peaks_per_cell_class.macs2_q0.01.csv')
##graph marker peaks
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, 
                              cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmap_cell_class_marker_peaks.macs2_q0.01.pdf", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)
##get tracks for marker genes and then print
p <- plotBrowserTrack(ArchRProj = proj2, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000)
plotPDF(plotList = p, name = "plot_tracks_for_marker_genes_cell_class.macs2_q0.01.pdf", 
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
proj2 <- addCoAccessibility(ArchRProj = proj2, reducedDims = "Harmony", maxDist = 300000)
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

##add p2g links and retreive, change resolution and returnloops depending on purpose
##changed max dist to 300k
proj2 <- addPeak2GeneLinks(ArchRProj = proj2,reducedDims = "Harmony", maxDist = 300000)
p2g <- getPeak2GeneLinks(ArchRProj = proj2, corCutOff = 0.45,
                         resolution = 10000, returnLoops = TRUE)
p2g
##plot tracks (need markerGenes in env)
p <- plotBrowserTrack(ArchRProj = proj2, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000,
                      loops = getPeak2GeneLinks(proj2))
plotPDF(plotList = p, name = "plot_tracks_marker_genes_Peak2Gene.pdf", 
        ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##plot heatmaps (plotting isn't working?)
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Clusters2")
p
plotPDF(plotList = p, name = "plot_heatmap_marker_genes_Peak2Gene.pdf", 
        ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

##get data for coacc/p2g project

##get co-accessabilty data and peak info
CoAc <- getCoAccessibility(ArchRProj = proj2, corCutOff = 0.4,
                           resolution = 10000, returnLoops = FALSE)
write.table(CoAc, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.co_acc.m2_q0.01.cor4.txt', row.names = F, sep="\t", quote = FALSE)
peaks.gr <- getPeakSet(proj2)
df.peaks.gr = data.frame(peaks.gr)
write.table(df.peaks.gr, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.peak_info.m2_q0.01.txt', row.names = T, sep="\t", quote = FALSE)
##get peak to gene data and rna index info
p2g <- getPeak2GeneLinks(ArchRProj = proj2, corCutOff = 0.4,
                         resolution = 10000, returnLoops = F)
write.table(p2g, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.p2g.m2_q0.01.cor4.txt', row.names = F, sep="\t", quote = FALSE)
rnaidx.info = metadata(p2g)[[2]]
rnaidx.df = data.frame(rnaidx.info)
write.table(rnaidx.df, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.m2_q0.01.rnaseq_info.txt', row.names = T, sep="\t", quote = FALSE)


##trajectory analysis

##looking at the rods during development
rod_trajectory <- c("C19", "C20", "C7", "C9", "C8")
rod_trajectory
proj2 <- addTrajectory(ArchRProj = proj2, name = "rodU", 
                       groupBy = "Clusters",trajectory = rod_trajectory, 
                       embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj2$rodU[!is.na(proj2$rodU)])
p <- plotTrajectory(proj2, trajectory = "rodU", colorBy = "cellColData", name = "rodU")
plotPDF(p, name = "plot_rod_traj_UMAP.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with RHO inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "RHO", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "RHO", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "rod_traj_RHO_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with SAG inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "SAG", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "SAG", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "rod_traj_SAG_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##rod pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj2, name = "rodU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj2, name = "rodU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj2, name = "rodU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj2, name = "rodU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "rod_traj_comparison_heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE, width = 6, height = 8)

##looking at the Mullers during development
muller_trajectory <- c("C19", "C20", "C10")
muller_trajectory
proj2 <- addTrajectory(ArchRProj = proj2, name = "mullerU", 
                       groupBy = "Clusters",trajectory = muller_trajectory, 
                       embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj2$mullerU[!is.na(proj2$mullerU)])
p <- plotTrajectory(proj2, trajectory = "mullerU", colorBy = "cellColData", name = "mullerU")
plotPDF(p, name = "plot_muller_traj_UMAP.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with SOX9 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "SOX9", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "SOX9", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "muller_traj_SOX9_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with SLC1A3 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "SLC1A3", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "SLC1A3", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "muller_traj_SLC1A3_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##rod pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj2, name = "mullerU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj2, name = "mullerU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj2, name = "mullerU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj2, name = "mullerU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "muller_traj_comparison_heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE, width = 6, height = 8)

##looking at the ganglions during development
ganglion_trajectory <- c("C19", "C20", "C3", "C17", "C18")
ganglion_trajectory
proj2 <- addTrajectory(ArchRProj = proj2, name = "ganglionU", 
                       groupBy = "Clusters",trajectory = ganglion_trajectory, 
                       embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj2$ganglionU[!is.na(proj2$ganglionU)])
p <- plotTrajectory(proj2, trajectory = "ganglionU", colorBy = "cellColData", name = "ganglionU")
plotPDF(p, name = "plot_ganglion_traj_UMAP.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with POU4F1 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "ganglionU", colorBy = "GeneScoreMatrix", name = "POU4F1", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "ganglionU", colorBy = "GeneIntegrationMatrix", name = "POU4F1", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "ganglion_traj_POU4F1_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with THY1 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "ganglionU", colorBy = "GeneScoreMatrix", name = "THY1", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "ganglionU", colorBy = "GeneIntegrationMatrix", name = "THY1", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "ganglion_traj_THY1_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##ganglion pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj2, name = "ganglionU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj2, name = "ganglionU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj2, name = "ganglionU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj2, name = "ganglionU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "ganglion_traj_comparison_heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE, width = 6, height = 8)

##looking at the amacrines during development
amacrine_trajectory <- c("C19", "C20", "C1", "C16", "C14")
amacrine_trajectory
proj2 <- addTrajectory(ArchRProj = proj2, name = "amacrineU", 
                       groupBy = "Clusters",trajectory = amacrine_trajectory, 
                       embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj2$amacrineU[!is.na(proj2$amacrineU)])
p <- plotTrajectory(proj2, trajectory = "amacrineU", colorBy = "cellColData", name = "amacrineU")
plotPDF(p, name = "plot_amacrine_traj_UMAP.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with TFAP2A inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "TFAP2A", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "TFAP2A", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "amacrine_traj_TFAP2A_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with TFAP2B inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "TFAP2B", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "TFAP2B", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "amacrine_traj_TFAP2B_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##amacrine pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj2, name = "amacrineU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj2, name = "amacrineU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj2, name = "amacrineU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj2, name = "amacrineU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "amacrine_traj_comparison_heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE, width = 6, height = 8)

##looking at the horizontals during development
horizontal_trajectory <- c("C19", "C20", "C1", "C15", "C13")
horizontal_trajectory
proj2 <- addTrajectory(ArchRProj = proj2, name = "horizontalU", 
                       groupBy = "Clusters",trajectory = horizontal_trajectory, 
                       embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj2$horizontalU[!is.na(proj2$horizontalU)])
p <- plotTrajectory(proj2, trajectory = "horizontalU", colorBy = "cellColData", name = "horizontalU")
plotPDF(p, name = "plot_horizontal_traj_UMAP.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with LHX1 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "LHX1", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "LHX1", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "horizontal_traj_LHX1_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with ONECUT2 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "horizontal_traj_ONECUT2_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##horizontal pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj2, name = "horizontalU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj2, name = "horizontalU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj2, name = "horizontalU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj2, name = "horizontalU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "horizontal_traj_comparison_heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE, width = 6, height = 8)


##looking at the bipolars during development
bipolar_trajectory <- c("C19", "C20", "C7", "C2", "C11", "C12")
bipolar_trajectory
proj2 <- addTrajectory(ArchRProj = proj2, name = "bipolarU", 
                       groupBy = "Clusters",trajectory = bipolar_trajectory, 
                       embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj2$bipolarU[!is.na(proj2$bipolarU)])
p <- plotTrajectory(proj2, trajectory = "bipolarU", colorBy = "cellColData", name = "bipolarU")
plotPDF(p, name = "plot_bipolar_traj_UMAP.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with VSX2 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "VSX2", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "VSX2", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "bipolar_traj_VSX2_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with LHX4 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "bipolarU", colorBy = "GeneScoreMatrix", name = "LHX4", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "bipolarU", colorBy = "GeneIntegrationMatrix", name = "LHX4", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "bipolar_traj_LHX4_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##bipolar pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj2, name = "bipolarU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj2, name = "bipolarU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj2, name = "bipolarU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj2, name = "bipolarU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "bipolar_traj_comparison_heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE, width = 6, height = 8)


##looking at the cones during development
cone_trajectory <- c("C19", "C20", "C7", "C5", "C6")
cone_trajectory
proj2 <- addTrajectory(ArchRProj = proj2, name = "coneU", 
                       groupBy = "Clusters",trajectory = cone_trajectory, 
                       embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj2$coneU[!is.na(proj2$coneU)])
p <- plotTrajectory(proj2, trajectory = "coneU", colorBy = "cellColData", name = "coneU")
plotPDF(p, name = "plot_cone_traj_UMAP.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with ARR3 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ARR3", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ARR3", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "cone_traj_ARR3_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with GNAT2 inferred and rnaseq expression
p1 <- plotTrajectory(proj2, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "GNAT2", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj2, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "GNAT2", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
plotPDF(p1,p2, name = "cone_traj_GNAT2_expression.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##cone pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj2, name = "coneU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj2, name = "coneU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj2, name = "coneU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj2, name = "coneU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "cone_traj_comparison_heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE, width = 6, height = 8)


##save/load archr project
saveArchRProject(ArchRProj = proj2)
proj2 = loadArchRProject(path = 'cluster_assignments')

##get cell class info on each cell class
ids = rownames(getCellColData(proj2))
cell_class = getCellColData(proj2)$Clusters2
cell.info = data.frame(ids,cell_class)
write.csv(cell.info, file='human_all.cell_class_info.092420.csv')

