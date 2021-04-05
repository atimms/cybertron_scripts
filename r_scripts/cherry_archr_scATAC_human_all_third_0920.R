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
library(ComplexHeatmap)
library(circlize)
set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")
##Set working directory and specific dir in active RSS where the data is sitting
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all_third')


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
                      outputDirectory = "human_all",
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
proj1 = loadArchRProject(path = 'human_all')

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
labelNew <- c('Early Progenitors', 'Ganglion Precursors', 'Developing Ganglions', 'AC/HC/GC Precursors', 'Late Progenitors', 'Mature Rods', 'Developing Rods', 'Photoreceptor/Bipolar Precursors', 'Developing Cones', 'Mature Horizontals', 'Mature Amacrines', 'Mature Bipolars', 'Mature Bipolars', 'Mature Mullers', 'Mature Cones', 'Mature Ganglions', 'Developing Amacrines', 'Developing Horizontals', 'Developing Bipolars', 'Developing Ganglions')
proj1$Clusters2 <- mapLabels(proj1$Clusters, newLabels = labelNew, oldLabels = labelOld)
p1 <- plotEmbedding(proj1, colorBy = "cellColData", name = "Clusters2", rastr = FALSE)
plotPDF(p1, name = "UMAP_cluster_names.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

#saveArchRProject(ArchRProj = proj1)
#proj1 = loadArchRProject(path = 'human_all')

##get cell class info on each cell class
ids = rownames(getCellColData(proj1))
cell_class = getCellColData(proj1)$Clusters2
cell.info = data.frame(ids,cell_class)
write.csv(cell.info, file='human_all.cell_class_info.092520.csv')

##I copied the folder to all_third_backup_092520, so if need to go back change name to all_third and pick up from here

##call peaks using min cells = 20 and q=0.01 for macs2
##change p2g and coac to 300k

##make pseudo-bulk replicates - change minCells to 20 (for mature ganglions)
proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "Clusters2", minCells = 20)
##call peaks with macs2
pathToMacs2 <- findMacs2()
proj1 <- addReproduciblePeakSet(ArchRProj = proj1, groupBy = "Clusters2", 
                                pathToMacs2 = pathToMacs2, cutOff = 0.01)
##get GRanges object
getPeakSet(proj1)
##add peak matrix
proj1 <- addPeakMatrix(proj1)
getAvailableMatrices(proj1)
##Identifying Marker peaks between cell classes
#Our scRNA labels
table(proj1$Clusters2)
#get peaks specific to cell class
markersPeaks <- getMarkerFeatures(ArchRProj = proj1, 
                                  useMatrix = "PeakMatrix", groupBy = "Clusters2",
                                  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
##all markers
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
write.csv(markerList, file='marker_list.peaks_per_cell_class.macs2_q0.01.csv')
##graph marker peaks
heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, 
                  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmap_cell_class_marker_peaks.macs2_q0.01.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)
##instead of drawing plot make a matrix
heatmap.matrix <- plotMarkerHeatmap(seMarker = markersPeaks, 
        cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE, returnMat = TRUE)
##get csv file from matrix
write.csv(heatmap.matrix, file='heatmap_cell_class_marker_peaks.matrix.csv')


##get tracks for marker genes and then print
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000)
plotPDF(plotList = p, name = "plot_tracks_for_marker_genes_cell_class.macs2_q0.01.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj1)
proj1 = loadArchRProject(path = 'human_all')


##Motif and Feature Enrichment
##add motif
proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")
##look for motif enrichment in all cell classes marker peaks
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = proj1,
                                   peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "heatmap_enriched_motifs_markers.pdf", width = 12, height = 9, ArchRProj = proj1, addDOC = FALSE)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 14, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "heatmap_enriched_motifs_markers.max14.pdf", width = 12, height = 9, ArchRProj = proj1, addDOC = FALSE)
##get matrix instead of map
motif_heatmap_matrix <- plotEnrichHeatmap(enrichMotifs, n = 14, transpose = TRUE, returnMatrix = T)
write.csv(motif_heatmap_matrix, file='motif_heatmap_matrix.max14.csv')
##modify file manually and then re draw
heatmapEM_new = read.csv('motif_heatmap_matrix.modified.csv', header=T, row.names=1)
heatmapEM_new_matrix = data.matrix(heatmapEM_new)
#pheatmap(heatmapEM_new, cluster_cols = F, cluster_rows = F)
color_for_hm = colorRamp2(seq(min(heatmapEM_new_matrix), max(heatmapEM_new_matrix), length = 3), c("grey", "blue", "black"))
Heatmap(heatmapEM_new_matrix, cluster_rows = F, cluster_columns = F, col = color_for_hm,
        name = "Norm. Enrichment −log10(P−adj) [0−Max]")
dev.copy2pdf(file="heatmap_enriched_motifs_markers.modified.pdf", width = 20)

##tim modified the order, so we can redraw
heatmapEM_new = read.csv('motif_heatmap_matrix.modified_TC.csv', header=T, row.names=1)
heatmapEM_new_matrix = data.matrix(heatmapEM_new)
#pheatmap(heatmapEM_new, cluster_cols = F, cluster_rows = F)
color_for_hm = colorRamp2(seq(min(heatmapEM_new_matrix), max(heatmapEM_new_matrix), length = 3), c("grey", "blue", "black"))
Heatmap(heatmapEM_new_matrix, cluster_rows = F, cluster_columns = F, col = color_for_hm,
        name = "Norm. Enrichment −log10(P−adj) [0−Max]")
dev.copy2pdf(file="heatmap_enriched_motifs_markers.modified_TC.pdf", width = 20)
dev.copy2pdf(file="heatmap_enriched_motifs_markers.modified_TC2.pdf", width = 16, height = 8)



##ChromVAR Deviatons Enrichment -- different way of looking at motif enrichment
##load peak annotation if not there
if("Motif" %ni% names(proj1@peakAnnotation)){
  proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")
}
##make background peaks
proj1 <- addBgdPeaks(proj1)
##compute per-cell deviations
proj1 <- addDeviationsMatrix(ArchRProj = proj1, 
                             peakAnnotation = "Motif", force = TRUE)
##get deviations and plot
plotVarDev <- getVarDeviations(proj1, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "variable_motif_deviation_scores.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
## extract a subset of motifs for downstream analysis (most variable)
motifs <- c("TAL", "CRX", "PITX", "GSC", "NEUROD1")
markerMotifs <- getFeatures(proj1, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
##get just z scores also remove motifs that shouldn't be included
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs
##get impute weights from before to smooth signal
p <- plotGroups(ArchRProj = proj1, groupBy = "Clusters2", colorBy = "MotifMatrix", 
                name = markerMotifs, imputeWeights = getImputeWeights(proj1))
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
plotPDF(p, name = "groups_deviations_imputation.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
##overlay zscores onto umap
p <- plotEmbedding(ArchRProj = proj1, colorBy = "MotifMatrix", 
                   name = sort(markerMotifs), embedding = "UMAP", imputeWeights = getImputeWeights(proj1))
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
plotPDF(p, name = "motif_umaps_zscores.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)


##compare with inferred gene expression umap
markerRNA <- getFeatures(proj1, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA
p <- plotEmbedding(ArchRProj = proj1, colorBy = "GeneScoreMatrix", name = sort(markerRNA), 
                   embedding = "UMAP",imputeWeights = getImputeWeights(proj1))
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
plotPDF(p, name = "motif_umaps_inferred_gene_expression.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)

##compare with gene expression from rnaseq umap
markerRNA <- getFeatures(proj1, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA
p <- plotEmbedding(ArchRProj = proj1, colorBy = "GeneIntegrationMatrix", name = sort(markerRNA), 
                   embedding = "UMAP",imputeWeights = getImputeWeights(proj1))
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
plotPDF(p, name = "motif_umaps_rnaseq_gene_expression.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)


##motif footprinting
##get motif postion
motifPositions <- getPositions(proj1)
motifPositions
##get motifs we're interested in
motifs <- c("TAL", "CRX", "PITX", "GSC", "NEUROD1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
##pseudo bulk if needed
proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "Clusters2")
##get footprints
seFoot <- getFootprints(ArchRProj = proj1, positions = motifPositions[markerMotifs], 
                        groupBy = "Clusters2")
##Subtracting the Tn5 Bias
plotFootprints(seFoot = seFoot, ArchRProj = proj1, normMethod = "Subtract",
               plotName = "footprints_subtract_bias.pdf", addDOC = FALSE, smoothWindow = 5)
##Dividing by the Tn5 Bias
plotFootprints(seFoot = seFoot, ArchRProj = proj1, normMethod = "Divide",
               plotName = "footprints_divide_bias.pdf", addDOC = FALSE, smoothWindow = 5)
##Footprinting Without Normalization for Tn5 Bias
plotFootprints(seFoot = seFoot, ArchRProj = proj1, normMethod = "None", 
               plotName = "footprints_no_norm.pdf", addDOC = FALSE, smoothWindow = 5)
##TSS insertion profile
seTSS <- getFootprints(ArchRProj = proj1, positions = GRangesList(TSS = getTSS(proj1)), 
                       groupBy = "Clusters2", flank = 2000)
plotFootprints(seFoot = seTSS, ArchRProj = proj1, normMethod = "None",
               plotName = "TSS_insertion_no_norm.pdf", addDOC = FALSE, flank = 2000, flankNorm = 100)


##CoAccessibility
proj1 <- addCoAccessibility(ArchRProj = proj1, reducedDims = "Harmony", maxDist = 300000)
##can change resolution to get more/less interactions
cA <- getCoAccessibility(ArchRProj = proj1, corCutOff = 0.5,
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
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000,
                      downstream = 50000, loops = getCoAccessibility(proj1))
plotPDF(plotList = p, name = "plot_tracks_marker_genes_with_CoAccessibility.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##peak2gene

##add p2g links and retreive, change resolution and returnloops depending on purpose
##changed max dist to 300k
proj1 <- addPeak2GeneLinks(ArchRProj = proj1,reducedDims = "Harmony", maxDist = 300000)
p2g <- getPeak2GeneLinks(ArchRProj = proj1, corCutOff = 0.45,
                         resolution = 10000, returnLoops = TRUE)
p2g
##plot tracks (need markerGenes in env)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "plot_tracks_marker_genes_Peak2Gene.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##plot heatmaps (plotting isn't working?)
p <- plotPeak2GeneHeatmap(ArchRProj = proj1, groupBy = "Clusters2")
p
plotPDF(plotList = p, name = "plot_heatmap_marker_genes_Peak2Gene.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##get data for coacc/p2g project

##get co-accessabilty data and peak info
CoAc <- getCoAccessibility(ArchRProj = proj1, corCutOff = 0.4,
                           resolution = 10000, returnLoops = FALSE)
write.table(CoAc, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.co_acc.m2_q0.01.cor4.txt', row.names = F, sep="\t", quote = FALSE)
peaks.gr <- getPeakSet(proj1)
df.peaks.gr = data.frame(peaks.gr)
write.table(df.peaks.gr, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.peak_info.m2_q0.01.txt', row.names = T, sep="\t", quote = FALSE)
##get peak to gene data and rna index info
p2g <- getPeak2GeneLinks(ArchRProj = proj1, corCutOff = 0.4,
                         resolution = 10000, returnLoops = F)
write.table(p2g, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.p2g.m2_q0.01.cor4.txt', row.names = F, sep="\t", quote = FALSE)
rnaidx.info = metadata(p2g)[[2]]
rnaidx.df = data.frame(rnaidx.info)
write.table(rnaidx.df, '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920/human_all.min_cell20.m2_q0.01.rnaseq_info.txt', row.names = T, sep="\t", quote = FALSE)

##gene score matrix and gene integration UMAPS for LINC00461 
genesToPlot  <- c('LINC00461')
p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "GeneScoreMatrix", name = genesToPlot, 
                    embedding = "UMAP",imputeWeights = getImputeWeights(proj1))
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "GeneIntegrationMatrix", name = genesToPlot, 
                    embedding = "UMAP",imputeWeights = getImputeWeights(proj1), log2Norm = T)

plotPDF(p1, p2, name = "LINC00461.gene_expression.UMAP.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)




##trajectory analysis

##looking at the rods during development
rod_trajectory <- c("C20", "C6", "C2", "C1")
#rod_trajectory <- c("C20", "C2", "C1")
#rod_trajectory <- c("C19","C20","C6", "C2", "C1")
rod_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "rodU", 
                       groupBy = "Clusters",trajectory = rod_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$rodU[!is.na(proj1$rodU)])
p <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "cellColData", name = "rodU")
#plotPDF(p, name = "rod_trajectory_UMAP.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
plotPDF(p, name = "rod_trajectory_UMAP.C20_6_2_1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##rod pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "rodU", useMatrix = "MotifMatrix", log2Norm = FALSE)
#p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0.9, labelTop = 10000, returnMat = T)
p1_rownames <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), returnMat = T)
p1_rownames = rownames(p1_rownames)
write.csv(p1_rownames, file='rod_trajectory.C20_6_2_1.motif_rownames.csv')
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
p1
trajGSM <- getTrajectory(ArchRProj = proj1, name = "rodU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2_rownames <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"), returnMat = T)
p2_rownames = rownames(p2_rownames)
write.csv(p2_rownames, file='rod_trajectory.C20_6_2_1.genescore_rownames.csv')
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
p2
p2n <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"), grepExclude = '-AS1')
p2n
trajGIM <- getTrajectory(ArchRProj = proj1, name = "rodU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3_rownames <- plotTrajectoryHeatmap(trajGIM, pal = paletteContinuous(set = "blueYellow"), returnMat = T)
p3_rownames = rownames(p3_rownames)
write.csv(p3_rownames, file='rod_trajectory.C20_6_2_1.geneintegration_rownames.csv')
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "rodU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4_rownames <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), returnMat = T)
p4_rownames = rownames(p4_rownames)
write.csv(p4_rownames, file='rod_trajectory.C20_6_2_1.peak_rownames.csv')
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
#plotPDF(p1, p2, p3, p4, name = "rod_trajectory_comparison_heatmaps.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
plotPDF(p1, p2, p3, p4, name = "rod_trajectory_comparison_heatmaps.C20_6_2_1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##label geneintegration with choosen genes
markersToLabel = c('chr11:CCND1', 'chr13:FLT1', 'chr8:CLU', 'chr4:SFRP2', 'chr17:SDK2', 'chr13:FGF9', 'chr4:GRIA2', 'chr13:DCT', 'chr6:EPHA7', 'chr2:NRXN1', 'chr12:NEUROD4', 'chr15:CHRNB4', 'chr2:CERKL', 'chr8:RP1', 'chr3:SAMD7', 'chr5:PDE6A', 'chr2:SAG', 'chr8:RIMS2', 'chr3:RHO', 'chr7:BBS9', 'chr6:AHI1', 'chr4:CNGA1', 'chrX:RS1', 'chr10:HTRA1', 'chr4:PROM1', 'chr11:ROM1', 'chr4:BBS7', 'chr1:ABCA4', 'chr6:PRPH2', 'chr17:AIPL1', 'chr17:RCVRN', 'chr6:TULP1', 'chr10:PCDH15', 'chr1:USH2A', 'chr11:MYO7A')
p5 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"), labelMarkers = markersToLabel, labelTop = 0)
p5
plotPDF(p1, p2, p3, p4, p5, name = "rod_trajectory_comparison_heatmaps.take2.C20_6_2_1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)

##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "rod_trajectory.C20_6_2_1.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "rod_trajectory.C20_6_2_1.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the cones during development
cone_trajectory <- c("C20", "C6", "C7", "C8")
cone_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "coneU", 
                       groupBy = "Clusters",trajectory = cone_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$coneU[!is.na(proj1$coneU)])
p <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "cellColData", name = "coneU")
plotPDF(p, name = "plot_cone_traj_UMAP.C20_6_7_8.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##cone pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "coneU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "coneU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "coneU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
p3_rownames <- plotTrajectoryHeatmap(trajGIM, pal = paletteContinuous(set = "blueYellow"), returnMat = T)
p3_rownames = rownames(p3_rownames)
write.csv(p3_rownames, file='cone_trajectory.C20_6_7_8.geneintegration_rownames.csv')
trajPM  <- getTrajectory(ArchRProj = proj1, name = "coneU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "cone_trajectory_comparison_heatmaps.C20_6_7_8.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "cone_trajectory.C20_6_7_8.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "cone_trajectory.C20_6_7_8.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the horizontals during development
horizontal_trajectory <- c("C20", "C5", "C15", "C13")
horizontal_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "horizontalU", 
                       groupBy = "Clusters",trajectory = horizontal_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$horizontalU[!is.na(proj1$horizontalU)])
p <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "cellColData", name = "horizontalU")
plotPDF(p, name = "plot_horizontal_traj_UMAP.C20_5_15_13.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##horizontal pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "horizontalU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "horizontalU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "horizontalU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "horizontalU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "horizontal_trajectory_comparison_heatmaps.C20_5_15_13.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "horizontal_trajectory.C20_5_15_13.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "horizontal_trajectory.C20_5_15_13.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the cone_bipolars during development
cone_bipolar_trajectory <- c("C20", "C6", "C11", "C10")
cone_bipolar_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "cone_bipolarU", 
                       groupBy = "Clusters",trajectory = cone_bipolar_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$cone_bipolarU[!is.na(proj1$cone_bipolarU)])
p <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "cellColData", name = "cone_bipolarU")
plotPDF(p, name = "plot_cone_bipolar_traj_UMAP.C20_6_11_10.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##cone_bipolar pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "cone_bipolarU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "cone_bipolarU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "cone_bipolarU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "cone_bipolarU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "cone_bipolar_trajectory_comparison_heatmaps.C20_6_11_10.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "cone_bipolar_trajectory.C20_6_11_10.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "cone_bipolar_trajectory.C20_6_11_10.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the rod_bipolars during development
rod_bipolar_trajectory <- c("C20", "C6", "C11", "C9")
rod_bipolar_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "rod_bipolarU", 
                       groupBy = "Clusters",trajectory = rod_bipolar_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$rod_bipolarU[!is.na(proj1$rod_bipolarU)])
p <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "cellColData", name = "rod_bipolarU")
plotPDF(p, name = "plot_rod_bipolar_traj_UMAP.C20_6_11_9.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##rod_bipolar pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "rod_bipolarU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "rod_bipolarU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "rod_bipolarU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "rod_bipolarU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "rod_bipolar_trajectory_comparison_heatmaps.C20_6_11_9.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "rod_bipolar_trajectory.C20_6_11_9.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "rod_bipolar_trajectory.C20_6_11_9.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the mullers during development
muller_trajectory <- c("C19", "C20", "C12")
muller_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "mullerU", 
                       groupBy = "Clusters",trajectory = muller_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$mullerU[!is.na(proj1$mullerU)])
p <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "cellColData", name = "mullerU")
plotPDF(p, name = "plot_muller_traj_UMAP.C19_20_12.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##muller pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "muller_trajectory_comparison_heatmaps.C19_20_12.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "muller_trajectory.C19_20_12.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "muller_trajectory.C19_20_12.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the mullers during development - take 2
muller_trajectory <- c("C20", "C12")
muller_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "mullerU", 
                       groupBy = "Clusters",trajectory = muller_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$mullerU[!is.na(proj1$mullerU)])
p <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "cellColData", name = "mullerU")
plotPDF(p, name = "plot_muller_traj_UMAP.C20_12.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##muller pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "mullerU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "muller_trajectory_comparison_heatmaps.C20_12.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "muller_trajectory.C20_12.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "muller_trajectory.C20_12.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)



##looking at the amacrines during development
amacrine_trajectory <- c("C20", "C5", "C16", "C14")
amacrine_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "amacrineU", 
                       groupBy = "Clusters",trajectory = amacrine_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$amacrineU[!is.na(proj1$amacrineU)])
p <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "cellColData", name = "amacrineU")
plotPDF(p, name = "plot_amacrine_traj_UMAP.C20_5_16_14.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##amacrine pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "amacrineU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "amacrineU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "amacrineU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "amacrineU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "amacrine_trajectory_comparison_heatmaps.C20_5_16_14.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "amacrine_trajectory.C20_5_16_14.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "amacrine_trajectory.C20_5_16_14.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the ganglions during development
ganglion_trajectory <- c("C19", "C4", "C3")
ganglion_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "ganglionU", 
                       groupBy = "Clusters",trajectory = ganglion_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$ganglionU[!is.na(proj1$ganglionU)])
p <- plotTrajectory(proj1, trajectory = "ganglionU", colorBy = "cellColData", name = "ganglionU")
plotPDF(p, name = "plot_ganglion_traj_UMAP.C19_4_3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##ganglion pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "ganglionU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "ganglionU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "ganglionU", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "ganglionU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "ganglion_trajectory_comparison_heatmaps.C19_4_3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "ganglion_trajectory.C19_4_3.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "ganglion_trajectory.C19_4_3.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)

##looking at the ganglions during development with different clusters
ganglion_trajectory2 <- c("C19", "C4", "C3", "C17", "C18")
ganglion_trajectory2
proj1 <- addTrajectory(ArchRProj = proj1, name = "ganglion2U", 
                       groupBy = "Clusters",trajectory = ganglion_trajectory2, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$ganglion2U[!is.na(proj1$ganglion2U)])
p <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "cellColData", name = "ganglion2U")
plotPDF(p, name = "plot_ganglion_traj_UMAP.C19_4_3_17_18.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##ganglion pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "ganglion2U", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "ganglion2U", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = proj1, name = "ganglion2U", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "ganglion2U", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, p4, name = "ganglion_trajectory2_comparison_heatmaps.C19_4_3_17_18.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
ht1 + ht2
plotPDF(ht1, ht2, name = "ganglion_trajectory2.C19_4_3_17_18.GSM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "ganglion_trajectory2.C19_4_3_17_18.GIM_MM.1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 8, height = 12)


##save/load archr project
saveArchRProject(ArchRProj = proj1)
proj1 = loadArchRProject(path = 'human_all')


##make new plots 1026......
genelist_102620 = c('RHO', 'SAG', 'ARR3', 'GNAT2', 'RBPMS', 'ELAVL4', 'GAP43', 'THY1', 'SLC1A3', 'GFAP', 'S100A1', 'CRX', 'NRL', 'NR2E3', 'LHX1', 'ONECUT2', 'OTX2', 'VSX2', 'LHX4', 'SOX9', 'SOX2', 'PAX6', 'TFAP2B', 'TFAP2A', 'POU4F2', 'POU4F1', 'MEIS1', 'MEIS2', 'MEIS3', 'PKNOX1', 'PKNOX2', 'PBX1', 'PBX2', 'PBX3')

##genescore i.e. expession from atac
p <- plotEmbedding(ArchRProj = proj1, colorBy = "GeneScoreMatrix", name = genelist_102620, 
  embedding = "UMAP", imputeWeights = getImputeWeights(proj1))
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
plotPDF(plotList = p, name = "GeneScore.UMAP.genelist102620.pdf", 
        ArchRProj = proj1,  addDOC = FALSE, width = 5, height = 5)

##geneintegration i.e. expession from rnaseq
p <- plotEmbedding(ArchRProj = proj1, colorBy = "GeneIntegrationMatrix", name = genelist_102620, 
                   embedding = "UMAP", imputeWeights = getImputeWeights(proj1), log2Norm = T)
plotPDF(plotList = p, name = "GeneIntegration.UMAP.genelist102720.pdf", 
        ArchRProj = proj1,  addDOC = FALSE, width = 5, height = 5)

##get tracks for marker genes and then print
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = genelist_102620, upstream = 50000, downstream = 50000)
plotPDF(plotList = p, name = "tracks.genelist102620.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##motif deviation
##ChromVAR Deviatons Enrichment -- different way of looking at motif enrichment
##load peak annotation if not there
if("Motif" %ni% names(proj1@peakAnnotation)){
  proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")
}
##make background peaks
proj1 <- addBgdPeaks(proj1)
##compute per-cell deviations
proj1 <- addDeviationsMatrix(ArchRProj = proj1, 
                             peakAnnotation = "Motif", force = TRUE)



## extract a subset of motifs for downstream analysis 
markerMotifs <- getFeatures(proj1, select = paste(genelist_102620, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
##get just z scores also remove motifs that shouldn't be included
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs
##get impute weights from before to smooth signal
p <- plotGroups(ArchRProj = proj1, groupBy = "Clusters2", colorBy = "MotifMatrix", 
                name = markerMotifs, imputeWeights = getImputeWeights(proj1))
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
plotPDF(p, name = "motif_deviations_imputation.genelist102620.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
##overlay zscores onto umap
p <- plotEmbedding(ArchRProj = proj1, colorBy = "MotifMatrix", 
                   name = sort(markerMotifs), embedding = "UMAP", imputeWeights = getImputeWeights(proj1))
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
plotPDF(p, name = "motif_umaps_zscores.genelist102620.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)

##coaccesssability
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = genelist_102620, upstream = 50000,
                      downstream = 50000,loops = getCoAccessibility(proj1))
plotPDF(plotList = p, name = "CoAccessibility_tracks.genelist102620.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
##peak2gene
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = genelist_102620, upstream = 50000, downstream = 50000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.genelist102620.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = genelist_102620, upstream = 10000, downstream = 10000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.10k.genelist102620.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)

p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = genelist_102620, upstream = 10000, downstream = 10000,
                      loops = getPeak2GeneLinks(proj1, corCutOff = 0.7))
plotPDF(plotList = p, name = "Peak2Gene_tracks.10k_0.7corr.genelist102620.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
##gene specific
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('RHO'), upstream = 10000, downstream = 10000,
                      loops = getPeak2GeneLinks(proj1, corCutOff = 0.8))
plotPDF(plotList = p, name = "Peak2Gene_tracks.10k_0.8corr.RHO.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('ABCA4'), upstream = 170000, downstream = 70000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.70_170k.ABCA4.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('ABCA4'), upstream = 170000, downstream = 70000,
                      loops = getPeak2GeneLinks(proj1, corCutOff = 0.7))
plotPDF(plotList = p, name = "Peak2Gene_tracks.70_170k_0.7corr.ABCA4.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('ATOH7'), upstream = 10000, downstream = 10000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.10k.ATOH7.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('ATOH7'), upstream = 10000, downstream = 60000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.10_160k.ATOH7.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('PAX6'), upstream = 200000, downstream = 200000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.200k.PAX6.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('LINC00461'), upstream = 160000, downstream = 20000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.20_160k.LINC00461.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)

p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('HTRA1'), upstream = 100000, downstream = 100000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.100k.HTRA1.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)

p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('PLEKHA1'), upstream = 100000, downstream = 100000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.100k.PLEKHA1.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)

p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('ARMS2'), upstream = 20000, downstream = 20000,
                      loops = getPeak2GeneLinks(proj1))
plotPDF(plotList = p, name = "Peak2Gene_tracks.20k.ARMS2.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)

p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = c('CRX'), upstream = 50000, downstream = 50000,
                      loops = getPeak2GeneLinks(proj1, corCutOff = 0.8))
plotPDF(plotList = p, name = "Peak2Gene_tracks.50k_cor0.8.CRX.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 7, height = 5)

##get genescore matrix for each cell class
seGS <- getGroupSE(proj1, useMatrix = "GeneScoreMatrix", groupBy = "Clusters2")
seGS
cellClass_gene_matrix = assays(seGS)$GeneScoreMatrix
genenames = rowData(seGS)$name
write.csv(cellClass_gene_matrix, file='cell_class_gene_matrix.102820.csv')
write.csv(genenames, file='cell_class_gene_matrix_genes.102820.csv')

##gene based trajectory analysis
##save/load archr project
#saveArchRProject(ArchRProj = proj1)
proj1 = loadArchRProject(path = 'human_all')

##rod data
##looking at the rods during development
rod_trajectory <- c("C20", "C6", "C2", "C1")
rod_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "rodU", 
                       groupBy = "Clusters",trajectory = rod_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "RHO", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "RHO", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_RHO.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_CRX.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_OTX2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "MEF2D", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "MEF2D", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_MEF2D.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "MEF2C", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "MEF2C", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_MEF2C.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "NRL", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "NRL", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_NRL.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "NR2E3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "NR2E3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_NR2E3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "SAG", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "SAG", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_SAG.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "HES1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "HES1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_HES1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "NOTCH1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "NOTCH1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_NOTCH1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneScoreMatrix", name = "PDRM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rodU", colorBy = "GeneIntegrationMatrix", name = "PDRM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_trajectory.C20_6_2_1.PT_PDRM1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##looking at the cones during development
cone_trajectory <- c("C20", "C6", "C7", "C8")
cone_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "coneU", 
                       groupBy = "Clusters",trajectory = cone_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ARR3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ARR3", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_ARR3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "NRL", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "NRL", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_NRL.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_CRX.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_OTX2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "OPN1SW", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "OPN1SW", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_OPN1SW.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "OPN1LW", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "OPN1LW", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_OPN1LW.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ONECUT1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ONECUT1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_ONECUT1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_ONECUT2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "ONECUT3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "ONECUT3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_ONECUT3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneScoreMatrix", name = "PDRM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "coneU", colorBy = "GeneIntegrationMatrix", name = "PDRM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_trajectory.C20_6_7_8.PT_PDRM1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##looking at the horizontals during development
horizontal_trajectory <- c("C20", "C5", "C15", "C13")
horizontal_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "horizontalU", 
                       groupBy = "Clusters",trajectory = horizontal_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "LHX1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "LHX1", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_LHX1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "ONECUT1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "ONECUT1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_ONECUT1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_ONECUT2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "ONECUT3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "ONECUT3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_ONECUT3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "POU4F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "POU4F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_POU4F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "CALB1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "CALB1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_CALB1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "PTF1A", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "PTF1A", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_PTF1A.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "TFAP2B", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "TFAP2B", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_TFAP2B.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "horizontalU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "horizontal_trajectory.C20_5_15_13.PT_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##looking at the cone_bipolars during development
cone_bipolar_trajectory <- c("C20", "C6", "C11", "C10")
cone_bipolar_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "cone_bipolarU", 
                       groupBy = "Clusters",trajectory = cone_bipolar_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_OTX2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_CRX.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_PRDM1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NEUROD4", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NEUROD4", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_NEUROD4.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "GRM6", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "GRM6", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_GRM6.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "TRPM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "TRPM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_TRPM1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIA", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIA", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_NFIA.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIB", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIB", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_NFIB.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIX", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIX", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_NFIX.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneScoreMatrix", name = "VSX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "cone_bipolarU", colorBy = "GeneIntegrationMatrix", name = "VSX2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "cone_bipolar_trajectory.C20_6_11_10.PT_VSX2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##looking at the rod_bipolars during development
rod_bipolar_trajectory <- c("C20", "C6", "C11", "C9")
rod_bipolar_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "rod_bipolarU", 
                       groupBy = "Clusters",trajectory = rod_bipolar_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "OTX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "OTX2", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_OTX2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "CRX", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "CRX", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_CRX.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "PRDM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "PRDM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_PRDM1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NEUROD4", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NEUROD4", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_NEUROD4.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "GRM6", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "GRM6", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_GRM6.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "TRPM1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "TRPM1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_TRPM1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIA", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIA", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_NFIA.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIB", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIB", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_NFIB.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "NFIX", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "NFIX", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_NFIX.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneScoreMatrix", name = "VSX2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "rod_bipolarU", colorBy = "GeneIntegrationMatrix", name = "VSX2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "rod_bipolar_trajectory.C20_6_11_9.PT_VSX2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##looking at the amacrines during development
amacrine_trajectory <- c("C20", "C6", "C11", "C9")
amacrine_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "amacrineU", 
                       groupBy = "Clusters",trajectory = amacrine_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "ONECUT1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "ONECUT1", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_ONECUT1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "ONECUT2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "ONECUT2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_ONECUT2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "ONECUT3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "ONECUT3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_ONECUT3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "PTF1A", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "PTF1A", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_PTF1A.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "FOXN4", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "FOXN4", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_FOXN4.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "PAX6", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "PAX6", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_PAX6.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "NEUROD1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "NEUROD1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_NEUROD1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MEIS1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MEIS1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_MEIS1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MEIS2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MEIS2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_MEIS2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MEIS3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MEIS3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_MEIS3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "GAD1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "GAD1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_GAD1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "GAD2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "GAD2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_GAD2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "SLC6A9", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "SLC6A9", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_SLC6A9.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "amacrineU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "amacrine_trajectory.C20_6_11_9.PT_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

ganglion_trajectory2 <- c("C19", "C4", "C3", "C17", "C18")
ganglion_trajectory2
proj1 <- addTrajectory(ArchRProj = proj1, name = "ganglion2U", 
                       groupBy = "Clusters",trajectory = ganglion_trajectory2, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU4F1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU4F1", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_POU4F1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU4F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU4F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_POU4F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU4F3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU4F3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_POU4F3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "EBF1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "EBF1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_EBF1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "EBF2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "EBF2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_EBF2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "EBF3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "EBF3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_EBF3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneScoreMatrix", name = "RBPMS", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "ganglion2U", colorBy = "GeneIntegrationMatrix", name = "RBPMS", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "ganglion_trajectory.C19_4_3_17_18.PT_RBPMS.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##looking at the mullers during development
muller_trajectory <- c("C19", "C20", "C12")
muller_trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "mullerU", 
                       groupBy = "Clusters",trajectory = muller_trajectory, 
                       reducedDims = "Harmony", embedding = "UMAP", force = TRUE)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "LINC00461", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "LINC00461", continuousSet = "blueYellow",  log2norm = TRUE)
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_LINC00461.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MIR9-2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MIR9-2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_MIR9-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "POU6F2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "POU6F2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_POU6F2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MEIS1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MEIS1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_MEIS1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MEIS2", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MEIS2", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_MEIS2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "MEIS3", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "MEIS3", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_MEIS3.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "HES1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "HES1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_HES1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "NOTCH1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "NOTCH1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_NOTCH1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "SOX9", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "SOX9", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_SOX9.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneScoreMatrix", name = "RLBP1", continuousSet = "horizonExtra",  log2norm = TRUE)
p2 <- plotTrajectory(proj1, trajectory = "mullerU", colorBy = "GeneIntegrationMatrix", name = "RLBP1", continuousSet = "blueYellow",  log2norm = TRUE)
plotPDF(p1,p2, name = "muller_trajectory.C19_20_12_RLBP1.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)






