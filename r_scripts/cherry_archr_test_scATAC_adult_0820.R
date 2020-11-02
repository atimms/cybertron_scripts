##install archr and related packages
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
#library(ArchR)
#ArchR::installExtraPackages()
#devtools::install_github("GreenleafLab/chromVARmotifs")

##get archr and set up env
#library(ggplot2)
library(pheatmap)
library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")
##Set working directory and specific dir in active RSS where the data is sitting
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_test_0820')

#inputFiles <- getTutorialData("Hematopoiesis")
#inputFiles

##input frag files
inputFiles = c('Hu5' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu5/outs/fragments.tsv.gz',
               'Hu7' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu7/outs/fragments.tsv.gz',
               'Hu8' = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu8/outs/fragments.tsv.gz')
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
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "human_adult",
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
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##save archr project
#saveArchRProject(ArchRProj = proj1, outputDirectory = "save_project_v1", load = FALSE)
##filter doublets, can adjust to filter more is needed
proj1 <- filterDoublets(ArchRProj = proj1)
##Dimensionality Reduction and Clustering
##lsi if samples are pretty similiar
##The most common parameters to tweak are iterations, varFeatures, and resolution
proj1 <- addIterativeLSI(ArchRProj = proj1, useMatrix = "TileMatrix", name = "IterativeLSI", 
  iterations = 4, clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1,0.2,0.4), sampleCells = 10000, n.start = 10), 
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
plotPDF(plotList = p, name = "plot_marker_genes_with_imputation.harmony.res0.2.csv.pdf", 
        ArchRProj = proj1,  addDOC = FALSE, width = 5, height = 5)
##get tracks for marker genes and then print
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters", 
  geneSymbol = markerGenes, upstream = 50000, downstream = 50000)
plotPDF(plotList = p, name = "plot_tracks_for_marker_genes.harmony.res0.2.csv.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


##add cluster2 to names the cluster
cM <- confusionMatrix(proj1$Clusters, proj1$cellNames)
labelOld <- rownames(cM)
labelNew <- c("horizontals", "rods", "amacrines", "bipolars", "bipolars", "bipolars", "amacrines", "bipolars", "mullers", "bipolars", "bipolars", "ganglions", "bipolars", "cones")
proj1$Clusters2 <- mapLabels(proj1$Clusters, newLabels = labelNew, oldLabels = labelOld)
p1 <- plotEmbedding(proj1, colorBy = "cellColData", name = "Clusters2")
plotPDF(p1, name = "UMAP_cluster_names.harmony.res0.2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##make pseudo-bulk replicates - need to change defaults?
proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "Clusters2")

##call peaks with macs2
pathToMacs2 <- findMacs2()
proj1 <- addReproduciblePeakSet(ArchRProj = proj1, 
  groupBy = "Clusters2",pathToMacs2 = pathToMacs2)
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
write.csv(markerList, file='marker_list_cell_class.harmony.res0.2.csv')
##graph marker peaks
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmap_cell_class_marker_peaks.harmony.res0.2.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)
##get tracks for marker genes and then print
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "Clusters2", 
                      geneSymbol = markerGenes, upstream = 50000, downstream = 50000)
plotPDF(plotList = p, name = "plot_tracks_for_marker_genes_cell_class.harmony.res0.2.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##pairwise test
markerTest <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "PeakMatrix",
  groupBy = "Clusters2",testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "rods", bgdGroups = "cones")
pv <- markerPlot(seMarker = markerTest, name = "rods", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv
plotPDF(pv, name = "rods_vs_cones_volcano.harmony.res0.2.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)

##Motif and Feature Enrichment
##add motif
proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")
##compare previous pairwise test (rods vs cones)
motifsUp <- peakAnnoEnrichment(seMarker = markerTest,
  ArchRProj = proj1, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
motifsUp
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
motifsDo <- peakAnnoEnrichment(seMarker = markerTest, ArchRProj = proj1,
  peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
plotPDF(ggUp, ggDo, name = "rods_vs_cones_enriched_motifs.harmony.res0.2.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
##look for motif enrichment in all cell classes marker peaks
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = proj1,
  peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "heatmap_enriched_motifs_markers.harmony.res0.2.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)

##alternative archr enrichemnst
##Encode TF Binding Sites
proj1 <- addArchRAnnotations(ArchRProj = proj1, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj1,
  peakAnnotation = "EncodeTFBS",cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEncode, name = "heatmap_enriched_EncodeTFBS.harmony.res0.2.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)
##codex TF Binding Sites
proj1 <- addArchRAnnotations(ArchRProj = proj1, collection = "Codex")
enrichCodex <- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = proj1,
  peakAnnotation = "Codex", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapCodex, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapCodex, name = "heatmap_enriched_CodexTFBS.harmony.res0.2.pdf", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)

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
plotPDF(plotVarDev, name = "variable_motif_deviation_scores.harmony.res0.2.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
## extract a subset of motifs for downstream analysis (most variable)
motifs <- c("JUN", "FOS", "BACH", "TAL", "SMARCC")
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
plotPDF(p, name = "groups_deviations_imputation.harmony.res0.2.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
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
plotPDF(p, name = "motif_umaps_zscores.harmony.res0.2.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)


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
plotPDF(p, name = "motif_umaps_inferred_gene_expression.harmony.res0.2.pdf", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)

##motif footprinting
##get motif postion
motifPositions <- getPositions(proj1)
motifPositions
##get motifs we're interested in
motifs <- c("JUN", "FOS", "BACH", "TAL", "SMARCC")
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
proj1 <- addCoAccessibility(ArchRProj = proj1, reducedDims = "IterativeLSI")
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
  downstream = 50000,loops = getCoAccessibility(proj1))
plotPDF(plotList = p, name = "plot_tracks_marker_genes_with_CoAccessibility.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##trajectory analysis
##looking at the bipolars
trajectory <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7")
trajectory
proj1 <- addTrajectory(ArchRProj = proj1, name = "Bipolar", 
  groupBy = "Clusters",trajectory = trajectory, 
  embedding = "UMAP", force = TRUE)
##each cell has a pseudotime value between 0 and 100
head(proj1$Bipolar[!is.na(proj1$Bipolar)])
p <- plotTrajectory(proj1, trajectory = "Bipolar", colorBy = "cellColData", name = "Bipolar")
plotPDF(p, name = "plot_bipolar-traj_UMAP.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##pseudotime compared with LHX4 inferred expression
p1 <- plotTrajectory(proj1, trajectory = "Bipolar", colorBy = "GeneScoreMatrix", name = "LHX4", continuousSet = "horizonExtra")
ggAlignPlots(p1[[1]], p1[[2]], type = "h")
plotPDF(p1, name = "bipolar-traj_LHX4_inferred_expression.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##bipolar pseudo-time heatmaps, comparing to motifs, gene score and peak accessibility
trajMM  <- getTrajectory(ArchRProj = proj1, name = "Bipolar", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj1, name = "Bipolar", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajPM  <- getTrajectory(ArchRProj = proj1, name = "Bipolar", useMatrix = "PeakMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p3, name = "bipolar_traj_heatmaps.pdf", ArchRProj = proj1, addDOC = FALSE, width = 6, height = 8)
##Integrative pseudo-time analyses
##compare genescore and motifs trajectories 
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]
##get significant objects, then foramt
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1, ht2, name = "bipolar_traj_sig_genescore_motif.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


##save archr project in dir first specified
saveArchRProject(ArchRProj = proj1)
proj1 = loadArchRProject(path = 'human_adult')



































###other versions/random stuff

##resoltion 0.4 version
proj1 <- addClusters(input = proj1, reducedDims = "Harmony",
                     method = "Seurat", name = "Clusters", resolution = 0.4, force = TRUE)
counts_cluster_sample = table(proj1$Clusters, proj1$Sample)
write.csv(counts_cluster_sample, file='counts_cluster_sample.harmony.res0.4.csv')
##umap using harmony
proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "Harmony", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "UMAP_sample_clusters.harmony.res0.4.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
##finding marker genes
markersGS <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "GeneScoreMatrix", 
                               groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6
write.csv(markerList, file='marker_list.harmony.res0.4.csv')



##investigate cluster4.... not needed again
proj1
proj1$Clusters
unique(proj1$Clusters)
idxSample <- BiocGenerics::which(proj1$Clusters %in% "C4")
cellsSample <- proj1$cellNames[idxSample]
proj_c4 = proj1[cellsSample, ]
idxSample <- BiocGenerics::which(proj1$Clusters %in% "C5")
cellsSample <- proj1$cellNames[idxSample]
proj_c5 = proj1[cellsSample, ]
idxSample <- BiocGenerics::which(proj1$Clusters %in% "C10")
cellsSample <- proj1$cellNames[idxSample]
proj_c10 = proj1[cellsSample, ]
mean(proj1$NucleosomeRatio)
mean(proj_c4$NucleosomeRatio)
mean(proj1$ReadsInBlacklist)
mean(proj_c4$ReadsInBlacklist)
df <- getCellColData(proj1, select = c("log10(nFrags)", "TSSEnrichment"))
p1 <- ggPoint(x = df[,1], y = df[,2], colorDensity = TRUE,
  continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment", xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))) + geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")
p1
df2 <- getCellColData(proj_c4, select = c("log10(nFrags)", "TSSEnrichment"))
p2 <- ggPoint(x = df2[,1], y = df2[,2], colorDensity = TRUE,
              continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments",
              ylabel = "TSS Enrichment", xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
              ylim = c(0, quantile(df[,2], probs = 0.99))) + geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")
p2
df3 <- getCellColData(proj_c5, select = c("log10(nFrags)", "TSSEnrichment"))
p3 <- ggPoint(x = df3[,1], y = df3[,2], colorDensity = TRUE,
              continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments",
              ylabel = "TSS Enrichment", xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
              ylim = c(0, quantile(df[,2], probs = 0.99))) + geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")
p3
df4 <- getCellColData(proj_c10, select = c("log10(nFrags)", "TSSEnrichment"))
p4 <- ggPoint(x = df4[,1], y = df4[,2], colorDensity = TRUE,
              continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments",
              ylabel = "TSS Enrichment", xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
              ylim = c(0, quantile(df[,2], probs = 0.99))) + geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")
p4
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p1, p3, p4, p2, cols=4)
dev.copy2pdf(file="cluster4_comparison.harmony.res0.4.pdf", width=20)
