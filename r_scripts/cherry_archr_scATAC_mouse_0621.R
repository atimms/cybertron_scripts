#sessionInfo()
#.libPaths()

##error using lsi, need archr 1.02
#devtools::install_github("GreenleafLab/ArchR", 
#                         ref="release_1.0.2", repos = BiocManager::repositories())
##also need harmony
#devtools::install_github("immunogenomics/harmony")
##for presto to get markers
#ArchR::installExtraPackages()

##archr set up by Marc using 4.03
#library(ggplot2)
#library(GenomicRanges)
#library(chromVARmotifs) do I need this?
#library(pheatmap)
library(ArchR)
#library(ComplexHeatmap)
#library(circlize)
#set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("mm10")
##Set working directory and specific dir in active RSS where the data is sitting
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/archr_analysis_0521')

##input frag files
inputFiles = c('E11' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/E11_fragments_cl.tsv.gz',
               'E12' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/E12_fragments_cl.tsv.gz',
               'E14' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/E14_fragments_cl.tsv.gz',
               'E16' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/E16_fragments_cl.tsv.gz',
               'E18' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/E18_fragments_cl.tsv.gz',
               'P0' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/P0_fragments_cl.tsv.gz',
               'P11' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/P11_fragments_cl.tsv.gz',
               'P14' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/P14_fragments_cl.tsv.gz',
               'P2' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/P2_fragments_cl.tsv.gz',
               'P5' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/P5_fragments_cl.tsv.gz',
               'P8' = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_fragments/P8_fragments_cl.tsv.gz')
inputFiles


##make arrow file
ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = names(inputFiles),
                               minTSS = 4, #Dont set this too high because you can always increase later
                               minFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE)
ArrowFiles
##doublet inferrance
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                               knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                               LSIMethod = 1)
##Creating An ArchRProject
proj1 <- ArchRProject(ArrowFiles = ArrowFiles, 
                      outputDirectory = "mouse_combined",
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

##filter doublets, can adjust to filter more is needed
proj1 <- filterDoublets(ArchRProj = proj1)

##add cell id as metadata
Meta_data <- read.table('/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scATACseq_annotation/cell_info.txt', header=T, sep='\t', comment.char = "", check.names = FALSE)
m = match(proj1$cellNames,Meta_data$cell_id)
proj1$celltypes = as.character(Meta_data$celltypes[m])

##need to remove na's in cell types and subet
idxPass <- which(!is.na(proj1$celltypes))
cellsPass <- proj1$cellNames[idxPass]
proj2 <- subsetArchRProject(ArchRProj = proj1, cells = cellsPass,
                                    outputDirectory = "mouse_subset", dropCells = F, force = TRUE )


##change sample cell from 10k to 50k?
proj2 <- addIterativeLSI(ArchRProj = proj2, useMatrix = "TileMatrix", name = "IterativeLSI", 
                         iterations = 4, clusterParams = list( #See Seurat::FindClusters
                           resolution = c(0.1,0.2,0.4), sampleCells = 50000, n.start = 10), 
                         varFeatures = 25000, dimsToUse = 1:30)
##use harmony if sample need more adjustment
proj2 <- addHarmony(ArchRProj = proj2, reducedDims = "IterativeLSI",
                    name = "Harmony", groupBy = "Sample")
##find cluster using seurat (can use LSI or harmony) -- force is needed is repeating
proj2 <- addClusters(input = proj2, reducedDims = "Harmony",
                     method = "Seurat", name = "Clusters", resolution = 0.2, force = TRUE)
counts_cluster_sample = table(proj2$Clusters, proj2$Sample)
write.csv(counts_cluster_sample, file='counts_cluster_sample.harmony.res0.2.csv')
##graph
cM <- confusionMatrix(paste0(proj2$Clusters), paste0(proj2$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
                        border_color = "black")
p
plotPDF(plotList = p, name = "cluster_sample_heatmap.harmony.res0.2.pdf", 
        ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

##umap using harmony
proj2 <- addUMAP(ArchRProj = proj2, reducedDims = "Harmony", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
ggAlignPlots(p1, p2,p3, type = "h")
plotPDF(p1,p2,p3, name = "UMAP_sample_clusters.harmony.res0.2.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

##umap using lsi
proj2 <- addUMAP(ArchRProj = proj2, reducedDims = "IterativeLSI", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
ggAlignPlots(p1, p2,p3, type = "h")
plotPDF(p1,p2,p3, name = "UMAP_sample_clusters.lsi.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)







##call peaks using celltypes
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "celltypes", minCells = 20)
pathToMacs2 <- findMacs2()
proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "celltypes", 
  pathToMacs2 = pathToMacs2, cutOff = 0.000001)

##Add the peak matrix
proj2 <- addPeakMatrix(proj2)
getAvailableMatrices(proj2) #check for matrix

#get peaks specific to cell class
markersPeaks <- getMarkerFeatures(ArchRProj = proj2, 
                                  useMatrix = "PeakMatrix", groupBy = "celltypes",
                                  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
##get peaks per cell class
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
write.table(markerList, "marker_list.peaks_per_cell_class.filtered.0621.txt", append = TRUE, sep = "\t")
markerList <- getMarkers(markersPeaks)
markerList
write.table(markerList, "marker_list.peaks_per_cell_class.0621.txt", append = TRUE, sep = "\t")

##graph marker peaks
heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, 
                                  cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmap_cell_class_marker_peaks.macs2_q0.000001.pdf", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)
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


##with get cell/peak matrix (not interesting)
#test <- getMatrixFromProject(proj1, useMatrix='PeakMatrix')
#assays(test)$PeakMatrix

##save/load archr project
saveArchRProject(ArchRProj = proj2)
proj2 = loadArchRProject(path = 'mouse_subset')

##integrate scRNASeq
##get seurat object
seRNA =  readRDS('/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/seurat_analysis_0621/mouse_harmony.rds')
seRNA
##Unconstrained Integration
proj2 <- addGeneIntegrationMatrix(ArchRProj = proj2, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix",
                                  reducedDims = "Harmony", seRNA = seRNA, addToArrow = FALSE, groupRNA = "seurat_clusters",
                                  nameCell = "predictedCell", nameGroup = "predictedGroup", nameScore = "predictedScore")
##make color palatte and then plot
pal <- paletteDiscrete(values = seRNA$seurat_clusters)
p1 <- plotEmbedding(proj2, colorBy = "cellColData", name = "predictedGroup", pal = pal, rastr = FALSE)
plotPDF(p1, name = "UMAP_clusters_rnaseq.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
##looks good so repeat but change addToArrow = T and force = T
#addArchRThreads(threads = 1) ##only seems to work with single thread 
proj2 <- addGeneIntegrationMatrix(ArchRProj = proj2, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix",
                                  reducedDims = "Harmony", seRNA = seRNA, addToArrow = TRUE, force = TRUE, groupRNA = "seurat_clusters",
                                  nameCell = "predictedCell", nameGroup = "predictedGroup", nameScore = "predictedScore")
#addArchRThreads(threads = 10) ##switch back to 10 threads
getAvailableMatrices(proj2) ##check we have GeneIntegrationMatrix

##save/load archr project
saveArchRProject(ArchRProj = proj2)
proj2 = loadArchRProject(path = 'mouse_subset')

##get peak to gene data and atac/rna index info - repeat and save in current dir
proj2 <- addPeak2GeneLinks(ArchRProj = proj2,reducedDims = "Harmony", maxDist = 300000)
p2g <- getPeak2GeneLinks(ArchRProj = proj2, corCutOff = 0.4,
                         resolution = 10000, returnLoops = F)
write.table(p2g, 'mouse_all.min_cell20.p2g.m2_q0.000001.cor4.txt', row.names = F, sep="\t", quote = FALSE)
rnaidx.info = metadata(p2g)[[2]]
rnaidx.df = data.frame(rnaidx.info)
write.table(rnaidx.df, 'mouse_all.min_cell20.m2_q0.000001.rnaseq_info.txt', row.names = T, sep="\t", quote = FALSE)
peaks.gr <- getPeakSet(proj2)
df.peaks.gr = data.frame(peaks.gr)
write.table(df.peaks.gr, 'mouse_all.min_cell20.peak_info.m2_q0.000001.txt', row.names = T, sep="\t", quote = FALSE)


##get bigwigs for each cell class
proj2 = loadArchRProject(path = 'mouse_subset')
getGroupBW(ArchRProj = proj2, groupBy = "celltypes")







##looks at peak calls (not used right now)
#https://github.com/GreenleafLab/ArchR/issues/445
#cone_peak_file = '/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/archr_analysis_0521/mouse_combined/PeakCalls/Cone-reproduciblePeaks.gr.rds'
#cone_peaks = readRDS(cone_peak_file)
#write.table(cone_peaks, "cone_peaks.txt", append = TRUE, sep = "\t")

#incorporate scRNA and get peak2gene data





