
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
                      outputDirectory = "mouse_all",
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
proj1 = loadArchRProject(path = 'mouse_all')

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
                     method = "Seurat", name = "Clusters", resolution = 1.6, force = TRUE)
counts_cluster_sample = table(proj1$Clusters, proj1$Sample)
write.csv(counts_cluster_sample, file='counts_cluster_sample.harmony.res1.6.csv')
##graph
cM <- confusionMatrix(paste0(proj1$Clusters), paste0(proj1$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
                        border_color = "black")
p
plotPDF(plotList = p, name = "cluster_sample_heatmap.harmony.res1.6.pdf", 
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##umap using harmony
proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "Harmony", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "UMAP_sample_clusters.harmony.res1.6.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##get barcode ids per cluster to compare to the original annotation
barcodes_cluster <- getCellColData(proj1, select = "Clusters")
barcodes_cluster
write.csv(barcodes_cluster, file='barcodes_cluster.LSI.res1.6.csv')

##finding marker genes
markersGS <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "GeneScoreMatrix", 
                               groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6
write.csv(markerList, file='marker_list.harmony.res1.6.csv')

##oriniginal i.e. resolution 0.2
# find cluster using seurat (can use LSI or harmony) -- force is needed is repeating
# proj1 <- addClusters(input = proj1, reducedDims = "Harmony",
#                      method = "Seurat", name = "Clusters", resolution = 0.2, force = TRUE)
# counts_cluster_sample = table(proj1$Clusters, proj1$Sample)
# write.csv(counts_cluster_sample, file='counts_cluster_sample.harmony.res0.2.csv')
# graph
# cM <- confusionMatrix(paste0(proj1$Clusters), paste0(proj1$Sample))
# cM <- cM / Matrix::rowSums(cM)
# p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
#                         border_color = "black")
# p
# plotPDF(plotList = p, name = "cluster_sample_heatmap.harmony.res0.2.pdf", 
#         ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
# 
# umap using harmony
# proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "Harmony", 
#                  name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
# p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
# p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
# ggAlignPlots(p1, p2, type = "h")
# plotPDF(p1,p2, name = "UMAP_sample_clusters.harmony.res0.2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
# 
# get barcode ids per cluster to compare to the original annotation
# barcodes_cluster <- getCellColData(proj1, select = "Clusters")
# barcodes_cluster
# write.csv(barcodes_cluster, file='barcodes_cluster.LSI.res0.2.csv')
# 
# finding marker genes
# markersGS <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "GeneScoreMatrix", 
#                                groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
#                                testMethod = "wilcoxon")
# markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
# markerList$C6
# write.csv(markerList, file='marker_list.harmony.res0.2.csv')


##save/load archr project
saveArchRProject(ArchRProj = proj1)
proj1 = loadArchRProject(path = 'mouse_all')








##add cluster2 to names the cluster
cM <- confusionMatrix(proj1$Clusters, proj1$cellNames)
labelOld <- rownames(cM)
##check the weird order and replicate
labelOld
labelNew <- c('Early Progenitors', 'Ganglion Precursors', 'Developing Ganglions', 'AC/HC/GC Precursors', 'Late Progenitors', 'Mature Rods', 'Developing Rods', 'Photoreceptor/Bipolar Precursors', 'Developing Cones', 'Mature Horizontals', 'Mature Amacrines', 'Mature Bipolars', 'Mature Bipolars', 'Mature Mullers', 'Mature Cones', 'Mature Ganglions', 'Developing Amacrines', 'Developing Horizontals', 'Developing Bipolars', 'Developing Ganglions')
proj1$Clusters2 <- mapLabels(proj1$Clusters, newLabels = labelNew, oldLabels = labelOld)
p1 <- plotEmbedding(proj1, colorBy = "cellColData", name = "Clusters2", rastr = FALSE)
plotPDF(p1, name = "UMAP_cluster_names.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


