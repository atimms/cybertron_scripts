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
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all_third')

##load the ArchRProject
proj1 = loadArchRProject(path = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/old/all_third_backup_092520/human_all')

##subset mature bipolars
idxPass <- which(proj1$Clusters2 == 'Mature Bipolars')
cellsPass <- proj1$cellNames[idxPass]
biopolar_proj <- subsetArchRProject(ArchRProj = proj1, cells = cellsPass,
  outputDirectory = "bipolar_subset")
##find cluster using seurat (can use LSI or harmony) -- force is needed is repeating
biopolar_proj <- addClusters(input = biopolar_proj, reducedDims = "Harmony",
                             method = "Seurat", name = "Clusters", resolution = 0.2, force = TRUE)
##graph
cM <- confusionMatrix(paste0(biopolar_proj$Clusters), paste0(biopolar_proj$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
                        border_color = "black")
p
plotPDF(plotList = p, name = "bipolars.cluster_sample_heatmap.harmony.res0.2.pdf", 
        ArchRProj = biopolar_proj, addDOC = FALSE, width = 5, height = 5)
counts_cluster_sample = table(biopolar_proj$Clusters, biopolar_proj$Sample)
write.csv(counts_cluster_sample, file='counts_cluster_sample.harmony.res0.2.csv')
##umap using harmony
biopolar_proj <- addUMAP(ArchRProj = biopolar_proj, reducedDims = "Harmony", 
                         name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
p1 <- plotEmbedding(ArchRProj = biopolar_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = biopolar_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "bipolar.UMAP_sample_clusters.harmony.res0.2.pdf", ArchRProj = biopolar_proj, addDOC = FALSE, width = 5, height = 5)
##finding marker genes
markersGS <- getMarkerFeatures(ArchRProj = biopolar_proj, useMatrix = "GeneScoreMatrix", 
                               groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.csv(markerList, file='biopolar.marker_list.harmony.res0.2.csv')

##harmony res 0.8
##find cluster using seurat (can use LSI or harmony) -- force is needed is repeating
biopolar_proj <- addClusters(input = biopolar_proj, reducedDims = "Harmony",
                             method = "Seurat", name = "Clusters", resolution = 0.8, force = TRUE)
##graph
cM <- confusionMatrix(paste0(biopolar_proj$Clusters), paste0(biopolar_proj$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
                        border_color = "black")
p
plotPDF(plotList = p, name = "bipolars.cluster_sample_heatmap.harmony.res0.8.pdf", 
        ArchRProj = biopolar_proj, addDOC = FALSE, width = 5, height = 5)
counts_cluster_sample = table(biopolar_proj$Clusters, biopolar_proj$Sample)
write.csv(counts_cluster_sample, file='counts_cluster_sample.harmony.res0.8.csv')
##umap using harmony
biopolar_proj <- addUMAP(ArchRProj = biopolar_proj, reducedDims = "Harmony", 
                         name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
p1 <- plotEmbedding(ArchRProj = biopolar_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = biopolar_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "bipolar.UMAP_sample_clusters.harmony.res0.8.pdf", ArchRProj = biopolar_proj, addDOC = FALSE, width = 5, height = 5)
##finding marker genes
markersGS <- getMarkerFeatures(ArchRProj = biopolar_proj, useMatrix = "GeneScoreMatrix", 
                               groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.csv(markerList, file='biopolar.marker_list.harmony.res0.8.csv')



##save/load project
saveArchRProject(ArchRProj = biopolar_proj)
biopolar_proj = loadArchRProject(path = 'bipolar_subset')
