# Set working directory 
setwd('/home/atimms/ngs_data/misc/cherry_scrublet_test_0820')

# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#Make sure Harmony library is installed
#library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)

#Import 10x data and convert each dataset to Seurat object (take from individual sample outputs, not cellranger's aggr output). Define sample with project= "samplename"
Hu37.data = Read10X(data.dir = './Hu37_scRNA')
Hu37 = CreateSeuratObject(counts = Hu37.data, project = "Hu37")
Hu37.doublets = read.table('./Hu37_scRNA/predicted_doublet_mask.txt')
Hu37$scrublet_doublet = unlist(Hu37.doublets)
head(Hu37[[]])
Hu37.noscrub = subset(Hu37, subset = scrublet_doublet == "False")
Hu7.data = Read10X(data.dir = './Hu7_scRNA')
Hu7 = CreateSeuratObject(counts = Hu7.data, project = "Hu7")
Hu7.doublets = read.table('./Hu7_scRNA/predicted_doublet_mask.txt')
Hu7$scrublet_doublet = unlist(Hu7.doublets)
head(Hu7[[]])
Hu7.noscrub = subset(Hu7, subset = scrublet_doublet == "False")
Hu5.data = Read10X(data.dir = './Hu5_scRNA')
Hu5 = CreateSeuratObject(counts = Hu5.data, project = "Hu5")
Hu5.doublets = read.table('./Hu5_scRNA/predicted_doublet_mask.txt')
Hu5$scrublet_doublet = unlist(Hu5.doublets)
head(Hu5[[]])
Hu5.noscrub = subset(Hu5, subset = scrublet_doublet == "False")

##merge data and check for all and no scrub
# Merge into one single Seurat object
human=merge(Hu37, y=c(Hu5,Hu7))
#Validate the merge by checking number of cells per group
table(human$orig.ident)
#Add sample and condition information explicitly into the metadata (as oppsoed to storing in 'orig.ident') for future downstream analysis
#Add sample info (simply copying 'orig.ident')
human$sample <- plyr::mapvalues(
  x = human$orig.ident, 
  from = c('Hu37', 'Hu5', 'Hu7'), 
  to = c('Hu37', 'Hu5', 'Hu7'))
#Validate new metadata columns by checking that number of cells per sample/phenotype adds up
table(human$sample)
#Store mitochondrial percentage in the Seurat object metadata
human[["percent.mt"]] <- PercentageFeatureSet(human, pattern = "^MT-")
#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
dev.copy2pdf(file="human.scrnaseq_adult_0820.qc.pdf", width=20)
#Filter the data
human <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)
##checking number of cells per group again
table(human$orig.ident)

##merge data and check for all and no scrub
# Merge into one single Seurat object
human.noscrub=merge(Hu37.noscrub, y=c(Hu5.noscrub,Hu7.noscrub))
#Validate the merge by checking number of cells per group
table(human.noscrub$orig.ident)
#Add sample and condition information explicitly into the metadata (as oppsoed to storing in 'orig.ident') for future downstream analysis
#Add sample info (simply copying 'orig.ident')
human.noscrub$sample <- plyr::mapvalues(
  x = human.noscrub$orig.ident, 
  from = c('Hu37', 'Hu5', 'Hu7'), 
  to = c('Hu37', 'Hu5', 'Hu7'))
#Validate new metadata columns by checking that number of cells per sample/phenotype adds up
table(human.noscrub$sample)
#Store mitochondrial percentage in the Seurat object metadata
human.noscrub[["percent.mt"]] <- PercentageFeatureSet(human.noscrub, pattern = "^MT-")
#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(human.noscrub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
dev.copy2pdf(file="human_scrubed.scrnaseq_adult_0820.qc.pdf", width=20)
#Filter the data
human.noscrub <- subset(human.noscrub, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)
##checking number of cells per group again
table(human.noscrub$orig.ident)

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human=NormalizeData(human)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)
human=ScaleData(human)
human = RunPCA(human)
DimPlot(human, reduction='pca', group.by='sample')
ElbowPlot(human, ndims=30)

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human.noscrub=NormalizeData(human.noscrub)
human.noscrub <- FindVariableFeatures(human.noscrub, selection.method = "vst", nfeatures = 2000)
human.noscrub=ScaleData(human.noscrub)
human.noscrub = RunPCA(human.noscrub)
DimPlot(human.noscrub, reduction='pca', group.by='sample')
ElbowPlot(human.noscrub, ndims=30)

#Save merged and filtered dataset as an R object to be able to re-load it for various future analyses without needing to perform the previous computations
#saveRDS(human, file = "human.rds")
#saveRDS(human.noscrub, file = "human.noscrub.rds")
##and load if needed
#human <- readRDS(file = "human.rds")
#human.noscrub <- readRDS(file = "human.noscrub.rds")

#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
human_harmony <- RunHarmony(object = human, group.by.vars = 'sample')
DimPlot(human_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="human.scrnaseq_adult_0820.harmony_plot.pdf", width = 20)
# Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
human_harmony <- RunUMAP(human_harmony, dims = 1:30, reduction = 'harmony')
human_harmony <- FindNeighbors(human_harmony, reduction = 'harmony', dims = 1:30)
##change resolutions - 0.4
human_harmony <- FindClusters(human_harmony, resolution = 0.4)
#Plot UMAP
DimPlot(human_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="human.scrnaseq_adult_0820.harmony.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="human.scrnaseq_adult_0820.harmony.UMAP.res0.4.pdf", width = 20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_harmony$seurat_clusters, human_harmony$sample)
write.csv(counts_cluster_sample, file='human.scrnaseq_adult_0820.harmony.counts_cluster_sample.res0.4.csv')
#Find all markers that define each cluster
human_harmony.markers <- FindAllMarkers(human_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_harmony.markers, file='human.scrnaseq_adult_0820.harmony.markers.res0.4.csv')

#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
human.noscrub_harmony <- RunHarmony(object = human.noscrub, group.by.vars = 'sample')
DimPlot(human.noscrub_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="human_scrubed.scrnaseq_adult_0820.harmony_plot.pdf", width = 20)
# Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
human.noscrub_harmony <- RunUMAP(human.noscrub_harmony, dims = 1:30, reduction = 'harmony')
human.noscrub_harmony <- FindNeighbors(human.noscrub_harmony, reduction = 'harmony', dims = 1:30)
##change resolutions - 0.4
human.noscrub_harmony <- FindClusters(human.noscrub_harmony, resolution = 0.4)
#Plot UMAP
DimPlot(human.noscrub_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="human_scrubed.scrnaseq_adult_0820.harmony.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human.noscrub_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="human_scrubed.scrnaseq_adult_0820.harmony.UMAP.res0.4.pdf", width = 20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human.noscrub_harmony$seurat_clusters, human.noscrub_harmony$sample)
write.csv(counts_cluster_sample, file='human_scrubed.scrnaseq_adult_0820.harmony.counts_cluster_sample.res0.4.csv')
#Find all markers that define each cluster
human.noscrub_harmony.markers <- FindAllMarkers(human.noscrub_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human.noscrub_harmony.markers, file='human_scrubed.scrnaseq_adult_0820.harmony.markers.res0.4.csv')

