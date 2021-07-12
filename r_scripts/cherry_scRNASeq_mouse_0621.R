# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/seurat_analysis_0621')

# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

##load data
load('/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/data/scRNAseq_matrix/Mat_RNA_list_seurat') 
### see the time point: ####
names(Mat_RNA_list_seurat)
#### and get the individual seurat objects and sample into metadata
#head(x = P2[[]])
E11 = Mat_RNA_list_seurat$E11
E11$sample <- plyr::mapvalues(x = E11$orig.ident, from = c('RNA'), to = c('E11'))
E12 = Mat_RNA_list_seurat$E12
E12$sample = E12$orig.ident
E14 = Mat_RNA_list_seurat$E14
E14$sample = E14$orig.ident
E16 = Mat_RNA_list_seurat$E16
E16$sample <- plyr::mapvalues(x = E16$orig.ident, from = c('RNA'), to = c('E16'))
E18 = Mat_RNA_list_seurat$E18
E18$sample = E18$orig.ident
P0 = Mat_RNA_list_seurat$P0
P0$sample <- plyr::mapvalues(x = P0$orig.ident, from = c('RNA'), to = c('P0'))
P2 = Mat_RNA_list_seurat$P2
P2$sample = P2$orig.ident
P5 = Mat_RNA_list_seurat$P5
head(x = P5[[]])
P5$sample <- plyr::mapvalues(x = P5$orig.ident, from = c('RNA'), to = c('P5'))
P8 = Mat_RNA_list_seurat$P8
P8$sample = P8$orig.ident

# Merge into one single Seurat object
mouse=merge(E11, y=c(E12,E14,E16,E18,P0,P2,P5,P8))

##rm individual data
rm(Mat_RNA_list_seurat,E11,E12,E14,E16,E18,P0,P2,P5,P8)

##check labels
table(mouse$sample)

#Store mitochondrial percentage in the Seurat object metadata
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^mt-")

#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1, group.by = 'sample')
dev.copy2pdf(file="mouse_scrnaseq_all_0621.qc.pdf", width=20)   
   
#Filter the data
mouse <- subset(mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
mouse=NormalizeData(mouse)
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
mouse=ScaleData(mouse)
mouse = RunPCA(mouse)
DimPlot(mouse, reduction='pca', group.by='sample')
dev.copy2pdf(file="mouse_scrnaseq_all_0620.pca.pdf", width=20)
ElbowPlot(mouse, ndims=30)
dev.copy2pdf(file="mouse_scrnaseq_all_0620.elbow_plot.pdf", width=20)

#Save merged and filtered dataset as an R object to be able to re-load it for various future analyses without needing to perform the previous computations
saveRDS(mouse, file = "mouse.rds")
##and load if needed
#mouse <- readRDS(file = "mouse.rds")

#Run SCTransform, set var.to.regress to percent.mt --- ended up using harmony
mouse_merged <- SCTransform(mouse, vars.to.regress = "percent.mt", verbose = FALSE)
##resolution 0.3 and 30 dims
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
mouse_merged <- RunPCA(mouse_merged, verbose = FALSE)
mouse_merged <- RunUMAP(mouse_merged, dims = 1:30, verbose = FALSE)
mouse_merged <- FindNeighbors(mouse_merged, dims = 1:30, verbose=FALSE)
mouse_merged <- FindClusters(mouse_merged, resolution = 0.3)
#Plot UMAP
DimPlot(mouse_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.merged.UMAP.res0.3.sample_split.pdf", width = 20)
DimPlot(mouse_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.merged.UMAP.res0.3.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(mouse_merged$seurat_clusters, mouse_merged$sample)
write.csv(counts_cluster_sample, file='mouse_scrnaseq.merged.counts_cluster_sample.res0.3.csv')
##save
saveRDS(mouse_merged, file = "mouse_merged.rds")

#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
mouse_harmony <- RunHarmony(object = mouse, group.by.vars = 'sample', kmeans_init_nstart=20, kmeans_init_iter_max=100)
DimPlot(mouse_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="mouse_scrnaseq.harmony.harmony_plot.pdf", width = 20)

# Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
mouse_harmony <- RunUMAP(mouse_harmony, dims = 1:30, reduction = 'harmony')
mouse_harmony <- FindNeighbors(mouse_harmony, reduction = 'harmony', dims = 1:30)

##change resolutions - 0.3
mouse_harmony <- FindClusters(mouse_harmony, resolution = 0.3)
#Plot UMAP
DimPlot(mouse_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.harmony.UMAP.res0.3.sample_split.pdf", width = 20)
DimPlot(mouse_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="mouse_scrnaseq.harmony.UMAP.res0.3.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(mouse_harmony$seurat_clusters, mouse_harmony$sample)
write.csv(counts_cluster_sample, file='mouse_scrnaseq.harmony.counts_cluster_sample.res0.3.csv')


#Find all markers that define each cluster
mouse_harmony.markers <- FindAllMarkers(mouse_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mouse_harmony.markers, file='mouse_scrnaseq_all_0620.harmony.markers.res0.3.csv')






##modify this... not used in the end
##load object
human_harmony <- readRDS(file = "./seurat_analysis/all/human_harmony.rds")
##remove cluster 5, and check umap
human_harmony_subset <- subset(human_harmony, subset = seurat_clusters == 5, invert = TRUE)
DimPlot(human_harmony, reduction = "umap", pt.size = 0.1, label = TRUE)
DimPlot(human_harmony_subset, reduction = "umap", pt.size = 0.1, label = TRUE)

##make marked umap
new.cluster.ids <- c('Rods', 'Progenitors', 'Ganglions', 'Amacrines', 'Ganglions', 'Progenitors', 'Horizontals', 'Cones', 'Mullers', 'Progenitors', 'Ganglions', 'Ganglion precursors', 'BC/Photoreceptor precursors', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Progenitors', 'Progenitors', 'Ganglions', 'Amacrines', 'Astrocytes', 'Microglia', 'Progenitors')
names(new.cluster.ids) <- levels(human_harmony_subset)
human_harmony_subset <- RenameIdents(human_harmony_subset, new.cluster.ids)
DimPlot(human_harmony_subset, reduction = "umap", pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.harmony.UMAP_celltypes.res0.4.pdf", width=20, height=20)

##get barcode info for each cluster
WhichCells(object = human_harmony_subset, idents = c('Rods', 'Progenitors', 'Ganglions', 'Amacrines', 'Horizontals', 'Cones', 'Mullers', 'Ganglion precursors', 'BC/Photoreceptor precursors', 'Bipolars', 'Astrocytes', 'Microglia'))
human_harmony_subset$celltype <- plyr::mapvalues(
  x = human_harmony_subset$seurat_clusters, 
  from = c('0', '1', '2', '3', '4', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25'), 
  to = c('Rods', 'Progenitors', 'Ganglions', 'Amacrines', 'Ganglions', 'Progenitors', 'Horizontals', 'Cones', 'Mullers', 'Progenitors', 'Ganglions', 'Ganglion precursors', 'BC/Photoreceptor precursors', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Progenitors', 'Progenitors', 'Ganglions', 'Amacrines', 'Astrocytes', 'Microglia', 'Progenitors')
)
head(human_harmony_subset$celltype)
colnames(human_harmony_subset)
cell.info = human_harmony_subset$celltype
write.csv(cell.info, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.res0.4.cellinfo.csv')

##get counts per sample/cell class
sample.celltype.info = table(human_harmony_subset$celltype, human_harmony_subset$orig.ident)
write.csv(sample.celltype.info, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.res0.4.sample_cell_info.csv')







##save
saveRDS(mouse_harmony, file = "mouse_harmony.rds")
##and load if needed
#mouse_harmony <- readRDS(file = "mouse_harmony.rds")


