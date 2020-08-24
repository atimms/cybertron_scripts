library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

workingDir = "/home/atimms/ngs_data/cellranger/kim_visium_0820";
setwd(workingDir);

##for each sample load and normalize i.e sctransform
##JAS_A1
##load data
JAS_A1 = Load10X_Spatial('/home/atimms/ngs_data/cellranger/kim_visium_0820/Results/JAS-A1_Visium_results')
##add new column to metadata named sample (as orig.ident is generic)
new_names <- JAS_A1$orig.ident
new_names <- gsub("SeuratProject", "JAS_A1", new_names)
JAS_A1$sample <- new_names
head(JAS_A1[[]])
##molecular counts / spot
plot1 <- VlnPlot(JAS_A1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(JAS_A1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.copy2pdf(file='kim_visium_0820.JAS_A1.counts_spot.pdf', width = 12, height = 8)
##remove cells with less than 100 counts
JAS_A1 <- subset(JAS_A1, subset = nCount_Spatial > 100)
##run sctransform
JAS_A1 <- SCTransform(JAS_A1, assay = "Spatial", verbose = FALSE)
##PCA, UMAP and Spatial Plot
JAS_A1 <- RunPCA(JAS_A1, assay = "SCT", verbose = FALSE)
JAS_A1 <- FindNeighbors(JAS_A1, reduction = "pca", dims = 1:30)
JAS_A1 <- FindClusters(JAS_A1, verbose = FALSE)
JAS_A1 <- RunUMAP(JAS_A1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(JAS_A1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(JAS_A1, label = TRUE, label.size = 3)
p1 + p2
dev.copy2pdf(file='kim_visium_0820.JAS_A1.UMAP_plus.pdf', width = 12, height = 8)
##find spatial variable features
#JAS_A1 <- FindSpatiallyVariableFeatures(JAS_A1, assay = "SCT", features = VariableFeatures(JAS_A1)[1:1000], 
#                                       selection.method = "markvariogram")
#top.features <- head(SpatiallyVariableFeatures(JAS_A1, selection.method = "markvariogram"), 6)
#SpatialFeaturePlot(JAS_A1, features = top.features, ncol = 3, alpha = c(0.1, 1))
#dev.copy2pdf(file='kim_visium_0820.JAS_A1.top6_variable_genes_plot.pdf', width = 12, height = 8)

##JAS_B1
##load data
JAS_B1 = Load10X_Spatial('/home/atimms/ngs_data/cellranger/kim_visium_0820/Results/JAS-B1_Visium_results')
##add new column to metadata named sample (as orig.ident is generic)
new_names <- JAS_B1$orig.ident
new_names <- gsub("SeuratProject", "JAS_B1", new_names)
JAS_B1$sample <- new_names
head(JAS_B1[[]])
##molecular counts / spot
plot1 <- VlnPlot(JAS_B1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(JAS_B1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.copy2pdf(file='kim_visium_0820.JAS_B1.counts_spot.pdf', width = 12, height = 8)
##remove cells with less than 100 counts
JAS_B1 <- subset(JAS_B1, subset = nCount_Spatial > 100)
##run sctransform
JAS_B1 <- SCTransform(JAS_B1, assay = "Spatial", verbose = FALSE)
##PCA, UMAP and Spatial Plot
JAS_B1 <- RunPCA(JAS_B1, assay = "SCT", verbose = FALSE)
JAS_B1 <- FindNeighbors(JAS_B1, reduction = "pca", dims = 1:30)
JAS_B1 <- FindClusters(JAS_B1, verbose = FALSE)
JAS_B1 <- RunUMAP(JAS_B1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(JAS_B1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(JAS_B1, label = TRUE, label.size = 3)
p1 + p2
dev.copy2pdf(file='kim_visium_0820.JAS_B1.UMAP_plus.pdf', width = 12, height = 8)
##find spatial variable features
#JAS_B1 <- FindSpatiallyVariableFeatures(JAS_B1, assay = "SCT", features = VariableFeatures(JAS_B1)[1:1000], 
#                                        selection.method = "markvariogram")
#top.features <- head(SpatiallyVariableFeatures(JAS_B1, selection.method = "markvariogram"), 6)
#SpatialFeaturePlot(JAS_B1, features = top.features, ncol = 3, alpha = c(0.1, 1))
#dev.copy2pdf(file='kim_visium_0820.JAS_B1.top6_variable_genes_plot.pdf', width = 12, height = 8)

##JAS_C1
##load data
JAS_C1 = Load10X_Spatial('/home/atimms/ngs_data/cellranger/kim_visium_0820/Results/JAS-C1_Visium_results')
##add new column to metadata named sample (as orig.ident is generic)
new_names <- JAS_C1$orig.ident
new_names <- gsub("SeuratProject", "JAS_C1", new_names)
JAS_C1$sample <- new_names
head(JAS_C1[[]])
##molecular counts / spot
plot1 <- VlnPlot(JAS_C1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(JAS_C1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.copy2pdf(file='kim_visium_0820.JAS_C1.counts_spot.pdf', width = 12, height = 8)
##remove cells with less than 100 counts
JAS_C1 <- subset(JAS_C1, subset = nCount_Spatial > 100)
##run sctransform
JAS_C1 <- SCTransform(JAS_C1, assay = "Spatial", verbose = FALSE)
##PCA, UMAP and Spatial Plot
JAS_C1 <- RunPCA(JAS_C1, assay = "SCT", verbose = FALSE)
JAS_C1 <- FindNeighbors(JAS_C1, reduction = "pca", dims = 1:30)
JAS_C1 <- FindClusters(JAS_C1, verbose = FALSE)
JAS_C1 <- RunUMAP(JAS_C1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(JAS_C1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(JAS_C1, label = TRUE, label.size = 3)
p1 + p2
dev.copy2pdf(file='kim_visium_0820.JAS_C1.UMAP_plus.pdf', width = 12, height = 8)
##find spatial variable features
#JAS_C1 <- FindSpatiallyVariableFeatures(JAS_C1, assay = "SCT", features = VariableFeatures(JAS_C1)[1:1000], 
#                                        selection.method = "markvariogram")
#top.features <- head(SpatiallyVariableFeatures(JAS_C1, selection.method = "markvariogram"), 6)
#SpatialFeaturePlot(JAS_C1, features = top.features, ncol = 3, alpha = c(0.1, 1))
#dev.copy2pdf(file='kim_visium_0820.JAS_C1.top6_variable_genes_plot.pdf', width = 12, height = 8)

##JAS_D1
##load data
JAS_D1 = Load10X_Spatial('/home/atimms/ngs_data/cellranger/kim_visium_0820/Results/JAS-D1_Visium_results')
##add new column to metadata named sample (as orig.ident is generic)
new_names <- JAS_D1$orig.ident
new_names <- gsub("SeuratProject", "JAS_D1", new_names)
JAS_D1$sample <- new_names
head(JAS_D1[[]])
##molecular counts / spot
plot1 <- VlnPlot(JAS_D1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(JAS_D1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.copy2pdf(file='kim_visium_0820.JAS_D1.counts_spot.pdf', width = 12, height = 8)
##remove cells with less than 100 counts
JAS_D1 <- subset(JAS_D1, subset = nCount_Spatial > 100)
##run sctransform
JAS_D1 <- SCTransform(JAS_D1, assay = "Spatial", verbose = FALSE)
##PCA, UMAP and Spatial Plot
JAS_D1 <- RunPCA(JAS_D1, assay = "SCT", verbose = FALSE)
JAS_D1 <- FindNeighbors(JAS_D1, reduction = "pca", dims = 1:30)
JAS_D1 <- FindClusters(JAS_D1, verbose = FALSE)
JAS_D1 <- RunUMAP(JAS_D1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(JAS_D1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(JAS_D1, label = TRUE, label.size = 3)
p1 + p2
dev.copy2pdf(file='kim_visium_0820.JAS_D1.UMAP_plus.pdf', width = 12, height = 8)
##find spatial variable features
#JAS_D1 <- FindSpatiallyVariableFeatures(JAS_D1, assay = "SCT", features = VariableFeatures(JAS_D1)[1:1000], 
#                                        selection.method = "markvariogram")
#top.features <- head(SpatiallyVariableFeatures(JAS_D1, selection.method = "markvariogram"), 6)
#SpatialFeaturePlot(JAS_D1, features = top.features, ncol = 3, alpha = c(0.1, 1))
#dev.copy2pdf(file='kim_visium_0820.JAS_D1.top6_variable_genes_plot.pdf', width = 12, height = 8)


##merge data
brain.merge <- merge(JAS_A1, y=c(JAS_B1, JAS_C1, JAS_D1), add.cell.ids = c("JAS_A1", "JAS_B1", "JAS_C1", "JAS_D1"))
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(JAS_A1), VariableFeatures(JAS_B1),VariableFeatures(JAS_C1), VariableFeatures(JAS_D1))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)
##plot data
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "sample"))
dev.copy2pdf(file='kim_visium_0820.merged.UMAP.pdf', width = 10, height = 6)
SpatialDimPlot(brain.merge)
dev.copy2pdf(file='kim_visium_0820.merged.spatial_plot.pdf', width = 10, height = 4)

#Find all markers that define each cluster, does this work
brain.merge.markers <- FindAllMarkers(brain.merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(brain.merge.markers, file='kim_visium_0820.merged.all_markers.csv')

#Save object
saveRDS(brain.merge, file = "brain_merged.rds")
