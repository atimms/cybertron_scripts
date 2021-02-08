library(dplyr)
library(Seurat)
library(patchwork)

setwd('/home/atimms/ngs_data/misc/cherry_mouse_scrnaseq_0121')

##need this to use umap
#reticulate::py_install(packages ='umap-learn')

##load data
mouse_data <- read.table("GSE63472_P14Retina_merged_digital_expression.txt",header=T,row.names = 1)

# Initialize the Seurat object with the raw (non-normalized data).
seurat_mat <- CreateSeuratObject(counts = mouse_data, project = "mouse_retina_GSE63472", min.cells = 3, min.features = 200)
seurat_mat


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter the data --- ? correct
seurat_mat <- subset(seurat_mat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

##normalize data
seurat_mat <- NormalizeData(seurat_mat)
##find variable genes, default is top 2k
seurat_mat <- FindVariableFeatures(seurat_mat, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_mat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_mat)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))
# Scaling the data
all.genes <- rownames(seurat_mat)
seurat_mat <- ScaleData(seurat_mat, features = all.genes)
#Perform linear dimensional reduction
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(object = seurat_mat))
# Examine and visualize PCA results a few different ways
print(seurat_mat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat_mat, dims = 1:2, reduction = "pca")
DimPlot(seurat_mat, reduction = "pca")

#ranking of principle components based on the percentage of variance explained by each one
ElbowPlot(seurat_mat)

#Cluster the cells, using 10 pcs in this case
seurat_mat <- FindNeighbors(seurat_mat, dims = 1:30)
seurat_mat <- FindClusters(seurat_mat, resolution = 0.4)
# Look at cluster IDs of the first 5 cells
head(Idents(seurat_mat), 5)


#Run non-linear dimensional reduction (UMAP/tSNE)
seurat_mat <- RunUMAP(seurat_mat, dims = 1:20)
DimPlot(seurat_mat, reduction = "umap", label = TRUE)
dev.copy2pdf(file="GSE63472_P14Retina.umap.pdf", width = 10, height = 10)
DimPlot(seurat_mat, reduction = "umap", group.by = "orig.ident", label = TRUE)
dev.copy2pdf(file="GSE63472_P14Retina.umap.sample_id.pdf", width = 10, height = 10)


# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_mat.markers <- FindAllMarkers(seurat_mat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(seurat_mat.markers, file="GSE63472_P14Retina.cluster_markers.csv")

#Dot plot - the size of the dot = % of cells and color represents the average expression
markers = c('LHX1', 'SLC17A6', 'PAX6', 'GAD1', 'SLC6A9', 'OPN1MW', 'VSX2', 'RLBP1', 'GFAP', 'PECAM1', 'KCNJ8', 'CX3CR1', 'HTRA1')
DotPlot(seurat_mat, features = markers) + RotatedAxis()
dev.copy2pdf(file="GSE63472_P14Retina.dotplot.pdf", width = 12, height = 6)

##name the clusters
new.cluster.ids <- c('Rods', 'Rods', 'Bipolars', 'Cones', 'Mullers', 'Amacrines', 'Bipolars', 'Bipolars', 'Amacrines', 'Bipolars', 'Rods', 'Rods', 'Bipolars', 'Rods', 'Amacrines', 'Amacrines', 'Ganglions', 'Rods', 'Rods', 'Rods', 'Vascular', 'Amacrines', 'Horizontals', 'Amacrines', 'Amacrines', 'Rods', 'Mullers', 'Amacrines', 'Microglia', 'Astrocytes')
names(new.cluster.ids) <- levels(seurat_mat)
seurat_mat <- RenameIdents(seurat_mat, new.cluster.ids)
DimPlot(seurat_mat, reduction = "umap", pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="GSE63472_P14Retina.umap.cell_class_names.pdf", width=10, height=10)
seurat_mat$celltype <- plyr::mapvalues(
  x = seurat_mat$seurat_clusters, 
  from = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29'), 
  to = c('Rods', 'Rods', 'Bipolars', 'Cones', 'Mullers', 'Amacrines', 'Bipolars', 'Bipolars', 'Amacrines', 'Bipolars', 'Rods', 'Rods', 'Bipolars', 'Rods', 'Amacrines', 'Amacrines', 'Ganglions', 'Rods', 'Rods', 'Rods', 'Vascular', 'Amacrines', 'Horizontals', 'Amacrines', 'Amacrines', 'Rods', 'Mullers', 'Amacrines', 'Microglia', 'Astrocytes')
)


##feature plot and violin for HTRA1
FeaturePlot(seurat_mat, features = c('HTRA1'))
dev.copy2pdf(file="GSE63472_P14Retina.HTRA1.FeaturePlot.pdf", width = 7, height =7)
VlnPlot(seurat_mat,c("HTRA1"), pt.size = 0)
dev.copy2pdf(file="GSE63472_P14Retina.HTRA1.ViolinPlot.pdf", width = 10, height =7)


##want to make the violin plot comparable to the humans data, so don't use all cell types and order
# Define an order of cluster identities and then factor the metadata column
my_levels <- c('Amacrines', 'Astrocytes', 'Bipolars', 'Cones', 'Ganglions', 'Horizontals', 'Microglia', 'Mullers', 'Rods')
seurat_mat$celltype <- factor(x = seurat_mat$celltype, levels = my_levels)
VlnPlot(seurat_mat,c("HTRA1"), pt.size = 0, group.by = "celltype", idents =c('Amacrines', 'Astrocytes', 'Bipolars', 'Cones', 'Ganglions', 'Horizontals', 'Microglia', 'Mullers', 'Rods'))
dev.copy2pdf(file="GSE63472_P14Retina.HTRA1.ViolinPlot_formatted.pdf", width = 10, height =7)
VlnPlot(seurat_mat,c("PLEKHA1"), pt.size = 0, group.by = "celltype", idents =c('Amacrines', 'Astrocytes', 'Bipolars', 'Cones', 'Ganglions', 'Horizontals', 'Microglia', 'Mullers', 'Rods'))
dev.copy2pdf(file="GSE63472_P14Retina.PLEKHA1.ViolinPlot_formatted.pdf", width = 10, height =7)


##save file
saveRDS(seurat_mat, file = "mouse_retina_GSE63472_0121.rds")

 ##and then load
seurat_mat <- readRDS(file ="mouse_retina_GSE63472_0121.rds")



