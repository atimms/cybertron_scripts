# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA')

# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Make sure Harmony library is installed
#library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)

##known markers
markers <- c('RCVRN', 'RHO', 'CRX', 'ARR3', 'GNAT2', 'VSX2', 'LHX4', 'TRPM1', 'GRM6', 'SLC1A3', 'RLBP1', 'PAX6', 'LHX1', 'ONECUT2', 'TFAP2B', 'GAD1', 'SLC6A9', 'RBPMS', 'NEFM', 'GFAP', 'CD74', 'P2RY12', 'BEST1', 'RPE65', 'SFRP2')
markers2 <- c('KCNQ5', 'PRDM1', 'DCT', 'SLC35F3', 'NOL4', 'MPPED2', 'WDR72', 'NLK', 'PDE1C', 'PEX5L')


#Import 10x data and convert each dataset to Seurat object (take from individual sample outputs, not cellranger's aggr output). Define sample with project= "samplename"
d53.data = Read10X(data.dir = './d53/outs/filtered_feature_bc_matrix/')
d53 = CreateSeuratObject(counts = d53.data, project = "d53", min.cells = 3, min.features = 200)
d59.data = Read10X(data.dir = './d59/outs/filtered_feature_bc_matrix/')
d59 = CreateSeuratObject(counts = d59.data, project = "d59", min.cells = 3, min.features = 200)
d74.data = Read10X(data.dir = './d74/outs/filtered_feature_bc_matrix/')
d74 = CreateSeuratObject(counts = d74.data, project = "d74", min.cells = 3, min.features = 200)
d78.data = Read10X(data.dir = './d78/outs/filtered_feature_bc_matrix/')
d78 = CreateSeuratObject(counts = d78.data, project = "d78", min.cells = 3, min.features = 200)
d113.data = Read10X(data.dir = './d113/outs/filtered_feature_bc_matrix/')
d113 = CreateSeuratObject(counts = d113.data, project = "d113", min.cells = 3, min.features = 200)
d132.data = Read10X(data.dir = './d132/outs/filtered_feature_bc_matrix/')
d132 = CreateSeuratObject(counts = d132.data, project = "d132", min.cells = 3, min.features = 200)
Hu37.data = Read10X(data.dir = './H37/outs/filtered_feature_bc_matrix/')
Hu37= CreateSeuratObject(counts = Hu37.data, project = "Hu37", min.cells = 3, min.features = 200)
Hu5.data = Read10X(data.dir = './Hu5/outs/filtered_feature_bc_matrix/')
Hu5= CreateSeuratObject(counts = Hu5.data, project = "Hu5", min.cells = 3, min.features = 200)
Hu7.data = Read10X(data.dir = './Hu7/outs/filtered_feature_bc_matrix/')
Hu7= CreateSeuratObject(counts = Hu7.data, project = "Hu7", min.cells = 3, min.features = 200)

# Merge into one single Seurat object
human=merge(d53, y=c(d59,d74,d78,d113,d132,Hu5,Hu7,Hu37))

#Validate the merge by checking number of cells per group
table(human$orig.ident)

##rm individual data
rm(d53,d59,d74,d78,d113,d132,Hu5,Hu7,Hu37,d53.data,d59.data,d74.data,d78.data,d113.data,d132.data,Hu5.data,Hu7.data,Hu37.data)

#Store mitochondrial percentage in the Seurat object metadata
human[["percent.mt"]] <- PercentageFeatureSet(human, pattern = "^MT-")

#Visualize QC metrics (will split by 'orig.ident')
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.qc.pdf", width=20)

#Add sample and condition information explicitly into the metadata (as oppsoed to storing in 'orig.ident') for future downstream analysis
#Add sample info (simply copying 'orig.ident')
human$sample <- plyr::mapvalues(
  x = human$orig.ident, 
  from = c('d53', 'd59', 'd74', 'd78', 'd113' , 'd132', 'Hu37', 'Hu5', 'Hu7'), 
  to = c('d53', 'd59', 'd74', 'd78', 'd113' , 'd132', 'Hu37', 'Hu5', 'Hu7')
)
#Add phenotype info (i.e. WT vs KO; allows for merging of different bioreps)
human$time <- plyr::mapvalues(
  x = human$orig.ident, 
  from = c('d53', 'd59', 'd74', 'd78', 'd113' , 'd132', 'Hu37', 'Hu5', 'Hu7'), 
  to = c('e50s', 'e50s', 'e70s', 'e70s','e100s', 'e100s', 'adult', 'adult', 'adult')
)

#Validate new metadata columns by checking that number of cells per sample/phenotype adds up
table(human$sample)
table(human$time)

#Run same QC metrics by new metadata columns to ensure it is the same as original QC metrics
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by='sample', ncol = 3, pt.size=0.1)
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by='time', ncol = 3, pt.size=0.1)

#Filter the data
human <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)

#Run standard Seurat workflow (Normalize, Find Variable Features, Scale, PCA) in order to check principal components between samples
human=NormalizeData(human)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)
human=ScaleData(human)
human = RunPCA(human)
DimPlot(human, reduction='pca', group.by='sample')
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.pca.pdf", width=20)
ElbowPlot(human, ndims=30)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.elbow_plot.pdf", width=20)

#If samples do not cluster together in PC space (which is the case here), then there is no need to run harmony (and likely no need to use the integrated method either); the merged analysis should do.
#Determine the number of dimensions to use in downstream analyses based on the point at which the Elbow Plot becomes flat (ok to be conservative)

#Save merged and filtered dataset as an R object to be able to re-load it for various future analyses without needing to perform the previous computations
saveRDS(human, file = "./seurat_analysis/all/human.rds")
##and load if needed
#human <- readRDS(file = "./seurat_analysis/all/human.rds")

#Run SCTransform, set var.to.regress to percent.mt
human_merged <- SCTransform(human, vars.to.regress = "percent.mt", verbose = FALSE)

##resolution 0.4 and 30 dims
#We can now run the standard Seurat workflow (PCA, UMAP, FindNeighbors, FindClusters). 
human_merged <- RunPCA(human_merged, verbose = FALSE)
human_merged <- RunUMAP(human_merged, dims = 1:30, verbose = FALSE)
human_merged <- FindNeighbors(human_merged, dims = 1:30, verbose=FALSE)
human_merged <- FindClusters(human_merged, resolution = 0.4)
#Plot UMAP
DimPlot(human_merged, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.merged.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human_merged, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.merged.UMAP.res0.4.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_merged$seurat_clusters, human_merged$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/all/human_scrnaseq_all_0620.merged.counts_cluster_sample.res0.4.csv')
##graph known markers
#Dot plot - the size of the dot = % of cells and color represents the average expression
DotPlot(human_merged, features = markers) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.merged.res0.4.known_markers.dotplot.pdf", width = 20)

#Find all markers that define each cluster -- not used
#human_merged.markers <- FindAllMarkers(human_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(human_merged.markers, file='./seurat_analysis/all/human_scrnaseq_all_0620.merged.markers.res0.4.csv')

##harmony analysis...
##load if needed
human <- readRDS(file = "./seurat_analysis/all/human.rds")

#Harmony needs RunPCA to have been performed, so make sure previous analysis steps have bene performed.
human_harmony <- RunHarmony(object = human, group.by.vars = 'sample')
DimPlot(human_harmony, reduction = 'harmony', group.by = 'sample')
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.harmony.harmony_plot.pdf", width = 20)

# Run standard Seurat workflow steps, but set reduction to "harmony" for UMAP and Neighbors
human_harmony <- RunUMAP(human_harmony, dims = 1:30, reduction = 'harmony')
human_harmony <- FindNeighbors(human_harmony, reduction = 'harmony', dims = 1:30)

##change resolutions - 0.4
human_harmony <- FindClusters(human_harmony, resolution = 0.4)
#Plot UMAP
DimPlot(human_harmony, reduction='umap', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.harmony.UMAP.res0.4.sample_split.pdf", width = 20)
DimPlot(human_harmony, reduction='umap', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.harmony.UMAP.res0.4.pdf", width=20)
#Identify the number of cells in each cluster between samples
counts_cluster_sample = table(human_harmony$seurat_clusters, human_harmony$sample)
write.csv(counts_cluster_sample, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.counts_cluster_sample.res0.4.csv')
counts_cluster_time = table(human_harmony$seurat_clusters, human_harmony$time)
write.csv(counts_cluster_time, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.counts_cluster_time.res0.4.csv')

#Dot plot - the size of the dot = % of cells and color represents the average expression
DotPlot(human_harmony, features = markers) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.harmony.res0.4.known_markers.dotplot.pdf", width = 20)
DotPlot(human_harmony, features = markers2) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.harmony.res0.4.new_markers.dotplot.pdf", width = 20)

#Find all markers that define each cluster
human_harmony.markers <- FindAllMarkers(human_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(human_harmony.markers, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.markers.res0.4.csv')

#Perform differential expression
human_harmony_time = human_harmony
Idents(human_harmony_time) <- 'time'
human_harmony.diffexp <- FindAllMarkers(human_harmony_time, min.pct=0.1, logfc.threshold = 0.25)
write.csv(human_harmony.diffexp, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.DE.res0.4.csv')

#Save object to avoid needing to re-run previous computations
saveRDS(human_harmony, file = "./seurat_analysis/all/human_harmony.rds")



##add cluster names and then extract cells and cluster ids

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

#Save object to avoid needing to re-run previous computations
saveRDS(human_harmony_subset, file = "./seurat_analysis/all/human_harmony_clusters_defined.rds")

##get gene expression by cell class i.e. cluster
human_harmony <- readRDS(file = "./seurat_analysis/all/human_harmony_clusters_defined.rds")
cluster.averages <- AverageExpression(human_harmony)
head(cluster.averages)
head(cluster.averages[["RNA"]])
write.csv(cluster.averages[["RNA"]], file='./seurat_analysis/all/human_scrnaseq_all_0620.cell_class_gene_matrix_112420.csv')

##get dotplots for marker genes 0121
human_harmony <- readRDS(file = "./seurat_analysis/all/human_harmony_clusters_defined.rds")

#Dot plot - the size of the dot = % of cells and color represents the average expression
##do 5x version with different genes
markers_a = c('PCK2', 'PCK1', 'PC', 'G6PC1', 'G6PC2', 'G6PC3', 'FBP1', 'FBP2', 'ME1', 'ME2', 'ME3')
DotPlot(human_harmony, features = markers_a) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.res0.4.markers_0121a.dotplot.pdf", width = 12, height = 6)
markers_b = c('SLC16A1', 'SLC16A7', 'SLC16A8', 'SLC16A3', 'SLC13A2', 'SLC13A3', 'SLC13A5', 'SLC25A10')
DotPlot(human_harmony, features = markers_b) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.res0.4.markers_0121b.dotplot.pdf", width = 12, height = 6)
markers_c = c('SLC1A3', 'GFAP', 'RHO', 'C1QA', 'GAD1', 'CABP5', 'OPN1LW', 'OPN1MW', 'RBPMS', 'THY1')
DotPlot(human_harmony, features = markers_c) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.res0.4.markers_0121c.dotplot.pdf", width = 12, height = 6)
markers_d = c('PGM1', 'PGM2', 'PGM3', 'PGM5', 'PYGM', 'PYGB', 'PYGL')
DotPlot(human_harmony, features = markers_d) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.res0.4.markers_0121d.dotplot.pdf", width = 12, height = 6)
markers_e = c('SLC25A31', 'SLC25A4', 'SLC25A5', 'SLC25A6', 'NDUFS1', 'NDUFS4', 'NDUFS8', 'SDHA', 'SDHB', 'SDHC', 'SDHD', 'COX4I1', 'COX4I2', 'COX7C', 'COX6B1', 'COX6B2', 'COX6A1', 'COX6A2', 'MT-ATP6', 'MT-ATP8', 'ATP5IF1', 'ATP5PO', 'SLC2A1')
DotPlot(human_harmony, features = markers_e) + RotatedAxis()
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.res0.4.markers_0121e.dotplot.pdf", width = 12, height = 6)






#>>>>>>>>>>>>>> egs
#Add celltype info for RNASeq analysis..... done 0820
##add cluster names and then extract cells and cluster ids
human_harmony <- readRDS(file = "./seurat_analysis/all/human_harmony.rds")
##make marked umap
new.cluster.ids <- c('Rods', 'Early Primary Progenitors', 'Ganglions', 'Amacrines', 'Ganglions', 'high mito rods', 'Amacrines', 'Cycling progenitors', 'Cones', 'Horizontals', 'Mullers', 'Cycling progenitors', 'Ganglions', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Amacrines', 'Ganglions', 'Early Primary Progenitors', 'Transitional 1', 'Ganglions', 'Astrocytes', 'Bipolars', 'bipolar/photoreceptor cell precursors', 'Microglia')
names(new.cluster.ids) <- levels(human_harmony)
human_harmony <- RenameIdents(human_harmony, new.cluster.ids)
DimPlot(human_harmony, reduction = "umap", pt.size = 0.1, label = TRUE)
dev.copy2pdf(file="./seurat_analysis/all/human_scrnaseq_all_0620.harmony.UMAP_celltypes.res0.4.pdf", width=20, height=20)

##get barcode info for each cluster
WhichCells(object = human_harmony, idents = c('Rods', 'Early Primary Progenitors', 'Ganglions', 'Amacrines', 'high mito rods', 'Cycling progenitors', 'Cones', 'Horizontals', 'Mullers', 'Bipolars', 'Transitional 1', 'Astrocytes', 'bipolar/photoreceptor cell precursors', 'Microglia'))
human_harmony$celltype <- plyr::mapvalues(
  x = human_harmony$seurat_clusters, 
  from = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26'), 
  to = c('Rods', 'Early Primary Progenitors', 'Ganglions', 'Amacrines', 'Ganglions', 'high mito rods', 'Amacrines', 'Cycling progenitors', 'Cones', 'Horizontals', 'Mullers', 'Cycling progenitors', 'Ganglions', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Bipolars', 'Amacrines', 'Ganglions', 'Early Primary Progenitors', 'Transitional 1', 'Ganglions', 'Astrocytes', 'Bipolars', 'bipolar/photoreceptor cell precursors', 'Microglia')
)
head(human_harmony$celltype)
colnames(human_harmony)
cell.info = human_harmony$celltype
write.csv(cell.info, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.res0.4.cellinfo.csv')

##get counts per sample/cell class
sample.celltype.info = table(human_harmony$celltype, human_harmony$orig.ident)
write.csv(sample.celltype.info, file='./seurat_analysis/all/human_scrnaseq_all_0620.harmony.res0.4.sample_cell_info.csv')

#Save object to avoid needing to re-run previous computations
#saveRDS(human_harmony, file = "./seurat_analysis/all/human_harmony.rds")
