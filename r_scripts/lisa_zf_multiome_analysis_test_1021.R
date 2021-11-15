# Seurat multiomic of zebrafish samples

##install zf annotation
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("BSgenome.Drerio.UCSC.danRer11")

# Load libraries
library(Seurat)
library(Signac)
library(harmony)
library(BSgenome.Drerio.UCSC.danRer11)
#library(dplyr)
library(ggplot2)

##analysis from Eric (uses the 2nd link):
# https://satijalab.org/signac/articles/pbmc_multiomic.html
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1

# Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 200 * 1024^3)

# move to working dir
workingDir = "/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis";
setwd(workingDir);

# Generate Seurat objects for each timepoint

## WT_1

# the 10x hdf5 file contains both data types.
WT_1.data <- Read10X_h5("/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/WT_1/outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- WT_1.data$`Gene Expression`
atac_counts <- WT_1.data$Peaks

# Create Seurat object
WT_1 <- CreateSeuratObject(counts = rna_counts, project = 'WT_1')
#WT_1[["percent.mt"]] <- PercentageFeatureSet(WT_1, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = BSgenome.Drerio.UCSC.danRer11)
annotations <- GRanges(seqinfo(BSgenome.Drerio.UCSC.danRer11))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "danRer11"

frag.file <- "/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/WT_1/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'danRer11',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
WT_1[["ATAC"]] <- chrom_assay

#Save rds objects prior to subsetting
saveRDS(WT_1, file = "WT_1.rds")

# Look at QC metrics
VlnPlot(WT_1, features = c("nCount_ATAC", "nCount_RNA"), ncol = 2,
        log = TRUE, pt.size = 0) + NoLegend()
dev.copy2pdf(file='WT_1.QC_metrics.pdf', width = 6, height = 4)

# Subset based on QC metrics
WT_1 <- subset(
  x = WT_1,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 2e2 &
    nCount_RNA < 1e4 &
    nCount_RNA > 2e2
)

## WT_2

# the 10x hdf5 file contains both data types.
WT_2.data <- Read10X_h5("/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/WT_2/outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- WT_2.data$`Gene Expression`
atac_counts <- WT_2.data$Peaks

# Create Seurat object
WT_2 <- CreateSeuratObject(counts = rna_counts, project = 'WT_2')
#WT_2[["percent.mt"]] <- PercentageFeatureSet(WT_2, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = BSgenome.Drerio.UCSC.danRer11)
annotations <- GRanges(seqinfo(BSgenome.Drerio.UCSC.danRer11))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "danRer11"

frag.file <- "/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/WT_2/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'danRer11',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
WT_2[["ATAC"]] <- chrom_assay

#Save rds objects prior to subsetting
saveRDS(WT_2, file = "WT_2.rds")

# Look at QC metrics
VlnPlot(WT_2, features = c("nCount_ATAC", "nCount_RNA"), ncol = 2,
        log = TRUE, pt.size = 0) + NoLegend()
dev.copy2pdf(file='WT_2.QC_metrics.pdf', width = 6, height = 4)

# Subset based on QC metrics
WT_2 <- subset(
  x = WT_2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 2e2 &
    nCount_RNA < 1e4 &
    nCount_RNA > 2e2
)

## dmd_1

# the 10x hdf5 file contains both data types.
dmd_1.data <- Read10X_h5("/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/dmd_1/outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- dmd_1.data$`Gene Expression`
atac_counts <- dmd_1.data$Peaks

# Create Seurat object
dmd_1 <- CreateSeuratObject(counts = rna_counts, project = 'dmd_1')
#dmd_1[["percent.mt"]] <- PercentageFeatureSet(dmd_1, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = BSgenome.Drerio.UCSC.danRer11)
annotations <- GRanges(seqinfo(BSgenome.Drerio.UCSC.danRer11))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "danRer11"

frag.file <- "/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/dmd_1/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'danRer11',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
dmd_1[["ATAC"]] <- chrom_assay

#Save rds objects prior to subsetting
saveRDS(dmd_1, file = "dmd_1.rds")

# Look at QC metrics
VlnPlot(dmd_1, features = c("nCount_ATAC", "nCount_RNA"), ncol = 2,
        log = TRUE, pt.size = 0) + NoLegend()
dev.copy2pdf(file='dmd_1.QC_metrics.pdf', width = 6, height = 4)

# Subset based on QC metrics
dmd_1 <- subset(
  x = dmd_1,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 2e2 &
    nCount_RNA < 1e4 &
    nCount_RNA > 2e2
)

## dmd_2

# the 10x hdf5 file contains both data types.
dmd_2.data <- Read10X_h5("/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/dmd_2/outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- dmd_2.data$`Gene Expression`
atac_counts <- dmd_2.data$Peaks

# Create Seurat object
dmd_2 <- CreateSeuratObject(counts = rna_counts, project = 'dmd_2')
#dmd_2[["percent.mt"]] <- PercentageFeatureSet(dmd_2, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = BSgenome.Drerio.UCSC.danRer11)
annotations <- GRanges(seqinfo(BSgenome.Drerio.UCSC.danRer11))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "danRer11"

frag.file <- "/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/modified_analysis/dmd_2/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'danRer11',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
dmd_2[["ATAC"]] <- chrom_assay

#Save rds objects prior to subsetting
saveRDS(dmd_2, file = "dmd_2.rds")

# Look at QC metrics
VlnPlot(dmd_2, features = c("nCount_ATAC", "nCount_RNA"), ncol = 2,
        log = TRUE, pt.size = 0) + NoLegend()
dev.copy2pdf(file='dmd_2.QC_metrics.pdf', width = 6, height = 4)

# Subset based on QC metrics
dmd_2 <- subset(
  x = dmd_2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 2e2 &
    nCount_RNA < 1e4 &
    nCount_RNA > 2e2
)

# Merge individual seurat objects
comb = merge(WT_1, 
             y = c(WT_2, dmd_1, dmd_2), 
             add.cell.ids = c('WT_1', 'WT_2', 'dmd_1', 'dmd_2')
)

#Add sample info (simply copying 'orig.ident')
comb$sample <- plyr::mapvalues(
  x = comb$orig.ident, 
  from = c('WT_1', 'WT_2', 'dmd_1', 'dmd_2'), 
  to = c('WT_1', 'WT_2', 'dmd_1', 'dmd_2'),
)
#Add genotype info
comb$genotype <- plyr::mapvalues(
  x = comb$orig.ident, 
  from = c('WT_1', 'WT_2', 'dmd_1', 'dmd_2'), 
  to = c('WT', 'WT', 'dmd', 'dmd'),
)

##remove individual sets
rm(dmd_1,dmd_1.data,dmd_2,dmd_2.data,WT_1,WT_1.data,WT_2,WT_2.data)

# Save object and load if needed
saveRDS(comb, file = "combined.rds")
comb <- readRDS(file ="combined.rds")


## RNA analysis
# Run SCTransform and PCA, then check samples in PCA space
DefaultAssay(comb) <- "RNA"
comb <- SCTransform(comb, verbose = FALSE)
comb <- RunPCA(comb, verbose = FALSE)
comb <- RunUMAP(comb, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
##view RNA umap and PCA plots
DimPlot(comb, reduction='umap.rna', split.by='sample', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="umap_rna.split_by_sample.pdf", width=20)
DimPlot(comb, reduction='umap.rna', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="umap_rna.by_sample.pdf", width=20)
DimPlot(comb, reduction='pca', group.by='sample')
dev.copy2pdf(file="pca_rna.pdf", width=20)
#Samples cluster together in PCA space, so probably no need for harmony on RNA

# Test harmony integration (RNA only)
#wk12_rna <- SCTransform(wk12, verbose = FALSE) %>% RunHarmony(wk12, group.by.vars = 'sample') %>% RunUMAP(dims = 1:50, reduction = 'harmony', reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
## Cells cluster completely separately in ATAC without harmony; so run harmony after SVD
DefaultAssay(comb) <- "ATAC"
comb <- RunTFIDF(comb)
comb <- FindTopFeatures(comb, min.cutoff = 'q0')
comb <- RunSVD(comb)
##checked umap and samples seperate so need to run harmony, which 
comb <- RunHarmony(comb, group.by.vars = 'sample', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
comb <- RunUMAP(comb, reduction = 'harmony', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(comb, reduction='umap.atac', pt.size=0.1, label = TRUE)
dev.copy2pdf(file="umap_atac.harmony.by_sample.pdf", width=20)

# WNN analysis
comb <- FindMultiModalNeighbors(comb, reduction.list = list("pca", "harmony"), dims.list = list(1:50, 2:50))
comb <- RunUMAP(comb, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
comb <- FindClusters(comb, graph.name = "wsnn", algorithm = 3, resolution = 0.4, verbose = FALSE)

# Visualize clustering based on RNA, ATAC, or WNN
p1 <- DimPlot(comb, reduction = "umap.rna", label = TRUE, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(comb, reduction = "umap.atac", label = TRUE, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(comb, reduction = "wnn.umap", label = TRUE, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.copy2pdf(file="umap.rna_atac_wnn.pdf", width=20)

# Save object and load if needed
saveRDS(comb, file = "combined.rds")
comb <- readRDS(file ="combined.rds")



# Test switching to RNA assay and normalizing (do not save object until you hear back from Seurat folks)
DefaultAssay(comb) = 'RNA'
comb = NormalizeData(comb)
##by cluster
comb.markers = FindAllMarkers(comb, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(comb.markers, file = 'RNA_cluster_markers.csv')
##by genotype
Idents(comb) = 'genotype'
comb.de = FindAllMarkers(comb, logfc.threshold = 0)
write.csv(comb.de, file = 'RNA_genotype_markers.csv')

# Run marker analysis differential expression on ATAC assay
DefaultAssay(comb) = 'ATAC'
Idents(comb) = 'seurat_clusters'
##by cluster
comb_atac.markers = FindAllMarkers(comb, only.pos = T, logfc.threshold = 0.25)
write.csv(comb_atac.markers, file = 'ATAC_cluster_markers.csv')
Idents(comb) = 'genotype'
comb_atac.de = FindAllMarkers(comb, logfc.threshold = 0)
write.csv(comb_atac.de, file = 'ATAC_genotype_markers.csv')


