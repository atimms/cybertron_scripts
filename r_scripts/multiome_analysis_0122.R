# Load libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)


workingDir = "/active/cherry_t/chao_multiome_0122/at_seurat_0122";
setwd(workingDir);

## Generate Seurat objects for each sample

## Sample C1
# Data from folder various sample folders --> the 10x hdf5 file contains both data types.
C1.data <- Read10X_h5("/active/cherry_t/chao_multiome_0122/C1/outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- C1.data$`Gene Expression`
atac_counts <- C1.data$Peaks

# Create Seurat objects
C1 <- CreateSeuratObject(counts = rna_counts, project = 'C1')
C1[["percent.mt"]] <- PercentageFeatureSet(C1, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "/active/cherry_t/chao_multiome_0122/C1/outs/atac_fragments.tsv.gz",
  min.cells = 10,
  annotation = annotations
)
C1[["ATAC"]] <- chrom_assay

#Save rds objects prior to subsetting
saveRDS(C1, file = "/active/cherry_t/Emily/multieome_analysis/C1.rds")

## Look at quality control (QC) metrics
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)

VlnPlot(C1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
