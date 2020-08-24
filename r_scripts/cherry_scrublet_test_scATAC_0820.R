##install signac
# Install bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
# To automatically install Bioconductor dependencies
#setRepositories(ind=1:2)
##install signac
#install.packages("Signac")
##Installing genome assembly and gene annotation 
#BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomeInfoDb)
library(ggplot2)
library(harmony)
library(future)
library(GenomicRanges)
set.seed(1234)
##use parallelization
plan("multisession", workers = 4) #use 4 cores
options(future.globals.maxSize = 80000 * 1024^2) # for 80 Gb RAM

##Set working directory and specific dir in active RSS where the data is sitting
setwd('/home/atimms/ngs_data/misc/cherry_scrublet_test_0820')

##Hu5
##load data
counts <- Read10X_h5("/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu5/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(file = "/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu5/outs/singlecell.csv",
  header = TRUE, row.names = 1)
##issue adding genome, 502 error
Hu5_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), genome = "hg38",
  fragments = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu5/outs/fragments.tsv.gz',
  min.cells = 1)
Hu5_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"),
                                  fragments = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu5/outs/fragments.tsv.gz',
                                  min.cells = 1)
Hu5 <- CreateSeuratObject(counts = Hu5_assay, assay = 'peaks',
  project = 'ATAC',meta.data = metadata)
##add gene annotations
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# add the gene information to the object
Annotation(Hu5) <- annotations

##QC metrics
Hu5 <- NucleosomeSignal(object = Hu5)
Hu5$nucleosome_group <- ifelse(Hu5$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Hu5, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
Hu5 <- TSSEnrichment(Hu5, fast = FALSE)
Hu5$high.tss <- ifelse(Hu5$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(Hu5, group.by = 'high.tss') + NoLegend()
Hu5$pct_reads_in_peaks <- Hu5$peak_region_fragments / Hu5$passed_filters * 100
Hu5$blacklist_ratio <- Hu5$blacklist_region_fragments / Hu5$peak_region_fragments

VlnPlot(object = Hu5, features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1, ncol = 5)

DepthCor(Hu5)


