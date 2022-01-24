# Set working directory to specific dir in active RSS
setwd('/active/cherry_t/Mactel_Mouse_SingleCell/5wk/mactel_mouse_scRNAseq')

# Load necessary libraries (dplyr v0.8.5; Seurat v3.1.1; patchwork v1.0.0; ggplot2 v3.3.0); make sure Seurat library has uwot installed (v0.1.4)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Adjust the maximum size of global objects (this may need to be increased later)
options(future.globals.maxSize = 8000 * 1024^2)

##get dotplots for marker genes 0122
markers_0122 = c('Slc16a1', 'Slc16a3', 'Slc16a7', 'Slc16a8')
mouse_harmony <- readRDS(file = "mouse.rds")
DotPlot(mouse_harmony, features = markers_0122) + RotatedAxis()
dev.copy2pdf(file="mactel_mouse_scRNAseq.SLC16_genes_0122.dotplot.pdf", width = 12, height = 6)

##change order
levels(mouse_harmony)
levels(mouse_harmony) <- c('Rods', 'Amacrine Cells', 'Bipolar Cells', 'Cones', 'Muller Glia', 'Horizontal Cells','RGCs', 'Microglia'  )
levels(mouse_harmony)
DotPlot(mouse_harmony, features = markers_0122) + RotatedAxis()
dev.copy2pdf(file="mactel_mouse_scRNAseq.SLC16_genes_0122.dotplot_chnaged.pdf", width = 12, height = 6)

