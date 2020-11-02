library(Seurat)
##2 dataset to integrate
t21_seurat <- readRDS(file = '/active/millen_k/scRNA-seq/DS_CBL/RNA3-005-nova_data/fltrd_cerebellum_005_cds.RDS')
cbl_integrated_clean <- readRDS(file = '/active/millen_k/scRNA-seq/CBL_dev_dataset/R/seurat_objects/cbl_integrated_cleanCC_200518.rds')


DefaultAssay(cbl_integrated_clean) <- 'RNA'

ref_list <- c(cbl_integrated_clean, t21_seurat)
anchors <- FindIntegrationAnchors(object.list = ref_list, dims = 1:50)
saveRDS(anchors, file = '/active/millen_k/scRNA-seq/DS_CBL/t21_cbl_dev_anchors.rds')
integration <- IntegrateData(anchorset = anchors, dims = 1:50)
DefaultAssay(integration) = 'integrated'
integration <- ScaleData(integration)
integration <- RunPCA(integration, npcs = 50)
integration <- RunUMAP(integration, dims = 1:50)
integration <- FindNeighbors(integration, dims = 1:50)
integration <- FindClusters(integration, resolution = 0.5)
saveRDS(integration, file = '/active/millen_k/scRNA-seq/DS_CBL/t21_cbl_dev_integration.rds')