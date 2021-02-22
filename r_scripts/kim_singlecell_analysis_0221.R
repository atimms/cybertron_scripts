workingDir = "~/singlecell";
setwd(workingDir);
getwd()

library(Seurat)
#load cleaned cbl dev dataset
cbl_dev_clean <- readRDS(file = '/active/millen_k/scRNA-seq/CBL_dev_dataset/R/seurat_objects/cbl_integrated_cleanCC_200518.rds')
#rename to use Zach's code with minimal other changes
cbl_dev_rna <- cbl_dev_clean
DefaultAssay(cbl_dev_rna) <- 'RNA'
cbl_dev_rna <- NormalizeData(cbl_dev_rna)
cbl_dev_rna <- SetIdent(cbl_dev_rna, value = 'figure_clusters')
#select only 100 cells per cell type
cbl_dev_rna_subset <- subset(cbl_dev_rna, max.cells.per.ident = 100)
avg_exp <- AverageExpression(cbl_dev_rna_subset)
write.csv(avg_exp$RNA, file="/active/millen_k/scRNA-seq/CBL_dev_dataset/R/cbl_dev_clean_100cells_ave.csv")

## Jun's brain data
jun_brain<-readRDS(file = '/active/millen_k/scRNA-seq/Juns_data/cerebrum_seurat_small.rds')
VlnPlot(jun_brain, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
jun_brain <- FindVariableFeatures(jun_brain, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(jun_brain), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(jun_brain)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# cluster cells
jun_brain <- FindNeighbors(jun_brain, dims = 1:10)
jun_brain <- FindClusters(jun_brain, resolution = 0.5)
# visualize
# install UMAP # reticulate::py_install(packages ='umap-learn')
# jun_brain <- RunUMAP(jun_brain, dims = 1:10) # didn't work to install
jun_brain <- RunTSNE(jun_brain, dims = 1:10)
DimPlot(jun_brain, reduction = "tsne")
TSNEPlot(jun_brain,do.label = TRUE,label.pt.size = 0.5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
jun_brain.markers <- FindAllMarkers(jun_brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
jun_brain.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
VlnPlot(jun_brain, features = c("MAST4", "PDGFRA"))

#try visium data
#output cells per sample
cbl_dev_clean_cells_per_sample<-table(cbl_dev_clean@active.ident, cbl_dev_clean@meta.data$sample_id_aligned)
write.csv(cbl_dev_clean_cells_per_sample, file = '/active/millen_k/kim_R/singlecell/cbl_dev_clean_cells_per_sample.csv')

#output metadata
data_to_write_out <- as.data.frame(as.matrix(cbl_dev_clean@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "/active/millen_k/kim_R/singlecell/cbl_dev_clean_meta.csv")

#cell type enrichment 2021-02-11
#using template code from Tom Nowakowski/Maddy
#try using CPM but that file is HUGE and takes forever just to load into R
#CPM file of 70K data is: /active/millen_k/kim_R/singlecell/cerebellum/cbl_integrated_cleanCC_200518_CPM.csv
cbl_data<-read.table('/active/millen_k/kim_R/singlecell/cerebellum/exprMatrix.tsv', row.names=1, header=T)
cbl_cluster<-read.table('/active/millen_k/kim_R/singlecell/cerebellum/meta_celltype_labels.txt', sep="\t", row.names=1, header=T)
cbl_tsne<-read.table('/active/millen_k/kim_R/singlecell/cerebellum/Seurat_tSNE.coords.tsv', header=T, row.names=1)

data<-log10(t(cbl_data[,row.names(cbl_tsne)])+1)
cluster<-cbl_cluster[row.names(cbl_tsne),1]
gene<-c("AHDC1", "AP1S2", "ARMC9", "AUTS2", "B3GALTL", "B4GALT1", "BCL11A", "BRAF", "BUB1B", "CASK", "CCDC22", "CDC42", "CDKN1C", "CHD7", "DDX3X", "DKC1", "DPYSL5", "EBF2", "EBF3", "FGFR1", "FOXC1", "FOXP1", "FZD3", "HRAS", "INSIG1", "KIAA0196", "KIF4A", "L1CAM", "LAMC1", "MACF1", "MAST1", "NID1", "OPHN1", "PDGFRB", "PITX2", "PMM2", "PPP1CB", "PTF1A", "PUS3", "RARS2", "SETD2", "SHANK3", "SPTAN1", "STXBP1", "TCF4", "TMLHE", "TUBA1A", "TUBB2A", "USP9X", "WDR37", "WNT1", "ZIC1", "ZIC4", "ZNF292")
types<-unique(cluster)
keys<-cluster
id<-0;

id2<-0;

for(i in 1: length(types)){
  a<-data[which(keys==types[i]),gene]
  b<-data[-which(keys==types[i]),gene]
  id[i]<-wilcox.test(a,b)[3]
  id2[i]<-log((length(which(a>0))/length(a))/(length(which(b>0))/length(b)))
}

G<-data.frame(cbind(types,unlist(id),id2))
G[,1]<-types
G[,2]<-p.adjust(G[,2], method = "bonferroni", n = length(types))

D<-data.frame(cbind(keys,data[,gene]));D[,1]<-keys
write.csv(G,"/active/millen_k/kim_R/singlecell/cerebellum/summary_BCL11A.csv")
#above didn't seem to work

#try escape?
#https://ncborcherding.github.io/vignettes/escape_vignette.html

#genelists to test for enrichment for paper revision
#my_gene_sets <- list(gene_set_1 = c("geneA", "geneB", "geneC"))
#6 lists to test below

list(CBLM_geneset = c(“AHDC1”, “AP1S2”, “ARMC9”, “AUTS2”, “B3GALTL”, “B4GALT1”, “BCL11A”, “BRAF”, “BUB1B”, “CASK”, “CCDC22”, “CDC42”, “CDKN1C”, “CHD7”, “DDX3X”, “DKC1”, “DPYSL5”, “EBF2”, “EBF3”, “FGFR1”, “FOXC1”, “FOXP1”, “FZD3”, “HRAS”, “INSIG1”, “KIAA0196”, “KIF4A”, “L1CAM”, “LAMC1”, “MACF1”, “MAST1”, “NID1”, “OPHN1”, “PDGFRB”, “PITX2”, “PMM2”, “PPP1CB”, “PTF1A”, “PUS3”, “RARS2”, “SETD2”, “SHANK3”, “SPTAN1”, “STXBP1”, “TCF4”, “TMLHE”, “TUBA1A”, “TUBB2A”, “USP9X”, “WDR37”, “WNT1”, “ZIC1”, “ZIC4”, “ZNF292”)) 
list(JS_geneset = c(“AHI1”, “ARL13B”, “ARMC9”, “B9D1”, “B9D2”, “C2CD3”, “C5orf42”, “CC2D2A”, “CEP104”, “CEP120”, “CEP164”, “CEP290”, “CEP41”, “CEP83”, “CLUAP1”, “CSPP1”, “CBY1”, “HYLS1”, “IFT74”, “IFT172”, “INPP5E”, “KIAA0556”, “KIAA0586”, “KIAA0753”, “KIF7”, “MKS1”, “NPHP1”, “NPHP4”, “OFD1”, “PDE6D”, “PIBF1”, “RPGRIP1L”, “TCTN1”, “TCTN2”, “TCTN3”, “TMEM107”, “TMEM138”, “TMEM216”, “TMEM231”, “TMEM237”, “TMEM67”, “TOGARAM1”))
list(ID_geneset = c(“ADNP”, “AHDC1”, “ANKRD11”, “ARHGEF9”, “ARID1B”, “ARIH1”, “ASXL1”, “ASXL3”, “AUTS2”, “BCL11A”, “BRAF”, “BRPF1”, “C9orf142”, “CA5B”, “CAMK2A”, “CASK”, “CBL”, “CDK13”, “CDKL5”, “CHAMP1”, “CHD2”, “CHD3”, “CHD4”, “CHMP2A”, “CIDEC”, “CLTC”, “CNKSR2”, “CNOT3”, “COL23A1”, “COL4A3BP”, “CREBBP”, “CSNK2A1”, “CSNK2B”, “CTCF”, “CTNNB1”, “CUL3”, “CYP27C1”, “DCX”, “DDX3X”, “DNM1”, “DNMT3A”, “DYRK1A”, “EBF3”, “EEF1A2”, “EFTUD2”, “EHMT1”, “EIF4A2”, “EIF4EBP1”, “EP300”, “FAM104A”, “FAM200A”, “FAM84A”, “FAM98A”, “FOSL2”, “FOXG1”, “FOXP1”, “FOXP2”, “GABRB2”, “GATAD2B”, “GFOD2”, “GNAI1”, “GNAO1”, “GOLPH3”, “GRB14”, “GRIN2A”, “GRIN2B”, “GUCA2A”, “HDAC8”, “HECW2”, “HIST1H1E”, “HIST1H2AC”, “HK1”, “HNRNPK”, “HNRNPU”, “IQSEC2”, “KANSL1”, “KAT6A”, “KAT6B”, “KCNB1”, “KCNH1”, “KCNK3”, “KCNQ2”, “KCNQ3”, “KDM6A”, “KIAA2022”, “KIF1A”, “KMT2A”, “MAP2K1”, “MAP4K4”, “MBD5”, “MECP2”, “MED13L”, “MEF2C”, “MEIS2”, “MRPL41”, “MSI1”, “MSL3”, “MTF1”, “MYT1L”, “NAA10”, “NALCN”, “NFIX”, “NSD1”, “PACS1”, “PAQR6”, “PDHA1”, “PDX1”, “PLAC8L1”, “PLEKHB2”, “POGZ”, “POU3F3”, “PPM1D”, “PPP1CB”, “PPP1R12A”, “PPP2R1A”, “PPP2R5D”, “PRKAR1A”, “PRKG1”, “PRPF40A”, “PTEN”, “PTPN11”, “PUF60”, “PUM2”, “PURA”, “QRICH1”, “RAB11A”, “RABGAP1L”, “RAC1”, “SATB2”, “SCN2A”, “SCN8A”, “SET”, “SETBP1”, “SETD5”, “SETDB2”, “SIAH1”, “SIN3A”, “SIX3”, “SLC12A2”, “SLC22A23”, “SLC2A1”, “SLC6A1”, “SMAD4”, “SMARCA2”, “SMC1A”, “SMPD2”, “SNAP25”, “SNX11”, “SOX11”, “SOX4”, “SOX5”, “STXBP1”, “SVIP”, “SYNCRIP”, “SYNGAP1”, “TAB2”, “TAOK1”, “TBL1XR1”, “TCF12”, “TCF20”, “TCF4”, “TEX12”, “TGFB2”, “TLK2”, “TMPRSS12”, “TNNI2”, “TNPO3”, “TRAF7”, “TRIM72”, “TRIP12”, “TUBA1A”, “U2AF2”, “UPF3B”, “USP50”, “USP9X”, “VAMP2”, “VEZF1”, “WAC”, “WDR26”, “WDR45”, “WHSC1”, “ZBTB10”, “ZBTB18”, “ZBTB20”, “ZC4H2”, “ZMYND11”))
list(ASD_geneset = c(“ACHE”, “ADCY3”, “ADNP”, “AGAP2”, “AKAP9”, “ANK2”, “APH1A”, “ARID1B”, “ASH1L”, “ASXL3”, “BCL11A”, “BTRC”, “C16orf13”, “CACNA2D3”, “CAPN12”, “CCIN”, “CCSER1”, “CHD2”, “CHD8”, “CIC”, “CLASP1”, “CMPK2”, “CNOT3”, “CTCF”, “CTTNBP2”, “CUL3”, “DDX3X”, “DIP2C”, “DNMT3A”, “DPP3”, “DSCAM”, “DYNC1H1”, “DYRK1A”, “ERBIN”, “ETFB”, “FAM47A”, “FAM98C”, “FBXO11”, “FOXP1”, “GABRB3”, “GIGYF1”, “GIMAP8”, “GRIA1”, “GRIN2B”, “ILF2”, “INTS6”, “KATNAL2”, “KDM5B”, “KDM6B”, “KIAA2022”, “KMT2C”, “KMT2E”, “KMT5B”, “MED13”, “MFRP”, “MLANA”, “MYO5A”, “MYT1L”, “NAA15”, “NCKAP1”, “NR3C2”, “NRXN1”, “NUAK1”, “P2RX5”, “PAX5”, “PCDH11X”, “PCM1”, “PHF2”, “PHF3”, “PLCD4”, “POGZ”, “PRB4”, “PRKAR1B”, “PTEN”, “PTK7”, “PTMS”, “PYHIN1”, “RANBP17”, “RAPGEF4”, “RIMS1”, “S100G”, “SCN2A”, “SETD5”, “SHANK2”, “SHANK3”, “SLC6A1”, “SMARCC2”, “SMURF1”, “SPAST”, “SRSF11”, “SUV420H1”, “SYNGAP1”, “TAF6”, “TBL1XR1”, “TBR1”, “TCF7L2”, “TLK2”, “TMEM39B”, “TNRC6B”, “TRIP12”, “TSPAN4”, “UBN2”, “UIMC1”, “USP45”, “WAC”, “WDFY3”, “ZC3H11A”, “ZNF559”))
list(SCA_geneset = c(“ADH1C”, “AFG3L2”, “ATN1”, “ATXN1”, “ATXN10”, “ATXN2”, “ATXN3”, “ATXN7”, “ATXN8”, “ATXN8OS”, “BEAN1”, “CACNA1A”, “CAPN1”, “DAB1”, “ELOVL4”, “ELOVL5”, “EP300”, “FAT2”, “FGF14”, “GBA”, “GJB1”, “GLUD2”, “ITPR1”, “KCNC3”, “KCND3”, “KCNJ10”, “KIF26B”, “MAPT”, “MYADM”, “NEFL”, “NOP56”, “PDYN”, “PLD3”, “PPP2R2B”, “PRKCG”, “SLC4A1”, “SPTBN2”, “SYNE1”, “SYT14”, “TBP”, “TGM6”, “TMEM240”, “TTC19”, “TTPA”))
list(ALZ_geneset = c(“MARCH10”, “ABCA7”, “AC099552.4”, “ACE”, “ACP2”, “ADAMTS20”, “AGFG2”, “ANKHD1”, “ANKHD1-EIF4EBP3”, “AP4M1”, “APBB3”, “APOE”, “APP”, “BCAM”, “C5orf60”, “C7orf43”, “CASS4”, “CBLC”, “CBX3”, “CCDC25”, “CCDC81”, “CD33”, “CEND1”, “CLCN1”, “CNPY4”, “CR1”, “CR1L”, “CR2”, “CYP3A43”, “DHX33”, “DTNBP1”, “EFCAB4A”, “EIF4EBP3”, “ELMO1”, “EPHX2”, “ESCO2”, “FAM209B”, “FNBP4”, “GAL3ST4”, “GAS2L2”, “GP1BA”, “GPD2”, “GRIN3B”, “HLA-DPA1”, “HSF5”, “IGHG3”, “IGIP”, “IMPA2”, “INCA1”, “INPP5D”, “KCNH6”, “KIF1C”, “KLHL40”, “LDB3”, “LILRB1”, “LPO”, “MGAT4B”, “MINK1”, “MS4A6A”, “MS4A7”, “MUC6”, “NGEF”, “NKG7”, “NME8”, “NRG2”, “NSF”, “NUP88”, “NYAP1”, “OPN5”, “OPRL1”, “OR2AE1”, “PBK”, “PCDHA1”, “PCDHA10”, “PCDHA12”, “PCDHA3”, “PCDHA5”, “PCDHA7”, “PCDHA9”, “PCOLCE”, “PILRA”, “PIP4K2C”, “PNPLA2”, “POLR2E”, “PSEN1”, “PSEN2”, “PTK2B”, “RHBDD1”, “RIN3”, “RPL9”, “RPS16”, “SAG”, “SCARA5”, “SIGLEC7”, “SIGLEC9”, “SIRPB1”, “SLC44A4”, “SLC4A9”, “SLC52A1”, “SNX1”, “SORL1”, “SPAG7”, “SPPL2A”, “SQSTM1”, “STAG3”, “STX18”, “STYX”, “SYTL2”, “TAMM41”, “TANC2”, “TBC1D9B”, “TFR2”, “TREM2”, “TREML4”, “USP6NL”, “USP8”, “ZCWPW1”, “ZNF594”, “ZNF609”, “ZNF655”))
                              