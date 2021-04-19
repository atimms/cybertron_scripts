##just use common library, so check which one to use and then set that parameter
#.libPaths()
#.libPaths( .libPaths()[2] )
##load libarys
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

workingDir = "/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321";
setwd(workingDir);

###cunn_flexcell_rnaseq_0321 --- all samples
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Gene, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0321.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.25_var_gene_clustering.pdf', width = 10, height = 10)


###cunn_flexcell_rnaseq_0321_cell_line --- all samples minus bone and suture
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_cell_line.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_cell_line.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 6, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321_cell_line.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321_cell_line.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Gene, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321_cell_line.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321_cell_line.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0321_cell_line.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321_cell_line.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321_cell_line.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.25_var_gene_clustering.pdf', width = 10, height = 10)

##clustering/heatmap from lists
genes = c('ABL1', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'AJUBA', 'APOD', 'ARHGAP6', 'ARHGEF7', 'BCAS3', 'BCL2', 'BCR', 'CAMSAP3', 'CLASP1', 'CLASP2', 'COL16A1', 'CORO1C', 'CORO2B', 'CTTN', 'DAPK3', 'DLC1', 'DMTN', 'DUSP22', 'DUSP3', 'EFNA5', 'EPB41L5', 'EPHA3', 'FAM107A', 'FERMT2', 'FMN1', 'GPM6B', 'GREM1', 'HRG', 'IQGAP1', 'ITGA2', 'ITGB1BP1', 'KDR', 'LDB1', 'LIMCH1', 'LIMS1', 'LRP1', 'MACF1', 'MAP4K4', 'MMP14', 'MYOC', 'NRP1', 'PDPK1', 'PEAK1', 'PHLDB2', 'PIP5K1A', 'PPM1F', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'RAC1', 'RCC2', 'RHOA', 'RHOD', 'ROCK1', 'ROCK2', 'S100A10', 'SDC4', 'SFRP1', 'SLC9A1', 'SLK', 'SMAD3', 'SORBS1', 'SRC', 'TAOK2', 'TEK', 'TESK2', 'THBS1', 'THSD1', 'THY1', 'TRIP6', 'TSC1', 'VCL', 'VEGFA', 'WDPCP', 'WHAMM', 'WNT4')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.focal_adhesion.heatmap.pdf', width = 9, height = 6)

genes = c('ABCB1', 'ACE', 'AGO3', 'ARIH2', 'ATXN1L', 'CCNE1', 'CD34', 'CITED1', 'CTC1', 'EIF2AK2', 'EPCAM', 'ETV6', 'FBLN1', 'FERMT1', 'FERMT2', 'FGF2', 'FZD1', 'GBA', 'GJA1', 'HMGA2', 'HMGB2', 'HNRNPU', 'KAT7', 'KDF1', 'KDM1A', 'KITLG', 'KRT18', 'LTBP3', 'MECOM', 'MIR16-1', 'MIR221', 'MIR222', 'MIR29B1', 'N4BP2L2', 'NANOG', 'NES', 'NF2', 'NR2E1', 'OVOL1', 'OVOL2', 'PDCD2', 'PDGFRA', 'PIM1', 'PTPRC', 'REST', 'RNF43', 'RUNX1', 'SFRP2', 'SIX2', 'SOX11', 'SOX17', 'SOX18', 'SOX5', 'SOX6', 'SOX9', 'TBX3', 'TERT', 'THPO', 'TRIM71', 'VEGFC', 'WNT1', 'WNT10B', 'WNT2B', 'WNT3', 'WNT5A', 'WNT7B', 'YAP1', 'YJEFN3', 'YTHDF2', 'ZFP36L1', 'ZNRF3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.msc_stem.heatmap.pdf', width = 9, height = 6)

genes = c('ABI2', 'ACTB', 'ADAM10', 'ADAM15', 'ADD1', 'AFDN', 'AHI1', 'AJAP1', 'AJM1', 'AJUBA', 'ALOX15B', 'ANG', 'ANXA1', 'ANXA2', 'APC', 'ARHGAP24', 'ARVCF', 'BAIAP2', 'BAIAP2L1', 'BMP6', 'BMPR2', 'CADM1', 'CADM2', 'CADM3', 'CAMSAP3', 'CCDC85A', 'CCDC85B', 'CCDC85C', 'CD99L2', 'CDC42', 'CDC42EP1', 'CDC42EP4', 'CDCA3', 'CDH1', 'CDH10', 'CDH11', 'CDH12', 'CDH13', 'CDH15', 'CDH17', 'CDH18', 'CDH19', 'CDH2', 'CDH20', 'CDH22', 'CDH24', 'CDH3', 'CDH4', 'CDH5', 'CDH6', 'CDH7', 'CDH8', 'CDH9', 'CDHR3', 'CEACAM1', 'CNN3', 'CRB1', 'CSK', 'CTNNA1', 'CTNNA2', 'CTNNA3', 'CTNNB1', 'CTNND1', 'CTNND2', 'CXADR', 'CYTH1', 'CYTH2', 'CYTH3', 'DAG1', 'DCHS1', 'DDX6', 'DLG1', 'DLG5', 'DLL1', 'DSC2', 'DSP', 'EFNB2', 'EIF4G2', 'EPB41L5', 'EPHA4', 'ESAM', 'FAT2', 'FERMT2', 'FLOT1', 'FLOT2', 'FMN1', 'FRMD4A', 'FRMD4B', 'FRMD5', 'FRS2', 'HIPK1', 'HMCN1', 'IGSF21', 'INAVA', 'ITGA6', 'JAG1', 'JAM3', 'JCAD', 'JUP', 'KIFC3', 'KLHL24', 'KRT18', 'LDB3', 'LIMD1', 'LIN7A', 'LIN7B', 'LIN7C', 'LYN', 'MAGI1', 'MPP4', 'MPP5', 'MPP7', 'MYH9', 'MYO1E', 'NDRG1', 'NECTIN1', 'NECTIN2', 'NECTIN3', 'NECTIN4', 'NEXN', 'NF2', 'NIBAN2', 'NOTCH1', 'NPHP1', 'NUMB', 'NUMBL', 'OXTR', 'PAK2', 'PAK4', 'PARD3', 'PARD3B', 'PARK7', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDLIM5', 'PDLIM7', 'PDZD11', 'PGM5', 'PIP5K1C', 'PKP1', 'PKP2', 'PKP3', 'PKP4', 'PLEKHA7', 'PLPP3', 'POF1B', 'PPP1CA', 'PPP1R9B', 'PRICKLE4', 'PTPN23', 'PTPRK', 'PTPRM', 'PVR', 'RAB10', 'RAMP2', 'RDX', 'RND1', 'S100A11', 'SCRIB', 'SDCBP', 'SH3BP1', 'SHROOM1', 'SHROOM2', 'SHROOM3', 'SHROOM4', 'SMAD7', 'SNAP23', 'SORBS1', 'SPTBN4', 'SRC', 'SSX2IP', 'STXBP6', 'SYNM', 'TBCD', 'TJP1', 'TJP2', 'TLN1', 'TMEM204', 'TMEM47', 'TMOD3', 'TNK2', 'TNKS1BP1', 'TRIM29', 'TRPV4', 'TSPAN33', 'VCL', 'VEGFA', 'VEZT', 'WNK3', 'WTIP', 'ZNF703', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.adherens.heatmap.pdf', width = 9, height = 6)

genes = c('ABL1', 'ACHE', 'ACVR1', 'ACVR2A', 'ACVR2B', 'ADAR', 'AKT1', 'ALPL', 'ALYREF', 'AMELX', 'ASF1A', 'ATF4', 'ATP5F1B', 'ATP6AP1', 'ATRAID', 'AXIN2', 'BAMBI', 'BCAP29', 'BCL2', 'BGLAP', 'BMP2', 'BMP3', 'BMP4', 'BMP6', 'BMP7', 'BMPR1A', 'BMPR1B', 'BMPR2', 'CAT', 'CBFB', 'CCDC47', 'CCL3', 'CCN1', 'CCN4', 'CDK6', 'CEBPA', 'CEBPB', 'CHRD', 'CITED1', 'CLEC5A', 'CLIC1', 'CLTC', 'COL1A1', 'COL6A1', 'CREB3L1', 'CRIM1', 'CTHRC1', 'CTNNBIP1', 'CYP24A1', 'DDR2', 'DDX21', 'DDX5', 'DHH', 'DHX9', 'DLX5', 'DNAI3', 'DNAJC13', 'EIF2AK2', 'EPHA2', 'FASN', 'FBL', 'FBN2', 'FBXO5', 'FERMT2', 'FFAR4', 'FGF23', 'FGFR2', 'FHL2', 'FIGNL1', 'FZD1', 'GATA1', 'GDF10', 'GDF2', 'GDPD2', 'GLI1', 'GLI2', 'GLI3', 'GNAS', 'GREM1', 'GTPBP4', 'H3-3A', 'H3-3B', 'HAND2', 'HDAC4', 'HDAC7', 'HEMGN', 'HGF', 'HNRNPC', 'HNRNPU', 'HOXA2', 'HPSE', 'HSD17B4', 'HSPE1', 'IARS1', 'IBSP', 'ID3', 'IFITM1', 'IFT80', 'IGF1', 'IGF2', 'IGFBP3', 'IGFBP5', 'IHH', 'IL6', 'IL6R', 'IL6ST', 'ITGA11', 'ITGAV', 'JAG1', 'JUNB', 'JUND', 'LEF1', 'LGR4', 'LIMD1', 'LOX', 'LRP3', 'LRP5', 'LRP5L', 'LTF', 'MEF2C', 'MEF2D', 'MEN1', 'MIR138-1', 'MIR200C', 'MIR208A', 'MIR20A', 'MIR20B', 'MIR21', 'MIR210', 'MIR27A', 'MIR29B1', 'MIR548D1', 'MIR665', 'MIR675', 'MIR9-1', 'MN1', 'MRC2', 'MSX2', 'MUS81', 'MYBBP1A', 'MYOC', 'NBR1', 'NELL1', 'NF1', 'NOCT', 'NOG', 'NOTCH1', 'NPNT', 'NPPC', 'NPR3', 'OSR2', 'OSTN', 'PDLIM7', 'PENK', 'PHB', 'PLXNB1', 'PPARG', 'PRKACA', 'PRKD1', 'PSMC2', 'PTCH1', 'PTH1R', 'PTHLH', 'PTK2', 'RANBP3L', 'RASSF2', 'RBMX', 'RDH14', 'REST', 'RHOA', 'RIOX1', 'RORB', 'RPS15', 'RRAS2', 'RRBP1', 'RSL1D1', 'RSPO2', 'RUNX2', 'SATB2', 'SEMA4D', 'SEMA7A', 'SFRP1', 'SFRP2', 'SHH', 'SHOX2', 'SIRT7', 'SKI', 'SMAD1', 'SMAD3', 'SMAD5', 'SMAD6', 'SMO', 'SMOC1', 'SNAI1', 'SNAI2', 'SND1', 'SNRNP200', 'SOX11', 'SOX2', 'SOX8', 'SOX9', 'SP7', 'SPP1', 'SUCO', 'SUFU', 'SYNCRIP', 'TCIRG1', 'TMEM119', 'TMEM64', 'TNC', 'TNF', 'TNN', 'TOB1', 'TP53INP2', 'TP63', 'TPM4', 'TRPM4', 'TWIST1', 'TWIST2', 'TWSG1', 'UCMA', 'UFL1', 'VCAN', 'VEGFC', 'WNT10B', 'WNT11', 'WNT3', 'WNT3A', 'WNT4', 'WNT7B', 'WWOX', 'WWTR1', 'YAP1', 'ZHX3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.osteo1.heatmap.pdf', width = 9, height = 6)

genes = c('A2M', 'ABL1', 'ACAN', 'ACER2', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'ADAM10', 'ADAM15', 'ADAM8', 'ADAM9', 'ADAMTS12', 'ADAMTS13', 'ADAMTS4', 'ADAMTS5', 'ADAMTS9', 'ADTRP', 'AEBP1', 'AGT', 'AJAP1', 'AJUBA', 'AMBN', 'AMELX', 'AMELY', 'ANGPTL3', 'ANGPTL7', 'ANTXR1', 'APOD', 'ARHGAP6', 'ARHGEF7', 'ASPN', 'ATP7A', 'BCAM', 'BCAS3', 'BCL2', 'BCL2L11', 'BCL6', 'BCR', 'BGN', 'BMP1', 'BMP2', 'BSG', 'BST1', 'CAMSAP3', 'CAPN1', 'CAPN2', 'CAPNS1', 'CAPNS2', 'CARMIL2', 'CASK', 'CCL21', 'CCL25', 'CCL28', 'CCN2', 'CCR7', 'CD34', 'CD36', 'CD3E', 'CD44', 'CD63', 'CD96', 'CDH13', 'CDK5', 'CDK6', 'CDKN2A', 'CEACAM6', 'CFLAR', 'CHADL', 'CIB1', 'CLASP1', 'CLASP2', 'CMA1', 'COL10A1', 'COL11A1', 'COL11A2', 'COL12A1', 'COL13A1', 'COL14A1', 'COL15A1', 'COL16A1', 'COL17A1', 'COL18A1', 'COL19A1', 'COL1A1', 'COL1A2', 'COL21A1', 'COL23A1', 'COL24A1', 'COL25A1', 'COL27A1', 'COL28A1', 'COL2A1', 'COL3A1', 'COL4A1', 'COL4A2', 'COL4A3', 'COL4A4', 'COL4A5', 'COL4A6', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2', 'COL6A3', 'COL6A5', 'COL6A6', 'COL7A1', 'COL8A1', 'COL8A2', 'COL9A1', 'COL9A2', 'COL9A3', 'COLGALT1', 'CORO1C', 'CORO2B', 'CPB2', 'CSF1', 'CST3', 'CTRB1', 'CTRB2', 'CTSG', 'CTSK', 'CTSL', 'CTSS', 'CTSV', 'CTTN', 'CX3CL1', 'DAG1', 'DAPK3', 'DCN', 'DDR1', 'DDR2', 'DEFB118', 'DISC1', 'DLC1', 'DMTN', 'DPP4', 'DUSP22', 'DUSP3', 'ECM2', 'EDA', 'EFEMP2', 'EFNA5', 'ELANE', 'ELN', 'EMILIN1', 'EMILIN2', 'EMILIN3', 'EMP2', 'ENAM', 'EPB41L5', 'EPDR1', 'EPHA1', 'EPHA3', 'ETS1', 'EXOC8', 'FAM107A', 'FAP', 'FBLN2', 'FBLN5', 'FBN1', 'FBN2', 'FERMT1', 'FERMT2', 'FERMT3', 'FGA', 'FGB', 'FGFR4', 'FGG', 'FKBP10', 'FLOT1', 'FMN1', 'FMOD', 'FN1', 'FREM1', 'FSCN1', 'FURIN', 'FUT1', 'GAS6', 'GFUS', 'GPM6B', 'GREM1', 'GSK3B', 'HAPLN1', 'HAPLN2', 'HAS1', 'HAS2', 'HAS3', 'HOXA7', 'HOXD3', 'HPN', 'HPSE', 'HRG', 'HSPG2', 'HTRA1', 'IHH', 'IL6', 'ILK', 'IQGAP1', 'ITGA1', 'ITGA10', 'ITGA11', 'ITGA2', 'ITGA2B', 'ITGA3', 'ITGA4', 'ITGA7', 'ITGA8', 'ITGAL', 'ITGAV', 'ITGB1', 'ITGB1BP1', 'ITGB2', 'ITGB3', 'ITGB4', 'ITGB5', 'ITGB6', 'ITGB7', 'ITGBL1', 'JAG1', 'JAM3', 'JUP', 'KDR', 'KIF9', 'KLK2', 'KLK4', 'KLK5', 'KLK7', 'KLKB1', 'L1CAM', 'LAMB1', 'LAMB2', 'LAMB3', 'LAMB4', 'LAMC1', 'LCP1', 'LDB1', 'LIMCH1', 'LIMS1', 'LOX', 'LRP1', 'LTBP3', 'LUM', 'LYPD3', 'LYPD5', 'LYVE1', 'MACF1', 'MADCAM1', 'MAP4K4', 'MELTF', 'MFAP4', 'MINK1', 'MIR192', 'MIR29B1', 'MIR29C', 'MIR92A1', 'MIR939', 'MIR98', 'MKLN1', 'MMP1', 'MMP10', 'MMP11', 'MMP12', 'MMP13', 'MMP14', 'MMP15', 'MMP16', 'MMP19', 'MMP2', 'MMP20', 'MMP3', 'MMP7', 'MMP8', 'MMP9', 'MSLN', 'MSLNL', 'MUC4', 'MYF5', 'MYH11', 'MYOC', 'NEXMIF', 'NF1', 'NF2', 'NID1', 'NID2', 'NOTCH1', 'NOXO1', 'NPNT', 'NRP1', 'NTN4', 'NTNG1', 'NTNG2', 'OGN', 'ONECUT1', 'ONECUT2', 'OTOA', 'PARVG', 'PDPK1', 'PDPN', 'PEAK1', 'PHLDB1', 'PHLDB2', 'PIK3CB', 'PIK3R1', 'PIP5K1A', 'PKD1', 'PKHD1', 'PLAU', 'PLEKHA2', 'PLET1', 'PLG', 'PLOD3', 'PODN', 'POSTN', 'PPFIA1', 'PPFIA2', 'PPM1F', 'PRELP', 'PRG2', 'PRG3', 'PRG4', 'PRKCZ', 'PRSS1', 'PRSS2', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'PXDN', 'PXN', 'QSOX1', 'RAC1', 'RAMP2', 'RASA1', 'RB1', 'RCC2', 'RGCC', 'RHOA', 'RHOD', 'RIC1', 'RIN2', 'ROCK1', 'ROCK2', 'RRAS', 'RUNX1', 'S100A10', 'SCUBE1', 'SCUBE3', 'SDC4', 'SEMA3E', 'SERPINE1', 'SFRP1', 'SGCE', 'SH3PXD2B', 'SIGLEC1', 'SKAP1', 'SLC2A10', 'SLC9A1', 'SLK', 'SMAD3', 'SMPD3', 'SNED1', 'SORBS1', 'SOX9', 'SRC', 'SRF', 'STATH', 'STRC', 'STRCP1', 'TAOK2', 'TCF15', 'TECTA', 'TEK', 'TESK2', 'TGFB1', 'THBS1', 'THBS3', 'THSD1', 'THSD4', 'THY1', 'TIAM1', 'TIE1', 'TIMM10B', 'TIMP1', 'TIMP2', 'TLL1', 'TLL2', 'TMEM8B', 'TMPRSS6', 'TNFRSF1A', 'TNFRSF1B', 'TNN', 'TNXB', 'TPSAB1', 'TRIP6', 'TRPM7', 'TSC1', 'TUFT1', 'UTRN', 'VCAM1', 'VCAN', 'VCL', 'VEGFA', 'VTN', 'VWA2', 'WASHC1', 'WDPCP', 'WHAMM', 'WNT4', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.matrix.heatmap.pdf', width = 9, height = 6)

genes = c('ACVRL1', 'ADAMTS12', 'ADAMTS7', 'AMELX', 'ANXA6', 'AXIN2', 'BMP2', 'BMP4', 'BMP6', 'BMPR1A', 'BMPR1B', 'BMPR2', 'BPNT2', 'CCN2', 'CCN3', 'CCN4', 'CHADL', 'CHST11', 'CHSY1', 'COL11A1', 'COL27A1', 'COL2A1', 'COMP', 'CREB3L2', 'CYTL1', 'DDR2', 'ECM1', 'EFEMP1', 'EIF2AK3', 'EXT1', 'FGF18', 'FGFR3', 'GDF5', 'GDF6', 'GLG1', 'GLI2', 'GLI3', 'GPLD1', 'GREM1', 'HMGA2', 'HOXA11', 'IFT80', 'IHH', 'LNPK', 'LOXL2', 'LTBP3', 'LTF', 'MAF', 'MAPK14', 'MATN1', 'MBOAT2', 'MDK', 'MEF2C', 'MEF2D', 'MEX3C', 'MIR21', 'MMP14', 'MMP16', 'MUSTN1', 'NFIB', 'NKX3-2', 'NPPC', 'OSR1', 'OSR2', 'PKDCC', 'PTH', 'PTH1R', 'PTHLH', 'PTPN11', 'RARB', 'RARG', 'RFLNA', 'RFLNB', 'RUNX1', 'RUNX2', 'RUNX3', 'SCIN', 'SCX', 'SERPINH1', 'SFRP2', 'SHOX2', 'SIRT6', 'SIX2', 'SLC39A14', 'SMAD3', 'SMAD7', 'SMPD3', 'SNAI2', 'SNX19', 'SOX5', 'SOX6', 'SOX9', 'STC1', 'SULF1', 'SULF2', 'TGFB1', 'TGFBI', 'TGFBR1', 'TGFBR2', 'TRIP11', 'TRPS1', 'TSKU', 'WNT10B', 'WNT2B', 'WNT5B', 'WNT7A', 'WNT9A', 'ZBTB16', 'ZNF219')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_cell_line.chondro.heatmap.pdf', width = 9, height = 6)



###cunn_flexcell_rnaseq_0321.CTRL
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_CTRL.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_CTRL.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.CTRL.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.CTRL.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.CTRL.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.CTRL.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.CTRL.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19911,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.CTRL.Strained_vs_Unstrained.csv")




###cunn_flexcell_rnaseq_0321.FLNA
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNA.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNA.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.FLNA.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.FLNA.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNA.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.FLNA.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNA.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19885,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.FLNA.Strained_vs_Unstrained.csv")


###cunn_flexcell_rnaseq_0321.FLNB
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNB.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNB.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.FLNB.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.FLNB.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNB.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.FLNB.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNB.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19444,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.FLNB.Strained_vs_Unstrained.csv")


###cunn_flexcell_rnaseq_0321.FLNC
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNC.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_FLNC.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.FLNC.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.FLNC.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNC.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.FLNC.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.FLNC.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:19586,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.FLNC.Strained_vs_Unstrained.csv")


###cunn_flexcell_rnaseq_0321.PIEZO1
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_PIEZO1.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_PIEZO1.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 3, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321.PIEZO1.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
#write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321.PIEZO1.deseq.vst_counts.csv")
##calculate sample distances from rlog 
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##and plot as heat map
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Sex, rld$sample_name ,  sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.PIEZO1.sample_heatmap.pdf', width = 7, height = 5)
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321.PIEZO1.Cohort_pca.pdf', width=6, height = 6)

##gene clustering
##take 25 most variable gene
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),25)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
newdf <- as.data.frame(colData(rld)[c("Condition", "sample_only", "Sex", "Cohort")])
pheatmap(mat, annotation_col=newdf, fontsize_row=8, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321.PIEZO1.25_var_gene_clustering.pdf', width = 10, height = 10)

##differential expression
##do the test
dds <- DESeq(dds)
##get results and summary
##this just gets the last test
(res <- results(dds))
##to get a specific test:
res2 <- results(dds, contrast=c("Condition", "Strained", "Unstrained"))
##get summary
summary(res2) #lots of significant genes
##save differentially expression results
##sort results by adjusted p-value
resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)
##save results as dataframe and take top 20k results, then write csv file
resOrdered2DF <- as.data.frame(resOrdered2)[1:20253,]
write.csv(resOrdered2DF, file="cunn_flexcell_rnaseq_0321.PIEZO1.Strained_vs_Unstrained.csv")



###cunn_flexcell_rnaseq_0321_ctl_flna
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_flna.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Gene + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
##and write to csv file
write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321_ctl_flna.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flna.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flna.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flna.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flna.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flna.Cohort_pca.pdf', width=6, height = 6)

##clustering/heatmap from lists
genes = c('ABL1', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'AJUBA', 'APOD', 'ARHGAP6', 'ARHGEF7', 'BCAS3', 'BCL2', 'BCR', 'CAMSAP3', 'CLASP1', 'CLASP2', 'COL16A1', 'CORO1C', 'CORO2B', 'CTTN', 'DAPK3', 'DLC1', 'DMTN', 'DUSP22', 'DUSP3', 'EFNA5', 'EPB41L5', 'EPHA3', 'FAM107A', 'FERMT2', 'FMN1', 'GPM6B', 'GREM1', 'HRG', 'IQGAP1', 'ITGA2', 'ITGB1BP1', 'KDR', 'LDB1', 'LIMCH1', 'LIMS1', 'LRP1', 'MACF1', 'MAP4K4', 'MMP14', 'MYOC', 'NRP1', 'PDPK1', 'PEAK1', 'PHLDB2', 'PIP5K1A', 'PPM1F', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'RAC1', 'RCC2', 'RHOA', 'RHOD', 'ROCK1', 'ROCK2', 'S100A10', 'SDC4', 'SFRP1', 'SLC9A1', 'SLK', 'SMAD3', 'SORBS1', 'SRC', 'TAOK2', 'TEK', 'TESK2', 'THBS1', 'THSD1', 'THY1', 'TRIP6', 'TSC1', 'VCL', 'VEGFA', 'WDPCP', 'WHAMM', 'WNT4')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.focal_adhesion.heatmap.pdf', width = 9, height = 6)

genes = c('ABCB1', 'ACE', 'AGO3', 'ARIH2', 'ATXN1L', 'CCNE1', 'CD34', 'CITED1', 'CTC1', 'EIF2AK2', 'EPCAM', 'ETV6', 'FBLN1', 'FERMT1', 'FERMT2', 'FGF2', 'FZD1', 'GBA', 'GJA1', 'HMGA2', 'HMGB2', 'HNRNPU', 'KAT7', 'KDF1', 'KDM1A', 'KITLG', 'KRT18', 'LTBP3', 'MECOM', 'MIR16-1', 'MIR221', 'MIR222', 'MIR29B1', 'N4BP2L2', 'NANOG', 'NES', 'NF2', 'NR2E1', 'OVOL1', 'OVOL2', 'PDCD2', 'PDGFRA', 'PIM1', 'PTPRC', 'REST', 'RNF43', 'RUNX1', 'SFRP2', 'SIX2', 'SOX11', 'SOX17', 'SOX18', 'SOX5', 'SOX6', 'SOX9', 'TBX3', 'TERT', 'THPO', 'TRIM71', 'VEGFC', 'WNT1', 'WNT10B', 'WNT2B', 'WNT3', 'WNT5A', 'WNT7B', 'YAP1', 'YJEFN3', 'YTHDF2', 'ZFP36L1', 'ZNRF3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.msc_stem.heatmap.pdf', width = 9, height = 6)

genes = c('ABI2', 'ACTB', 'ADAM10', 'ADAM15', 'ADD1', 'AFDN', 'AHI1', 'AJAP1', 'AJM1', 'AJUBA', 'ALOX15B', 'ANG', 'ANXA1', 'ANXA2', 'APC', 'ARHGAP24', 'ARVCF', 'BAIAP2', 'BAIAP2L1', 'BMP6', 'BMPR2', 'CADM1', 'CADM2', 'CADM3', 'CAMSAP3', 'CCDC85A', 'CCDC85B', 'CCDC85C', 'CD99L2', 'CDC42', 'CDC42EP1', 'CDC42EP4', 'CDCA3', 'CDH1', 'CDH10', 'CDH11', 'CDH12', 'CDH13', 'CDH15', 'CDH17', 'CDH18', 'CDH19', 'CDH2', 'CDH20', 'CDH22', 'CDH24', 'CDH3', 'CDH4', 'CDH5', 'CDH6', 'CDH7', 'CDH8', 'CDH9', 'CDHR3', 'CEACAM1', 'CNN3', 'CRB1', 'CSK', 'CTNNA1', 'CTNNA2', 'CTNNA3', 'CTNNB1', 'CTNND1', 'CTNND2', 'CXADR', 'CYTH1', 'CYTH2', 'CYTH3', 'DAG1', 'DCHS1', 'DDX6', 'DLG1', 'DLG5', 'DLL1', 'DSC2', 'DSP', 'EFNB2', 'EIF4G2', 'EPB41L5', 'EPHA4', 'ESAM', 'FAT2', 'FERMT2', 'FLOT1', 'FLOT2', 'FMN1', 'FRMD4A', 'FRMD4B', 'FRMD5', 'FRS2', 'HIPK1', 'HMCN1', 'IGSF21', 'INAVA', 'ITGA6', 'JAG1', 'JAM3', 'JCAD', 'JUP', 'KIFC3', 'KLHL24', 'KRT18', 'LDB3', 'LIMD1', 'LIN7A', 'LIN7B', 'LIN7C', 'LYN', 'MAGI1', 'MPP4', 'MPP5', 'MPP7', 'MYH9', 'MYO1E', 'NDRG1', 'NECTIN1', 'NECTIN2', 'NECTIN3', 'NECTIN4', 'NEXN', 'NF2', 'NIBAN2', 'NOTCH1', 'NPHP1', 'NUMB', 'NUMBL', 'OXTR', 'PAK2', 'PAK4', 'PARD3', 'PARD3B', 'PARK7', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDLIM5', 'PDLIM7', 'PDZD11', 'PGM5', 'PIP5K1C', 'PKP1', 'PKP2', 'PKP3', 'PKP4', 'PLEKHA7', 'PLPP3', 'POF1B', 'PPP1CA', 'PPP1R9B', 'PRICKLE4', 'PTPN23', 'PTPRK', 'PTPRM', 'PVR', 'RAB10', 'RAMP2', 'RDX', 'RND1', 'S100A11', 'SCRIB', 'SDCBP', 'SH3BP1', 'SHROOM1', 'SHROOM2', 'SHROOM3', 'SHROOM4', 'SMAD7', 'SNAP23', 'SORBS1', 'SPTBN4', 'SRC', 'SSX2IP', 'STXBP6', 'SYNM', 'TBCD', 'TJP1', 'TJP2', 'TLN1', 'TMEM204', 'TMEM47', 'TMOD3', 'TNK2', 'TNKS1BP1', 'TRIM29', 'TRPV4', 'TSPAN33', 'VCL', 'VEGFA', 'VEZT', 'WNK3', 'WTIP', 'ZNF703', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.adherens.heatmap.pdf', width = 9, height = 6)

genes = c('ABL1', 'ACHE', 'ACVR1', 'ACVR2A', 'ACVR2B', 'ADAR', 'AKT1', 'ALPL', 'ALYREF', 'AMELX', 'ASF1A', 'ATF4', 'ATP5F1B', 'ATP6AP1', 'ATRAID', 'AXIN2', 'BAMBI', 'BCAP29', 'BCL2', 'BGLAP', 'BMP2', 'BMP3', 'BMP4', 'BMP6', 'BMP7', 'BMPR1A', 'BMPR1B', 'BMPR2', 'CAT', 'CBFB', 'CCDC47', 'CCL3', 'CCN1', 'CCN4', 'CDK6', 'CEBPA', 'CEBPB', 'CHRD', 'CITED1', 'CLEC5A', 'CLIC1', 'CLTC', 'COL1A1', 'COL6A1', 'CREB3L1', 'CRIM1', 'CTHRC1', 'CTNNBIP1', 'CYP24A1', 'DDR2', 'DDX21', 'DDX5', 'DHH', 'DHX9', 'DLX5', 'DNAI3', 'DNAJC13', 'EIF2AK2', 'EPHA2', 'FASN', 'FBL', 'FBN2', 'FBXO5', 'FERMT2', 'FFAR4', 'FGF23', 'FGFR2', 'FHL2', 'FIGNL1', 'FZD1', 'GATA1', 'GDF10', 'GDF2', 'GDPD2', 'GLI1', 'GLI2', 'GLI3', 'GNAS', 'GREM1', 'GTPBP4', 'H3-3A', 'H3-3B', 'HAND2', 'HDAC4', 'HDAC7', 'HEMGN', 'HGF', 'HNRNPC', 'HNRNPU', 'HOXA2', 'HPSE', 'HSD17B4', 'HSPE1', 'IARS1', 'IBSP', 'ID3', 'IFITM1', 'IFT80', 'IGF1', 'IGF2', 'IGFBP3', 'IGFBP5', 'IHH', 'IL6', 'IL6R', 'IL6ST', 'ITGA11', 'ITGAV', 'JAG1', 'JUNB', 'JUND', 'LEF1', 'LGR4', 'LIMD1', 'LOX', 'LRP3', 'LRP5', 'LRP5L', 'LTF', 'MEF2C', 'MEF2D', 'MEN1', 'MIR138-1', 'MIR200C', 'MIR208A', 'MIR20A', 'MIR20B', 'MIR21', 'MIR210', 'MIR27A', 'MIR29B1', 'MIR548D1', 'MIR665', 'MIR675', 'MIR9-1', 'MN1', 'MRC2', 'MSX2', 'MUS81', 'MYBBP1A', 'MYOC', 'NBR1', 'NELL1', 'NF1', 'NOCT', 'NOG', 'NOTCH1', 'NPNT', 'NPPC', 'NPR3', 'OSR2', 'OSTN', 'PDLIM7', 'PENK', 'PHB', 'PLXNB1', 'PPARG', 'PRKACA', 'PRKD1', 'PSMC2', 'PTCH1', 'PTH1R', 'PTHLH', 'PTK2', 'RANBP3L', 'RASSF2', 'RBMX', 'RDH14', 'REST', 'RHOA', 'RIOX1', 'RORB', 'RPS15', 'RRAS2', 'RRBP1', 'RSL1D1', 'RSPO2', 'RUNX2', 'SATB2', 'SEMA4D', 'SEMA7A', 'SFRP1', 'SFRP2', 'SHH', 'SHOX2', 'SIRT7', 'SKI', 'SMAD1', 'SMAD3', 'SMAD5', 'SMAD6', 'SMO', 'SMOC1', 'SNAI1', 'SNAI2', 'SND1', 'SNRNP200', 'SOX11', 'SOX2', 'SOX8', 'SOX9', 'SP7', 'SPP1', 'SUCO', 'SUFU', 'SYNCRIP', 'TCIRG1', 'TMEM119', 'TMEM64', 'TNC', 'TNF', 'TNN', 'TOB1', 'TP53INP2', 'TP63', 'TPM4', 'TRPM4', 'TWIST1', 'TWIST2', 'TWSG1', 'UCMA', 'UFL1', 'VCAN', 'VEGFC', 'WNT10B', 'WNT11', 'WNT3', 'WNT3A', 'WNT4', 'WNT7B', 'WWOX', 'WWTR1', 'YAP1', 'ZHX3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.osteo1.heatmap.pdf', width = 9, height = 6)

genes = c('A2M', 'ABL1', 'ACAN', 'ACER2', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'ADAM10', 'ADAM15', 'ADAM8', 'ADAM9', 'ADAMTS12', 'ADAMTS13', 'ADAMTS4', 'ADAMTS5', 'ADAMTS9', 'ADTRP', 'AEBP1', 'AGT', 'AJAP1', 'AJUBA', 'AMBN', 'AMELX', 'AMELY', 'ANGPTL3', 'ANGPTL7', 'ANTXR1', 'APOD', 'ARHGAP6', 'ARHGEF7', 'ASPN', 'ATP7A', 'BCAM', 'BCAS3', 'BCL2', 'BCL2L11', 'BCL6', 'BCR', 'BGN', 'BMP1', 'BMP2', 'BSG', 'BST1', 'CAMSAP3', 'CAPN1', 'CAPN2', 'CAPNS1', 'CAPNS2', 'CARMIL2', 'CASK', 'CCL21', 'CCL25', 'CCL28', 'CCN2', 'CCR7', 'CD34', 'CD36', 'CD3E', 'CD44', 'CD63', 'CD96', 'CDH13', 'CDK5', 'CDK6', 'CDKN2A', 'CEACAM6', 'CFLAR', 'CHADL', 'CIB1', 'CLASP1', 'CLASP2', 'CMA1', 'COL10A1', 'COL11A1', 'COL11A2', 'COL12A1', 'COL13A1', 'COL14A1', 'COL15A1', 'COL16A1', 'COL17A1', 'COL18A1', 'COL19A1', 'COL1A1', 'COL1A2', 'COL21A1', 'COL23A1', 'COL24A1', 'COL25A1', 'COL27A1', 'COL28A1', 'COL2A1', 'COL3A1', 'COL4A1', 'COL4A2', 'COL4A3', 'COL4A4', 'COL4A5', 'COL4A6', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2', 'COL6A3', 'COL6A5', 'COL6A6', 'COL7A1', 'COL8A1', 'COL8A2', 'COL9A1', 'COL9A2', 'COL9A3', 'COLGALT1', 'CORO1C', 'CORO2B', 'CPB2', 'CSF1', 'CST3', 'CTRB1', 'CTRB2', 'CTSG', 'CTSK', 'CTSL', 'CTSS', 'CTSV', 'CTTN', 'CX3CL1', 'DAG1', 'DAPK3', 'DCN', 'DDR1', 'DDR2', 'DEFB118', 'DISC1', 'DLC1', 'DMTN', 'DPP4', 'DUSP22', 'DUSP3', 'ECM2', 'EDA', 'EFEMP2', 'EFNA5', 'ELANE', 'ELN', 'EMILIN1', 'EMILIN2', 'EMILIN3', 'EMP2', 'ENAM', 'EPB41L5', 'EPDR1', 'EPHA1', 'EPHA3', 'ETS1', 'EXOC8', 'FAM107A', 'FAP', 'FBLN2', 'FBLN5', 'FBN1', 'FBN2', 'FERMT1', 'FERMT2', 'FERMT3', 'FGA', 'FGB', 'FGFR4', 'FGG', 'FKBP10', 'FLOT1', 'FMN1', 'FMOD', 'FN1', 'FREM1', 'FSCN1', 'FURIN', 'FUT1', 'GAS6', 'GFUS', 'GPM6B', 'GREM1', 'GSK3B', 'HAPLN1', 'HAPLN2', 'HAS1', 'HAS2', 'HAS3', 'HOXA7', 'HOXD3', 'HPN', 'HPSE', 'HRG', 'HSPG2', 'HTRA1', 'IHH', 'IL6', 'ILK', 'IQGAP1', 'ITGA1', 'ITGA10', 'ITGA11', 'ITGA2', 'ITGA2B', 'ITGA3', 'ITGA4', 'ITGA7', 'ITGA8', 'ITGAL', 'ITGAV', 'ITGB1', 'ITGB1BP1', 'ITGB2', 'ITGB3', 'ITGB4', 'ITGB5', 'ITGB6', 'ITGB7', 'ITGBL1', 'JAG1', 'JAM3', 'JUP', 'KDR', 'KIF9', 'KLK2', 'KLK4', 'KLK5', 'KLK7', 'KLKB1', 'L1CAM', 'LAMB1', 'LAMB2', 'LAMB3', 'LAMB4', 'LAMC1', 'LCP1', 'LDB1', 'LIMCH1', 'LIMS1', 'LOX', 'LRP1', 'LTBP3', 'LUM', 'LYPD3', 'LYPD5', 'LYVE1', 'MACF1', 'MADCAM1', 'MAP4K4', 'MELTF', 'MFAP4', 'MINK1', 'MIR192', 'MIR29B1', 'MIR29C', 'MIR92A1', 'MIR939', 'MIR98', 'MKLN1', 'MMP1', 'MMP10', 'MMP11', 'MMP12', 'MMP13', 'MMP14', 'MMP15', 'MMP16', 'MMP19', 'MMP2', 'MMP20', 'MMP3', 'MMP7', 'MMP8', 'MMP9', 'MSLN', 'MSLNL', 'MUC4', 'MYF5', 'MYH11', 'MYOC', 'NEXMIF', 'NF1', 'NF2', 'NID1', 'NID2', 'NOTCH1', 'NOXO1', 'NPNT', 'NRP1', 'NTN4', 'NTNG1', 'NTNG2', 'OGN', 'ONECUT1', 'ONECUT2', 'OTOA', 'PARVG', 'PDPK1', 'PDPN', 'PEAK1', 'PHLDB1', 'PHLDB2', 'PIK3CB', 'PIK3R1', 'PIP5K1A', 'PKD1', 'PKHD1', 'PLAU', 'PLEKHA2', 'PLET1', 'PLG', 'PLOD3', 'PODN', 'POSTN', 'PPFIA1', 'PPFIA2', 'PPM1F', 'PRELP', 'PRG2', 'PRG3', 'PRG4', 'PRKCZ', 'PRSS1', 'PRSS2', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'PXDN', 'PXN', 'QSOX1', 'RAC1', 'RAMP2', 'RASA1', 'RB1', 'RCC2', 'RGCC', 'RHOA', 'RHOD', 'RIC1', 'RIN2', 'ROCK1', 'ROCK2', 'RRAS', 'RUNX1', 'S100A10', 'SCUBE1', 'SCUBE3', 'SDC4', 'SEMA3E', 'SERPINE1', 'SFRP1', 'SGCE', 'SH3PXD2B', 'SIGLEC1', 'SKAP1', 'SLC2A10', 'SLC9A1', 'SLK', 'SMAD3', 'SMPD3', 'SNED1', 'SORBS1', 'SOX9', 'SRC', 'SRF', 'STATH', 'STRC', 'STRCP1', 'TAOK2', 'TCF15', 'TECTA', 'TEK', 'TESK2', 'TGFB1', 'THBS1', 'THBS3', 'THSD1', 'THSD4', 'THY1', 'TIAM1', 'TIE1', 'TIMM10B', 'TIMP1', 'TIMP2', 'TLL1', 'TLL2', 'TMEM8B', 'TMPRSS6', 'TNFRSF1A', 'TNFRSF1B', 'TNN', 'TNXB', 'TPSAB1', 'TRIP6', 'TRPM7', 'TSC1', 'TUFT1', 'UTRN', 'VCAM1', 'VCAN', 'VCL', 'VEGFA', 'VTN', 'VWA2', 'WASHC1', 'WDPCP', 'WHAMM', 'WNT4', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.matrix.heatmap.pdf', width = 9, height = 6)

genes = c('ACVRL1', 'ADAMTS12', 'ADAMTS7', 'AMELX', 'ANXA6', 'AXIN2', 'BMP2', 'BMP4', 'BMP6', 'BMPR1A', 'BMPR1B', 'BMPR2', 'BPNT2', 'CCN2', 'CCN3', 'CCN4', 'CHADL', 'CHST11', 'CHSY1', 'COL11A1', 'COL27A1', 'COL2A1', 'COMP', 'CREB3L2', 'CYTL1', 'DDR2', 'ECM1', 'EFEMP1', 'EIF2AK3', 'EXT1', 'FGF18', 'FGFR3', 'GDF5', 'GDF6', 'GLG1', 'GLI2', 'GLI3', 'GPLD1', 'GREM1', 'HMGA2', 'HOXA11', 'IFT80', 'IHH', 'LNPK', 'LOXL2', 'LTBP3', 'LTF', 'MAF', 'MAPK14', 'MATN1', 'MBOAT2', 'MDK', 'MEF2C', 'MEF2D', 'MEX3C', 'MIR21', 'MMP14', 'MMP16', 'MUSTN1', 'NFIB', 'NKX3-2', 'NPPC', 'OSR1', 'OSR2', 'PKDCC', 'PTH', 'PTH1R', 'PTHLH', 'PTPN11', 'RARB', 'RARG', 'RFLNA', 'RFLNB', 'RUNX1', 'RUNX2', 'RUNX3', 'SCIN', 'SCX', 'SERPINH1', 'SFRP2', 'SHOX2', 'SIRT6', 'SIX2', 'SLC39A14', 'SMAD3', 'SMAD7', 'SMPD3', 'SNAI2', 'SNX19', 'SOX5', 'SOX6', 'SOX9', 'STC1', 'SULF1', 'SULF2', 'TGFB1', 'TGFBI', 'TGFBR1', 'TGFBR2', 'TRIP11', 'TRPS1', 'TSKU', 'WNT10B', 'WNT2B', 'WNT5B', 'WNT7A', 'WNT9A', 'ZBTB16', 'ZNF219')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flna.chondro.heatmap.pdf', width = 9, height = 6)

###cunn_flexcell_rnaseq_0321_ctl_flnb --- all samples minus bone and suture
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_flnb.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_flnb.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Gene + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321_ctl_flnb.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321_ctl_flnb.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnb.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnb.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnb.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnb.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnb.Cohort_pca.pdf', width=6, height = 6)

##clustering/heatmap from lists
genes = c('ABL1', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'AJUBA', 'APOD', 'ARHGAP6', 'ARHGEF7', 'BCAS3', 'BCL2', 'BCR', 'CAMSAP3', 'CLASP1', 'CLASP2', 'COL16A1', 'CORO1C', 'CORO2B', 'CTTN', 'DAPK3', 'DLC1', 'DMTN', 'DUSP22', 'DUSP3', 'EFNA5', 'EPB41L5', 'EPHA3', 'FAM107A', 'FERMT2', 'FMN1', 'GPM6B', 'GREM1', 'HRG', 'IQGAP1', 'ITGA2', 'ITGB1BP1', 'KDR', 'LDB1', 'LIMCH1', 'LIMS1', 'LRP1', 'MACF1', 'MAP4K4', 'MMP14', 'MYOC', 'NRP1', 'PDPK1', 'PEAK1', 'PHLDB2', 'PIP5K1A', 'PPM1F', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'RAC1', 'RCC2', 'RHOA', 'RHOD', 'ROCK1', 'ROCK2', 'S100A10', 'SDC4', 'SFRP1', 'SLC9A1', 'SLK', 'SMAD3', 'SORBS1', 'SRC', 'TAOK2', 'TEK', 'TESK2', 'THBS1', 'THSD1', 'THY1', 'TRIP6', 'TSC1', 'VCL', 'VEGFA', 'WDPCP', 'WHAMM', 'WNT4')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnb.focal_adhesion.heatmap.pdf', width = 9, height = 6)

genes = c('ABCB1', 'ACE', 'AGO3', 'ARIH2', 'ATXN1L', 'CCNE1', 'CD34', 'CITED1', 'CTC1', 'EIF2AK2', 'EPCAM', 'ETV6', 'FBLN1', 'FERMT1', 'FERMT2', 'FGF2', 'FZD1', 'GBA', 'GJA1', 'HMGA2', 'HMGB2', 'HNRNPU', 'KAT7', 'KDF1', 'KDM1A', 'KITLG', 'KRT18', 'LTBP3', 'MECOM', 'MIR16-1', 'MIR221', 'MIR222', 'MIR29B1', 'N4BP2L2', 'NANOG', 'NES', 'NF2', 'NR2E1', 'OVOL1', 'OVOL2', 'PDCD2', 'PDGFRA', 'PIM1', 'PTPRC', 'REST', 'RNF43', 'RUNX1', 'SFRP2', 'SIX2', 'SOX11', 'SOX17', 'SOX18', 'SOX5', 'SOX6', 'SOX9', 'TBX3', 'TERT', 'THPO', 'TRIM71', 'VEGFC', 'WNT1', 'WNT10B', 'WNT2B', 'WNT3', 'WNT5A', 'WNT7B', 'YAP1', 'YJEFN3', 'YTHDF2', 'ZFP36L1', 'ZNRF3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnb.msc_stem.heatmap.pdf', width = 9, height = 6)

genes = c('ABI2', 'ACTB', 'ADAM10', 'ADAM15', 'ADD1', 'AFDN', 'AHI1', 'AJAP1', 'AJM1', 'AJUBA', 'ALOX15B', 'ANG', 'ANXA1', 'ANXA2', 'APC', 'ARHGAP24', 'ARVCF', 'BAIAP2', 'BAIAP2L1', 'BMP6', 'BMPR2', 'CADM1', 'CADM2', 'CADM3', 'CAMSAP3', 'CCDC85A', 'CCDC85B', 'CCDC85C', 'CD99L2', 'CDC42', 'CDC42EP1', 'CDC42EP4', 'CDCA3', 'CDH1', 'CDH10', 'CDH11', 'CDH12', 'CDH13', 'CDH15', 'CDH17', 'CDH18', 'CDH19', 'CDH2', 'CDH20', 'CDH22', 'CDH24', 'CDH3', 'CDH4', 'CDH5', 'CDH6', 'CDH7', 'CDH8', 'CDH9', 'CDHR3', 'CEACAM1', 'CNN3', 'CRB1', 'CSK', 'CTNNA1', 'CTNNA2', 'CTNNA3', 'CTNNB1', 'CTNND1', 'CTNND2', 'CXADR', 'CYTH1', 'CYTH2', 'CYTH3', 'DAG1', 'DCHS1', 'DDX6', 'DLG1', 'DLG5', 'DLL1', 'DSC2', 'DSP', 'EFNB2', 'EIF4G2', 'EPB41L5', 'EPHA4', 'ESAM', 'FAT2', 'FERMT2', 'FLOT1', 'FLOT2', 'FMN1', 'FRMD4A', 'FRMD4B', 'FRMD5', 'FRS2', 'HIPK1', 'HMCN1', 'IGSF21', 'INAVA', 'ITGA6', 'JAG1', 'JAM3', 'JCAD', 'JUP', 'KIFC3', 'KLHL24', 'KRT18', 'LDB3', 'LIMD1', 'LIN7A', 'LIN7B', 'LIN7C', 'LYN', 'MAGI1', 'MPP4', 'MPP5', 'MPP7', 'MYH9', 'MYO1E', 'NDRG1', 'NECTIN1', 'NECTIN2', 'NECTIN3', 'NECTIN4', 'NEXN', 'NF2', 'NIBAN2', 'NOTCH1', 'NPHP1', 'NUMB', 'NUMBL', 'OXTR', 'PAK2', 'PAK4', 'PARD3', 'PARD3B', 'PARK7', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDLIM5', 'PDLIM7', 'PDZD11', 'PGM5', 'PIP5K1C', 'PKP1', 'PKP2', 'PKP3', 'PKP4', 'PLEKHA7', 'PLPP3', 'POF1B', 'PPP1CA', 'PPP1R9B', 'PRICKLE4', 'PTPN23', 'PTPRK', 'PTPRM', 'PVR', 'RAB10', 'RAMP2', 'RDX', 'RND1', 'S100A11', 'SCRIB', 'SDCBP', 'SH3BP1', 'SHROOM1', 'SHROOM2', 'SHROOM3', 'SHROOM4', 'SMAD7', 'SNAP23', 'SORBS1', 'SPTBN4', 'SRC', 'SSX2IP', 'STXBP6', 'SYNM', 'TBCD', 'TJP1', 'TJP2', 'TLN1', 'TMEM204', 'TMEM47', 'TMOD3', 'TNK2', 'TNKS1BP1', 'TRIM29', 'TRPV4', 'TSPAN33', 'VCL', 'VEGFA', 'VEZT', 'WNK3', 'WTIP', 'ZNF703', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnb.adherens.heatmap.pdf', width = 9, height = 6)

genes = c('ABL1', 'ACHE', 'ACVR1', 'ACVR2A', 'ACVR2B', 'ADAR', 'AKT1', 'ALPL', 'ALYREF', 'AMELX', 'ASF1A', 'ATF4', 'ATP5F1B', 'ATP6AP1', 'ATRAID', 'AXIN2', 'BAMBI', 'BCAP29', 'BCL2', 'BGLAP', 'BMP2', 'BMP3', 'BMP4', 'BMP6', 'BMP7', 'BMPR1A', 'BMPR1B', 'BMPR2', 'CAT', 'CBFB', 'CCDC47', 'CCL3', 'CCN1', 'CCN4', 'CDK6', 'CEBPA', 'CEBPB', 'CHRD', 'CITED1', 'CLEC5A', 'CLIC1', 'CLTC', 'COL1A1', 'COL6A1', 'CREB3L1', 'CRIM1', 'CTHRC1', 'CTNNBIP1', 'CYP24A1', 'DDR2', 'DDX21', 'DDX5', 'DHH', 'DHX9', 'DLX5', 'DNAI3', 'DNAJC13', 'EIF2AK2', 'EPHA2', 'FASN', 'FBL', 'FBN2', 'FBXO5', 'FERMT2', 'FFAR4', 'FGF23', 'FGFR2', 'FHL2', 'FIGNL1', 'FZD1', 'GATA1', 'GDF10', 'GDF2', 'GDPD2', 'GLI1', 'GLI2', 'GLI3', 'GNAS', 'GREM1', 'GTPBP4', 'H3-3A', 'H3-3B', 'HAND2', 'HDAC4', 'HDAC7', 'HEMGN', 'HGF', 'HNRNPC', 'HNRNPU', 'HOXA2', 'HPSE', 'HSD17B4', 'HSPE1', 'IARS1', 'IBSP', 'ID3', 'IFITM1', 'IFT80', 'IGF1', 'IGF2', 'IGFBP3', 'IGFBP5', 'IHH', 'IL6', 'IL6R', 'IL6ST', 'ITGA11', 'ITGAV', 'JAG1', 'JUNB', 'JUND', 'LEF1', 'LGR4', 'LIMD1', 'LOX', 'LRP3', 'LRP5', 'LRP5L', 'LTF', 'MEF2C', 'MEF2D', 'MEN1', 'MIR138-1', 'MIR200C', 'MIR208A', 'MIR20A', 'MIR20B', 'MIR21', 'MIR210', 'MIR27A', 'MIR29B1', 'MIR548D1', 'MIR665', 'MIR675', 'MIR9-1', 'MN1', 'MRC2', 'MSX2', 'MUS81', 'MYBBP1A', 'MYOC', 'NBR1', 'NELL1', 'NF1', 'NOCT', 'NOG', 'NOTCH1', 'NPNT', 'NPPC', 'NPR3', 'OSR2', 'OSTN', 'PDLIM7', 'PENK', 'PHB', 'PLXNB1', 'PPARG', 'PRKACA', 'PRKD1', 'PSMC2', 'PTCH1', 'PTH1R', 'PTHLH', 'PTK2', 'RANBP3L', 'RASSF2', 'RBMX', 'RDH14', 'REST', 'RHOA', 'RIOX1', 'RORB', 'RPS15', 'RRAS2', 'RRBP1', 'RSL1D1', 'RSPO2', 'RUNX2', 'SATB2', 'SEMA4D', 'SEMA7A', 'SFRP1', 'SFRP2', 'SHH', 'SHOX2', 'SIRT7', 'SKI', 'SMAD1', 'SMAD3', 'SMAD5', 'SMAD6', 'SMO', 'SMOC1', 'SNAI1', 'SNAI2', 'SND1', 'SNRNP200', 'SOX11', 'SOX2', 'SOX8', 'SOX9', 'SP7', 'SPP1', 'SUCO', 'SUFU', 'SYNCRIP', 'TCIRG1', 'TMEM119', 'TMEM64', 'TNC', 'TNF', 'TNN', 'TOB1', 'TP53INP2', 'TP63', 'TPM4', 'TRPM4', 'TWIST1', 'TWIST2', 'TWSG1', 'UCMA', 'UFL1', 'VCAN', 'VEGFC', 'WNT10B', 'WNT11', 'WNT3', 'WNT3A', 'WNT4', 'WNT7B', 'WWOX', 'WWTR1', 'YAP1', 'ZHX3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnb.osteo1.heatmap.pdf', width = 9, height = 6)

genes = c('A2M', 'ABL1', 'ACAN', 'ACER2', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'ADAM10', 'ADAM15', 'ADAM8', 'ADAM9', 'ADAMTS12', 'ADAMTS13', 'ADAMTS4', 'ADAMTS5', 'ADAMTS9', 'ADTRP', 'AEBP1', 'AGT', 'AJAP1', 'AJUBA', 'AMBN', 'AMELX', 'AMELY', 'ANGPTL3', 'ANGPTL7', 'ANTXR1', 'APOD', 'ARHGAP6', 'ARHGEF7', 'ASPN', 'ATP7A', 'BCAM', 'BCAS3', 'BCL2', 'BCL2L11', 'BCL6', 'BCR', 'BGN', 'BMP1', 'BMP2', 'BSG', 'BST1', 'CAMSAP3', 'CAPN1', 'CAPN2', 'CAPNS1', 'CAPNS2', 'CARMIL2', 'CASK', 'CCL21', 'CCL25', 'CCL28', 'CCN2', 'CCR7', 'CD34', 'CD36', 'CD3E', 'CD44', 'CD63', 'CD96', 'CDH13', 'CDK5', 'CDK6', 'CDKN2A', 'CEACAM6', 'CFLAR', 'CHADL', 'CIB1', 'CLASP1', 'CLASP2', 'CMA1', 'COL10A1', 'COL11A1', 'COL11A2', 'COL12A1', 'COL13A1', 'COL14A1', 'COL15A1', 'COL16A1', 'COL17A1', 'COL18A1', 'COL19A1', 'COL1A1', 'COL1A2', 'COL21A1', 'COL23A1', 'COL24A1', 'COL25A1', 'COL27A1', 'COL28A1', 'COL2A1', 'COL3A1', 'COL4A1', 'COL4A2', 'COL4A3', 'COL4A4', 'COL4A5', 'COL4A6', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2', 'COL6A3', 'COL6A5', 'COL6A6', 'COL7A1', 'COL8A1', 'COL8A2', 'COL9A1', 'COL9A2', 'COL9A3', 'COLGALT1', 'CORO1C', 'CORO2B', 'CPB2', 'CSF1', 'CST3', 'CTRB1', 'CTRB2', 'CTSG', 'CTSK', 'CTSL', 'CTSS', 'CTSV', 'CTTN', 'CX3CL1', 'DAG1', 'DAPK3', 'DCN', 'DDR1', 'DDR2', 'DEFB118', 'DISC1', 'DLC1', 'DMTN', 'DPP4', 'DUSP22', 'DUSP3', 'ECM2', 'EDA', 'EFEMP2', 'EFNA5', 'ELANE', 'ELN', 'EMILIN1', 'EMILIN2', 'EMILIN3', 'EMP2', 'ENAM', 'EPB41L5', 'EPDR1', 'EPHA1', 'EPHA3', 'ETS1', 'EXOC8', 'FAM107A', 'FAP', 'FBLN2', 'FBLN5', 'FBN1', 'FBN2', 'FERMT1', 'FERMT2', 'FERMT3', 'FGA', 'FGB', 'FGFR4', 'FGG', 'FKBP10', 'FLOT1', 'FMN1', 'FMOD', 'FN1', 'FREM1', 'FSCN1', 'FURIN', 'FUT1', 'GAS6', 'GFUS', 'GPM6B', 'GREM1', 'GSK3B', 'HAPLN1', 'HAPLN2', 'HAS1', 'HAS2', 'HAS3', 'HOXA7', 'HOXD3', 'HPN', 'HPSE', 'HRG', 'HSPG2', 'HTRA1', 'IHH', 'IL6', 'ILK', 'IQGAP1', 'ITGA1', 'ITGA10', 'ITGA11', 'ITGA2', 'ITGA2B', 'ITGA3', 'ITGA4', 'ITGA7', 'ITGA8', 'ITGAL', 'ITGAV', 'ITGB1', 'ITGB1BP1', 'ITGB2', 'ITGB3', 'ITGB4', 'ITGB5', 'ITGB6', 'ITGB7', 'ITGBL1', 'JAG1', 'JAM3', 'JUP', 'KDR', 'KIF9', 'KLK2', 'KLK4', 'KLK5', 'KLK7', 'KLKB1', 'L1CAM', 'LAMB1', 'LAMB2', 'LAMB3', 'LAMB4', 'LAMC1', 'LCP1', 'LDB1', 'LIMCH1', 'LIMS1', 'LOX', 'LRP1', 'LTBP3', 'LUM', 'LYPD3', 'LYPD5', 'LYVE1', 'MACF1', 'MADCAM1', 'MAP4K4', 'MELTF', 'MFAP4', 'MINK1', 'MIR192', 'MIR29B1', 'MIR29C', 'MIR92A1', 'MIR939', 'MIR98', 'MKLN1', 'MMP1', 'MMP10', 'MMP11', 'MMP12', 'MMP13', 'MMP14', 'MMP15', 'MMP16', 'MMP19', 'MMP2', 'MMP20', 'MMP3', 'MMP7', 'MMP8', 'MMP9', 'MSLN', 'MSLNL', 'MUC4', 'MYF5', 'MYH11', 'MYOC', 'NEXMIF', 'NF1', 'NF2', 'NID1', 'NID2', 'NOTCH1', 'NOXO1', 'NPNT', 'NRP1', 'NTN4', 'NTNG1', 'NTNG2', 'OGN', 'ONECUT1', 'ONECUT2', 'OTOA', 'PARVG', 'PDPK1', 'PDPN', 'PEAK1', 'PHLDB1', 'PHLDB2', 'PIK3CB', 'PIK3R1', 'PIP5K1A', 'PKD1', 'PKHD1', 'PLAU', 'PLEKHA2', 'PLET1', 'PLG', 'PLOD3', 'PODN', 'POSTN', 'PPFIA1', 'PPFIA2', 'PPM1F', 'PRELP', 'PRG2', 'PRG3', 'PRG4', 'PRKCZ', 'PRSS1', 'PRSS2', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'PXDN', 'PXN', 'QSOX1', 'RAC1', 'RAMP2', 'RASA1', 'RB1', 'RCC2', 'RGCC', 'RHOA', 'RHOD', 'RIC1', 'RIN2', 'ROCK1', 'ROCK2', 'RRAS', 'RUNX1', 'S100A10', 'SCUBE1', 'SCUBE3', 'SDC4', 'SEMA3E', 'SERPINE1', 'SFRP1', 'SGCE', 'SH3PXD2B', 'SIGLEC1', 'SKAP1', 'SLC2A10', 'SLC9A1', 'SLK', 'SMAD3', 'SMPD3', 'SNED1', 'SORBS1', 'SOX9', 'SRC', 'SRF', 'STATH', 'STRC', 'STRCP1', 'TAOK2', 'TCF15', 'TECTA', 'TEK', 'TESK2', 'TGFB1', 'THBS1', 'THBS3', 'THSD1', 'THSD4', 'THY1', 'TIAM1', 'TIE1', 'TIMM10B', 'TIMP1', 'TIMP2', 'TLL1', 'TLL2', 'TMEM8B', 'TMPRSS6', 'TNFRSF1A', 'TNFRSF1B', 'TNN', 'TNXB', 'TPSAB1', 'TRIP6', 'TRPM7', 'TSC1', 'TUFT1', 'UTRN', 'VCAM1', 'VCAN', 'VCL', 'VEGFA', 'VTN', 'VWA2', 'WASHC1', 'WDPCP', 'WHAMM', 'WNT4', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnb.matrix.heatmap.pdf', width = 9, height = 6)

genes = c('ACVRL1', 'ADAMTS12', 'ADAMTS7', 'AMELX', 'ANXA6', 'AXIN2', 'BMP2', 'BMP4', 'BMP6', 'BMPR1A', 'BMPR1B', 'BMPR2', 'BPNT2', 'CCN2', 'CCN3', 'CCN4', 'CHADL', 'CHST11', 'CHSY1', 'COL11A1', 'COL27A1', 'COL2A1', 'COMP', 'CREB3L2', 'CYTL1', 'DDR2', 'ECM1', 'EFEMP1', 'EIF2AK3', 'EXT1', 'FGF18', 'FGFR3', 'GDF5', 'GDF6', 'GLG1', 'GLI2', 'GLI3', 'GPLD1', 'GREM1', 'HMGA2', 'HOXA11', 'IFT80', 'IHH', 'LNPK', 'LOXL2', 'LTBP3', 'LTF', 'MAF', 'MAPK14', 'MATN1', 'MBOAT2', 'MDK', 'MEF2C', 'MEF2D', 'MEX3C', 'MIR21', 'MMP14', 'MMP16', 'MUSTN1', 'NFIB', 'NKX3-2', 'NPPC', 'OSR1', 'OSR2', 'PKDCC', 'PTH', 'PTH1R', 'PTHLH', 'PTPN11', 'RARB', 'RARG', 'RFLNA', 'RFLNB', 'RUNX1', 'RUNX2', 'RUNX3', 'SCIN', 'SCX', 'SERPINH1', 'SFRP2', 'SHOX2', 'SIRT6', 'SIX2', 'SLC39A14', 'SMAD3', 'SMAD7', 'SMPD3', 'SNAI2', 'SNX19', 'SOX5', 'SOX6', 'SOX9', 'STC1', 'SULF1', 'SULF2', 'TGFB1', 'TGFBI', 'TGFBR1', 'TGFBR2', 'TRIP11', 'TRPS1', 'TSKU', 'WNT10B', 'WNT2B', 'WNT5B', 'WNT7A', 'WNT9A', 'ZBTB16', 'ZNF219')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnb.chondro.heatmap.pdf', width = 9, height = 6)

###cunn_flexcell_rnaseq_0321_ctl_flnc --- all samples minus bone and suture
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_flnc.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_flnc.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Gene + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321_ctl_flnc.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321_ctl_flnc.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnc.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnc.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnc.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnc.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_flnc.Cohort_pca.pdf', width=6, height = 6)

##clustering/heatmap from lists
genes = c('ABL1', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'AJUBA', 'APOD', 'ARHGAP6', 'ARHGEF7', 'BCAS3', 'BCL2', 'BCR', 'CAMSAP3', 'CLASP1', 'CLASP2', 'COL16A1', 'CORO1C', 'CORO2B', 'CTTN', 'DAPK3', 'DLC1', 'DMTN', 'DUSP22', 'DUSP3', 'EFNA5', 'EPB41L5', 'EPHA3', 'FAM107A', 'FERMT2', 'FMN1', 'GPM6B', 'GREM1', 'HRG', 'IQGAP1', 'ITGA2', 'ITGB1BP1', 'KDR', 'LDB1', 'LIMCH1', 'LIMS1', 'LRP1', 'MACF1', 'MAP4K4', 'MMP14', 'MYOC', 'NRP1', 'PDPK1', 'PEAK1', 'PHLDB2', 'PIP5K1A', 'PPM1F', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'RAC1', 'RCC2', 'RHOA', 'RHOD', 'ROCK1', 'ROCK2', 'S100A10', 'SDC4', 'SFRP1', 'SLC9A1', 'SLK', 'SMAD3', 'SORBS1', 'SRC', 'TAOK2', 'TEK', 'TESK2', 'THBS1', 'THSD1', 'THY1', 'TRIP6', 'TSC1', 'VCL', 'VEGFA', 'WDPCP', 'WHAMM', 'WNT4')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnc.focal_adhesion.heatmap.pdf', width = 9, height = 6)

genes = c('ABCB1', 'ACE', 'AGO3', 'ARIH2', 'ATXN1L', 'CCNE1', 'CD34', 'CITED1', 'CTC1', 'EIF2AK2', 'EPCAM', 'ETV6', 'FBLN1', 'FERMT1', 'FERMT2', 'FGF2', 'FZD1', 'GBA', 'GJA1', 'HMGA2', 'HMGB2', 'HNRNPU', 'KAT7', 'KDF1', 'KDM1A', 'KITLG', 'KRT18', 'LTBP3', 'MECOM', 'MIR16-1', 'MIR221', 'MIR222', 'MIR29B1', 'N4BP2L2', 'NANOG', 'NES', 'NF2', 'NR2E1', 'OVOL1', 'OVOL2', 'PDCD2', 'PDGFRA', 'PIM1', 'PTPRC', 'REST', 'RNF43', 'RUNX1', 'SFRP2', 'SIX2', 'SOX11', 'SOX17', 'SOX18', 'SOX5', 'SOX6', 'SOX9', 'TBX3', 'TERT', 'THPO', 'TRIM71', 'VEGFC', 'WNT1', 'WNT10B', 'WNT2B', 'WNT3', 'WNT5A', 'WNT7B', 'YAP1', 'YJEFN3', 'YTHDF2', 'ZFP36L1', 'ZNRF3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnc.msc_stem.heatmap.pdf', width = 9, height = 6)

genes = c('ABI2', 'ACTB', 'ADAM10', 'ADAM15', 'ADD1', 'AFDN', 'AHI1', 'AJAP1', 'AJM1', 'AJUBA', 'ALOX15B', 'ANG', 'ANXA1', 'ANXA2', 'APC', 'ARHGAP24', 'ARVCF', 'BAIAP2', 'BAIAP2L1', 'BMP6', 'BMPR2', 'CADM1', 'CADM2', 'CADM3', 'CAMSAP3', 'CCDC85A', 'CCDC85B', 'CCDC85C', 'CD99L2', 'CDC42', 'CDC42EP1', 'CDC42EP4', 'CDCA3', 'CDH1', 'CDH10', 'CDH11', 'CDH12', 'CDH13', 'CDH15', 'CDH17', 'CDH18', 'CDH19', 'CDH2', 'CDH20', 'CDH22', 'CDH24', 'CDH3', 'CDH4', 'CDH5', 'CDH6', 'CDH7', 'CDH8', 'CDH9', 'CDHR3', 'CEACAM1', 'CNN3', 'CRB1', 'CSK', 'CTNNA1', 'CTNNA2', 'CTNNA3', 'CTNNB1', 'CTNND1', 'CTNND2', 'CXADR', 'CYTH1', 'CYTH2', 'CYTH3', 'DAG1', 'DCHS1', 'DDX6', 'DLG1', 'DLG5', 'DLL1', 'DSC2', 'DSP', 'EFNB2', 'EIF4G2', 'EPB41L5', 'EPHA4', 'ESAM', 'FAT2', 'FERMT2', 'FLOT1', 'FLOT2', 'FMN1', 'FRMD4A', 'FRMD4B', 'FRMD5', 'FRS2', 'HIPK1', 'HMCN1', 'IGSF21', 'INAVA', 'ITGA6', 'JAG1', 'JAM3', 'JCAD', 'JUP', 'KIFC3', 'KLHL24', 'KRT18', 'LDB3', 'LIMD1', 'LIN7A', 'LIN7B', 'LIN7C', 'LYN', 'MAGI1', 'MPP4', 'MPP5', 'MPP7', 'MYH9', 'MYO1E', 'NDRG1', 'NECTIN1', 'NECTIN2', 'NECTIN3', 'NECTIN4', 'NEXN', 'NF2', 'NIBAN2', 'NOTCH1', 'NPHP1', 'NUMB', 'NUMBL', 'OXTR', 'PAK2', 'PAK4', 'PARD3', 'PARD3B', 'PARK7', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDLIM5', 'PDLIM7', 'PDZD11', 'PGM5', 'PIP5K1C', 'PKP1', 'PKP2', 'PKP3', 'PKP4', 'PLEKHA7', 'PLPP3', 'POF1B', 'PPP1CA', 'PPP1R9B', 'PRICKLE4', 'PTPN23', 'PTPRK', 'PTPRM', 'PVR', 'RAB10', 'RAMP2', 'RDX', 'RND1', 'S100A11', 'SCRIB', 'SDCBP', 'SH3BP1', 'SHROOM1', 'SHROOM2', 'SHROOM3', 'SHROOM4', 'SMAD7', 'SNAP23', 'SORBS1', 'SPTBN4', 'SRC', 'SSX2IP', 'STXBP6', 'SYNM', 'TBCD', 'TJP1', 'TJP2', 'TLN1', 'TMEM204', 'TMEM47', 'TMOD3', 'TNK2', 'TNKS1BP1', 'TRIM29', 'TRPV4', 'TSPAN33', 'VCL', 'VEGFA', 'VEZT', 'WNK3', 'WTIP', 'ZNF703', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnc.adherens.heatmap.pdf', width = 9, height = 6)

genes = c('ABL1', 'ACHE', 'ACVR1', 'ACVR2A', 'ACVR2B', 'ADAR', 'AKT1', 'ALPL', 'ALYREF', 'AMELX', 'ASF1A', 'ATF4', 'ATP5F1B', 'ATP6AP1', 'ATRAID', 'AXIN2', 'BAMBI', 'BCAP29', 'BCL2', 'BGLAP', 'BMP2', 'BMP3', 'BMP4', 'BMP6', 'BMP7', 'BMPR1A', 'BMPR1B', 'BMPR2', 'CAT', 'CBFB', 'CCDC47', 'CCL3', 'CCN1', 'CCN4', 'CDK6', 'CEBPA', 'CEBPB', 'CHRD', 'CITED1', 'CLEC5A', 'CLIC1', 'CLTC', 'COL1A1', 'COL6A1', 'CREB3L1', 'CRIM1', 'CTHRC1', 'CTNNBIP1', 'CYP24A1', 'DDR2', 'DDX21', 'DDX5', 'DHH', 'DHX9', 'DLX5', 'DNAI3', 'DNAJC13', 'EIF2AK2', 'EPHA2', 'FASN', 'FBL', 'FBN2', 'FBXO5', 'FERMT2', 'FFAR4', 'FGF23', 'FGFR2', 'FHL2', 'FIGNL1', 'FZD1', 'GATA1', 'GDF10', 'GDF2', 'GDPD2', 'GLI1', 'GLI2', 'GLI3', 'GNAS', 'GREM1', 'GTPBP4', 'H3-3A', 'H3-3B', 'HAND2', 'HDAC4', 'HDAC7', 'HEMGN', 'HGF', 'HNRNPC', 'HNRNPU', 'HOXA2', 'HPSE', 'HSD17B4', 'HSPE1', 'IARS1', 'IBSP', 'ID3', 'IFITM1', 'IFT80', 'IGF1', 'IGF2', 'IGFBP3', 'IGFBP5', 'IHH', 'IL6', 'IL6R', 'IL6ST', 'ITGA11', 'ITGAV', 'JAG1', 'JUNB', 'JUND', 'LEF1', 'LGR4', 'LIMD1', 'LOX', 'LRP3', 'LRP5', 'LRP5L', 'LTF', 'MEF2C', 'MEF2D', 'MEN1', 'MIR138-1', 'MIR200C', 'MIR208A', 'MIR20A', 'MIR20B', 'MIR21', 'MIR210', 'MIR27A', 'MIR29B1', 'MIR548D1', 'MIR665', 'MIR675', 'MIR9-1', 'MN1', 'MRC2', 'MSX2', 'MUS81', 'MYBBP1A', 'MYOC', 'NBR1', 'NELL1', 'NF1', 'NOCT', 'NOG', 'NOTCH1', 'NPNT', 'NPPC', 'NPR3', 'OSR2', 'OSTN', 'PDLIM7', 'PENK', 'PHB', 'PLXNB1', 'PPARG', 'PRKACA', 'PRKD1', 'PSMC2', 'PTCH1', 'PTH1R', 'PTHLH', 'PTK2', 'RANBP3L', 'RASSF2', 'RBMX', 'RDH14', 'REST', 'RHOA', 'RIOX1', 'RORB', 'RPS15', 'RRAS2', 'RRBP1', 'RSL1D1', 'RSPO2', 'RUNX2', 'SATB2', 'SEMA4D', 'SEMA7A', 'SFRP1', 'SFRP2', 'SHH', 'SHOX2', 'SIRT7', 'SKI', 'SMAD1', 'SMAD3', 'SMAD5', 'SMAD6', 'SMO', 'SMOC1', 'SNAI1', 'SNAI2', 'SND1', 'SNRNP200', 'SOX11', 'SOX2', 'SOX8', 'SOX9', 'SP7', 'SPP1', 'SUCO', 'SUFU', 'SYNCRIP', 'TCIRG1', 'TMEM119', 'TMEM64', 'TNC', 'TNF', 'TNN', 'TOB1', 'TP53INP2', 'TP63', 'TPM4', 'TRPM4', 'TWIST1', 'TWIST2', 'TWSG1', 'UCMA', 'UFL1', 'VCAN', 'VEGFC', 'WNT10B', 'WNT11', 'WNT3', 'WNT3A', 'WNT4', 'WNT7B', 'WWOX', 'WWTR1', 'YAP1', 'ZHX3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnc.osteo1.heatmap.pdf', width = 9, height = 6)

genes = c('A2M', 'ABL1', 'ACAN', 'ACER2', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'ADAM10', 'ADAM15', 'ADAM8', 'ADAM9', 'ADAMTS12', 'ADAMTS13', 'ADAMTS4', 'ADAMTS5', 'ADAMTS9', 'ADTRP', 'AEBP1', 'AGT', 'AJAP1', 'AJUBA', 'AMBN', 'AMELX', 'AMELY', 'ANGPTL3', 'ANGPTL7', 'ANTXR1', 'APOD', 'ARHGAP6', 'ARHGEF7', 'ASPN', 'ATP7A', 'BCAM', 'BCAS3', 'BCL2', 'BCL2L11', 'BCL6', 'BCR', 'BGN', 'BMP1', 'BMP2', 'BSG', 'BST1', 'CAMSAP3', 'CAPN1', 'CAPN2', 'CAPNS1', 'CAPNS2', 'CARMIL2', 'CASK', 'CCL21', 'CCL25', 'CCL28', 'CCN2', 'CCR7', 'CD34', 'CD36', 'CD3E', 'CD44', 'CD63', 'CD96', 'CDH13', 'CDK5', 'CDK6', 'CDKN2A', 'CEACAM6', 'CFLAR', 'CHADL', 'CIB1', 'CLASP1', 'CLASP2', 'CMA1', 'COL10A1', 'COL11A1', 'COL11A2', 'COL12A1', 'COL13A1', 'COL14A1', 'COL15A1', 'COL16A1', 'COL17A1', 'COL18A1', 'COL19A1', 'COL1A1', 'COL1A2', 'COL21A1', 'COL23A1', 'COL24A1', 'COL25A1', 'COL27A1', 'COL28A1', 'COL2A1', 'COL3A1', 'COL4A1', 'COL4A2', 'COL4A3', 'COL4A4', 'COL4A5', 'COL4A6', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2', 'COL6A3', 'COL6A5', 'COL6A6', 'COL7A1', 'COL8A1', 'COL8A2', 'COL9A1', 'COL9A2', 'COL9A3', 'COLGALT1', 'CORO1C', 'CORO2B', 'CPB2', 'CSF1', 'CST3', 'CTRB1', 'CTRB2', 'CTSG', 'CTSK', 'CTSL', 'CTSS', 'CTSV', 'CTTN', 'CX3CL1', 'DAG1', 'DAPK3', 'DCN', 'DDR1', 'DDR2', 'DEFB118', 'DISC1', 'DLC1', 'DMTN', 'DPP4', 'DUSP22', 'DUSP3', 'ECM2', 'EDA', 'EFEMP2', 'EFNA5', 'ELANE', 'ELN', 'EMILIN1', 'EMILIN2', 'EMILIN3', 'EMP2', 'ENAM', 'EPB41L5', 'EPDR1', 'EPHA1', 'EPHA3', 'ETS1', 'EXOC8', 'FAM107A', 'FAP', 'FBLN2', 'FBLN5', 'FBN1', 'FBN2', 'FERMT1', 'FERMT2', 'FERMT3', 'FGA', 'FGB', 'FGFR4', 'FGG', 'FKBP10', 'FLOT1', 'FMN1', 'FMOD', 'FN1', 'FREM1', 'FSCN1', 'FURIN', 'FUT1', 'GAS6', 'GFUS', 'GPM6B', 'GREM1', 'GSK3B', 'HAPLN1', 'HAPLN2', 'HAS1', 'HAS2', 'HAS3', 'HOXA7', 'HOXD3', 'HPN', 'HPSE', 'HRG', 'HSPG2', 'HTRA1', 'IHH', 'IL6', 'ILK', 'IQGAP1', 'ITGA1', 'ITGA10', 'ITGA11', 'ITGA2', 'ITGA2B', 'ITGA3', 'ITGA4', 'ITGA7', 'ITGA8', 'ITGAL', 'ITGAV', 'ITGB1', 'ITGB1BP1', 'ITGB2', 'ITGB3', 'ITGB4', 'ITGB5', 'ITGB6', 'ITGB7', 'ITGBL1', 'JAG1', 'JAM3', 'JUP', 'KDR', 'KIF9', 'KLK2', 'KLK4', 'KLK5', 'KLK7', 'KLKB1', 'L1CAM', 'LAMB1', 'LAMB2', 'LAMB3', 'LAMB4', 'LAMC1', 'LCP1', 'LDB1', 'LIMCH1', 'LIMS1', 'LOX', 'LRP1', 'LTBP3', 'LUM', 'LYPD3', 'LYPD5', 'LYVE1', 'MACF1', 'MADCAM1', 'MAP4K4', 'MELTF', 'MFAP4', 'MINK1', 'MIR192', 'MIR29B1', 'MIR29C', 'MIR92A1', 'MIR939', 'MIR98', 'MKLN1', 'MMP1', 'MMP10', 'MMP11', 'MMP12', 'MMP13', 'MMP14', 'MMP15', 'MMP16', 'MMP19', 'MMP2', 'MMP20', 'MMP3', 'MMP7', 'MMP8', 'MMP9', 'MSLN', 'MSLNL', 'MUC4', 'MYF5', 'MYH11', 'MYOC', 'NEXMIF', 'NF1', 'NF2', 'NID1', 'NID2', 'NOTCH1', 'NOXO1', 'NPNT', 'NRP1', 'NTN4', 'NTNG1', 'NTNG2', 'OGN', 'ONECUT1', 'ONECUT2', 'OTOA', 'PARVG', 'PDPK1', 'PDPN', 'PEAK1', 'PHLDB1', 'PHLDB2', 'PIK3CB', 'PIK3R1', 'PIP5K1A', 'PKD1', 'PKHD1', 'PLAU', 'PLEKHA2', 'PLET1', 'PLG', 'PLOD3', 'PODN', 'POSTN', 'PPFIA1', 'PPFIA2', 'PPM1F', 'PRELP', 'PRG2', 'PRG3', 'PRG4', 'PRKCZ', 'PRSS1', 'PRSS2', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'PXDN', 'PXN', 'QSOX1', 'RAC1', 'RAMP2', 'RASA1', 'RB1', 'RCC2', 'RGCC', 'RHOA', 'RHOD', 'RIC1', 'RIN2', 'ROCK1', 'ROCK2', 'RRAS', 'RUNX1', 'S100A10', 'SCUBE1', 'SCUBE3', 'SDC4', 'SEMA3E', 'SERPINE1', 'SFRP1', 'SGCE', 'SH3PXD2B', 'SIGLEC1', 'SKAP1', 'SLC2A10', 'SLC9A1', 'SLK', 'SMAD3', 'SMPD3', 'SNED1', 'SORBS1', 'SOX9', 'SRC', 'SRF', 'STATH', 'STRC', 'STRCP1', 'TAOK2', 'TCF15', 'TECTA', 'TEK', 'TESK2', 'TGFB1', 'THBS1', 'THBS3', 'THSD1', 'THSD4', 'THY1', 'TIAM1', 'TIE1', 'TIMM10B', 'TIMP1', 'TIMP2', 'TLL1', 'TLL2', 'TMEM8B', 'TMPRSS6', 'TNFRSF1A', 'TNFRSF1B', 'TNN', 'TNXB', 'TPSAB1', 'TRIP6', 'TRPM7', 'TSC1', 'TUFT1', 'UTRN', 'VCAM1', 'VCAN', 'VCL', 'VEGFA', 'VTN', 'VWA2', 'WASHC1', 'WDPCP', 'WHAMM', 'WNT4', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnc.matrix.heatmap.pdf', width = 9, height = 6)

genes = c('ACVRL1', 'ADAMTS12', 'ADAMTS7', 'AMELX', 'ANXA6', 'AXIN2', 'BMP2', 'BMP4', 'BMP6', 'BMPR1A', 'BMPR1B', 'BMPR2', 'BPNT2', 'CCN2', 'CCN3', 'CCN4', 'CHADL', 'CHST11', 'CHSY1', 'COL11A1', 'COL27A1', 'COL2A1', 'COMP', 'CREB3L2', 'CYTL1', 'DDR2', 'ECM1', 'EFEMP1', 'EIF2AK3', 'EXT1', 'FGF18', 'FGFR3', 'GDF5', 'GDF6', 'GLG1', 'GLI2', 'GLI3', 'GPLD1', 'GREM1', 'HMGA2', 'HOXA11', 'IFT80', 'IHH', 'LNPK', 'LOXL2', 'LTBP3', 'LTF', 'MAF', 'MAPK14', 'MATN1', 'MBOAT2', 'MDK', 'MEF2C', 'MEF2D', 'MEX3C', 'MIR21', 'MMP14', 'MMP16', 'MUSTN1', 'NFIB', 'NKX3-2', 'NPPC', 'OSR1', 'OSR2', 'PKDCC', 'PTH', 'PTH1R', 'PTHLH', 'PTPN11', 'RARB', 'RARG', 'RFLNA', 'RFLNB', 'RUNX1', 'RUNX2', 'RUNX3', 'SCIN', 'SCX', 'SERPINH1', 'SFRP2', 'SHOX2', 'SIRT6', 'SIX2', 'SLC39A14', 'SMAD3', 'SMAD7', 'SMPD3', 'SNAI2', 'SNX19', 'SOX5', 'SOX6', 'SOX9', 'STC1', 'SULF1', 'SULF2', 'TGFB1', 'TGFBI', 'TGFBR1', 'TGFBR2', 'TRIP11', 'TRPS1', 'TSKU', 'WNT10B', 'WNT2B', 'WNT5B', 'WNT7A', 'WNT9A', 'ZBTB16', 'ZNF219')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_flnc.chondro.heatmap.pdf', width = 9, height = 6)


###cunn_flexcell_rnaseq_0321_ctl_piezo1 --- all samples minus bone and suture
##read in count and metadata
countData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_piezo1.star_fc.counts.txt', header=T, row.names=1)
colData1 <- read.table('cunn_flexcell_rnaseq_0321_ctl_piezo1.star_fc.metadata.txt', header=T, row.names=1)
head(countData1)
head(colData1)
##convert column of numbers to a factor
colData1$Cohort <- as.factor(colData1$Cohort)

##add to deseq, give countdata and metadata and then design information i.e. info on sample types
dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, ~Sex + Gene + Condition)
dds

##remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 4, ]
nrow(dds)

##write normalized counts
dds_norm <- estimateSizeFactors(dds)
count_data <- counts(dds_norm, normalized=TRUE)
write.csv(count_data, file="cunn_flexcell_rnaseq_0321_ctl_piezo1.norm_counts.csv")
##vst transform data -- new version
rld <- vst(dds, blind=FALSE)
##check
head(assay(rld), 3)
head(assay(dds),3)
write.csv(assay(rld), file="cunn_flexcell_rnaseq_0321_ctl_piezo1.vst_counts.csv")
##principal components analysis 
plotPCA(rld, intgroup = c("sample_only"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_piezo1.sample_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Condition"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_piezo1.Condition_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Gene"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_piezo1.Gene_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Sex"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_piezo1.Sex_pca.pdf', width=6, height = 6)
plotPCA(rld, intgroup = c("Cohort"))
ggsave('cunn_flexcell_rnaseq_0321_ctl_piezo1.Cohort_pca.pdf', width=6, height = 6)

##clustering/heatmap from lists
genes = c('ABL1', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'AJUBA', 'APOD', 'ARHGAP6', 'ARHGEF7', 'BCAS3', 'BCL2', 'BCR', 'CAMSAP3', 'CLASP1', 'CLASP2', 'COL16A1', 'CORO1C', 'CORO2B', 'CTTN', 'DAPK3', 'DLC1', 'DMTN', 'DUSP22', 'DUSP3', 'EFNA5', 'EPB41L5', 'EPHA3', 'FAM107A', 'FERMT2', 'FMN1', 'GPM6B', 'GREM1', 'HRG', 'IQGAP1', 'ITGA2', 'ITGB1BP1', 'KDR', 'LDB1', 'LIMCH1', 'LIMS1', 'LRP1', 'MACF1', 'MAP4K4', 'MMP14', 'MYOC', 'NRP1', 'PDPK1', 'PEAK1', 'PHLDB2', 'PIP5K1A', 'PPM1F', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'RAC1', 'RCC2', 'RHOA', 'RHOD', 'ROCK1', 'ROCK2', 'S100A10', 'SDC4', 'SFRP1', 'SLC9A1', 'SLK', 'SMAD3', 'SORBS1', 'SRC', 'TAOK2', 'TEK', 'TESK2', 'THBS1', 'THSD1', 'THY1', 'TRIP6', 'TSC1', 'VCL', 'VEGFA', 'WDPCP', 'WHAMM', 'WNT4')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_piezo1.focal_adhesion.heatmap.pdf', width = 9, height = 6)

genes = c('ABCB1', 'ACE', 'AGO3', 'ARIH2', 'ATXN1L', 'CCNE1', 'CD34', 'CITED1', 'CTC1', 'EIF2AK2', 'EPCAM', 'ETV6', 'FBLN1', 'FERMT1', 'FERMT2', 'FGF2', 'FZD1', 'GBA', 'GJA1', 'HMGA2', 'HMGB2', 'HNRNPU', 'KAT7', 'KDF1', 'KDM1A', 'KITLG', 'KRT18', 'LTBP3', 'MECOM', 'MIR16-1', 'MIR221', 'MIR222', 'MIR29B1', 'N4BP2L2', 'NANOG', 'NES', 'NF2', 'NR2E1', 'OVOL1', 'OVOL2', 'PDCD2', 'PDGFRA', 'PIM1', 'PTPRC', 'REST', 'RNF43', 'RUNX1', 'SFRP2', 'SIX2', 'SOX11', 'SOX17', 'SOX18', 'SOX5', 'SOX6', 'SOX9', 'TBX3', 'TERT', 'THPO', 'TRIM71', 'VEGFC', 'WNT1', 'WNT10B', 'WNT2B', 'WNT3', 'WNT5A', 'WNT7B', 'YAP1', 'YJEFN3', 'YTHDF2', 'ZFP36L1', 'ZNRF3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_piezo1.msc_stem.heatmap.pdf', width = 9, height = 6)

genes = c('ABI2', 'ACTB', 'ADAM10', 'ADAM15', 'ADD1', 'AFDN', 'AHI1', 'AJAP1', 'AJM1', 'AJUBA', 'ALOX15B', 'ANG', 'ANXA1', 'ANXA2', 'APC', 'ARHGAP24', 'ARVCF', 'BAIAP2', 'BAIAP2L1', 'BMP6', 'BMPR2', 'CADM1', 'CADM2', 'CADM3', 'CAMSAP3', 'CCDC85A', 'CCDC85B', 'CCDC85C', 'CD99L2', 'CDC42', 'CDC42EP1', 'CDC42EP4', 'CDCA3', 'CDH1', 'CDH10', 'CDH11', 'CDH12', 'CDH13', 'CDH15', 'CDH17', 'CDH18', 'CDH19', 'CDH2', 'CDH20', 'CDH22', 'CDH24', 'CDH3', 'CDH4', 'CDH5', 'CDH6', 'CDH7', 'CDH8', 'CDH9', 'CDHR3', 'CEACAM1', 'CNN3', 'CRB1', 'CSK', 'CTNNA1', 'CTNNA2', 'CTNNA3', 'CTNNB1', 'CTNND1', 'CTNND2', 'CXADR', 'CYTH1', 'CYTH2', 'CYTH3', 'DAG1', 'DCHS1', 'DDX6', 'DLG1', 'DLG5', 'DLL1', 'DSC2', 'DSP', 'EFNB2', 'EIF4G2', 'EPB41L5', 'EPHA4', 'ESAM', 'FAT2', 'FERMT2', 'FLOT1', 'FLOT2', 'FMN1', 'FRMD4A', 'FRMD4B', 'FRMD5', 'FRS2', 'HIPK1', 'HMCN1', 'IGSF21', 'INAVA', 'ITGA6', 'JAG1', 'JAM3', 'JCAD', 'JUP', 'KIFC3', 'KLHL24', 'KRT18', 'LDB3', 'LIMD1', 'LIN7A', 'LIN7B', 'LIN7C', 'LYN', 'MAGI1', 'MPP4', 'MPP5', 'MPP7', 'MYH9', 'MYO1E', 'NDRG1', 'NECTIN1', 'NECTIN2', 'NECTIN3', 'NECTIN4', 'NEXN', 'NF2', 'NIBAN2', 'NOTCH1', 'NPHP1', 'NUMB', 'NUMBL', 'OXTR', 'PAK2', 'PAK4', 'PARD3', 'PARD3B', 'PARK7', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDLIM5', 'PDLIM7', 'PDZD11', 'PGM5', 'PIP5K1C', 'PKP1', 'PKP2', 'PKP3', 'PKP4', 'PLEKHA7', 'PLPP3', 'POF1B', 'PPP1CA', 'PPP1R9B', 'PRICKLE4', 'PTPN23', 'PTPRK', 'PTPRM', 'PVR', 'RAB10', 'RAMP2', 'RDX', 'RND1', 'S100A11', 'SCRIB', 'SDCBP', 'SH3BP1', 'SHROOM1', 'SHROOM2', 'SHROOM3', 'SHROOM4', 'SMAD7', 'SNAP23', 'SORBS1', 'SPTBN4', 'SRC', 'SSX2IP', 'STXBP6', 'SYNM', 'TBCD', 'TJP1', 'TJP2', 'TLN1', 'TMEM204', 'TMEM47', 'TMOD3', 'TNK2', 'TNKS1BP1', 'TRIM29', 'TRPV4', 'TSPAN33', 'VCL', 'VEGFA', 'VEZT', 'WNK3', 'WTIP', 'ZNF703', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_piezo1.adherens.heatmap.pdf', width = 9, height = 6)

genes = c('ABL1', 'ACHE', 'ACVR1', 'ACVR2A', 'ACVR2B', 'ADAR', 'AKT1', 'ALPL', 'ALYREF', 'AMELX', 'ASF1A', 'ATF4', 'ATP5F1B', 'ATP6AP1', 'ATRAID', 'AXIN2', 'BAMBI', 'BCAP29', 'BCL2', 'BGLAP', 'BMP2', 'BMP3', 'BMP4', 'BMP6', 'BMP7', 'BMPR1A', 'BMPR1B', 'BMPR2', 'CAT', 'CBFB', 'CCDC47', 'CCL3', 'CCN1', 'CCN4', 'CDK6', 'CEBPA', 'CEBPB', 'CHRD', 'CITED1', 'CLEC5A', 'CLIC1', 'CLTC', 'COL1A1', 'COL6A1', 'CREB3L1', 'CRIM1', 'CTHRC1', 'CTNNBIP1', 'CYP24A1', 'DDR2', 'DDX21', 'DDX5', 'DHH', 'DHX9', 'DLX5', 'DNAI3', 'DNAJC13', 'EIF2AK2', 'EPHA2', 'FASN', 'FBL', 'FBN2', 'FBXO5', 'FERMT2', 'FFAR4', 'FGF23', 'FGFR2', 'FHL2', 'FIGNL1', 'FZD1', 'GATA1', 'GDF10', 'GDF2', 'GDPD2', 'GLI1', 'GLI2', 'GLI3', 'GNAS', 'GREM1', 'GTPBP4', 'H3-3A', 'H3-3B', 'HAND2', 'HDAC4', 'HDAC7', 'HEMGN', 'HGF', 'HNRNPC', 'HNRNPU', 'HOXA2', 'HPSE', 'HSD17B4', 'HSPE1', 'IARS1', 'IBSP', 'ID3', 'IFITM1', 'IFT80', 'IGF1', 'IGF2', 'IGFBP3', 'IGFBP5', 'IHH', 'IL6', 'IL6R', 'IL6ST', 'ITGA11', 'ITGAV', 'JAG1', 'JUNB', 'JUND', 'LEF1', 'LGR4', 'LIMD1', 'LOX', 'LRP3', 'LRP5', 'LRP5L', 'LTF', 'MEF2C', 'MEF2D', 'MEN1', 'MIR138-1', 'MIR200C', 'MIR208A', 'MIR20A', 'MIR20B', 'MIR21', 'MIR210', 'MIR27A', 'MIR29B1', 'MIR548D1', 'MIR665', 'MIR675', 'MIR9-1', 'MN1', 'MRC2', 'MSX2', 'MUS81', 'MYBBP1A', 'MYOC', 'NBR1', 'NELL1', 'NF1', 'NOCT', 'NOG', 'NOTCH1', 'NPNT', 'NPPC', 'NPR3', 'OSR2', 'OSTN', 'PDLIM7', 'PENK', 'PHB', 'PLXNB1', 'PPARG', 'PRKACA', 'PRKD1', 'PSMC2', 'PTCH1', 'PTH1R', 'PTHLH', 'PTK2', 'RANBP3L', 'RASSF2', 'RBMX', 'RDH14', 'REST', 'RHOA', 'RIOX1', 'RORB', 'RPS15', 'RRAS2', 'RRBP1', 'RSL1D1', 'RSPO2', 'RUNX2', 'SATB2', 'SEMA4D', 'SEMA7A', 'SFRP1', 'SFRP2', 'SHH', 'SHOX2', 'SIRT7', 'SKI', 'SMAD1', 'SMAD3', 'SMAD5', 'SMAD6', 'SMO', 'SMOC1', 'SNAI1', 'SNAI2', 'SND1', 'SNRNP200', 'SOX11', 'SOX2', 'SOX8', 'SOX9', 'SP7', 'SPP1', 'SUCO', 'SUFU', 'SYNCRIP', 'TCIRG1', 'TMEM119', 'TMEM64', 'TNC', 'TNF', 'TNN', 'TOB1', 'TP53INP2', 'TP63', 'TPM4', 'TRPM4', 'TWIST1', 'TWIST2', 'TWSG1', 'UCMA', 'UFL1', 'VCAN', 'VEGFC', 'WNT10B', 'WNT11', 'WNT3', 'WNT3A', 'WNT4', 'WNT7B', 'WWOX', 'WWTR1', 'YAP1', 'ZHX3')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_piezo1.osteo1.heatmap.pdf', width = 9, height = 6)

genes = c('A2M', 'ABL1', 'ACAN', 'ACER2', 'ACTG1', 'ACTN1', 'ACTN2', 'ACTN3', 'ACVRL1', 'ADAM10', 'ADAM15', 'ADAM8', 'ADAM9', 'ADAMTS12', 'ADAMTS13', 'ADAMTS4', 'ADAMTS5', 'ADAMTS9', 'ADTRP', 'AEBP1', 'AGT', 'AJAP1', 'AJUBA', 'AMBN', 'AMELX', 'AMELY', 'ANGPTL3', 'ANGPTL7', 'ANTXR1', 'APOD', 'ARHGAP6', 'ARHGEF7', 'ASPN', 'ATP7A', 'BCAM', 'BCAS3', 'BCL2', 'BCL2L11', 'BCL6', 'BCR', 'BGN', 'BMP1', 'BMP2', 'BSG', 'BST1', 'CAMSAP3', 'CAPN1', 'CAPN2', 'CAPNS1', 'CAPNS2', 'CARMIL2', 'CASK', 'CCL21', 'CCL25', 'CCL28', 'CCN2', 'CCR7', 'CD34', 'CD36', 'CD3E', 'CD44', 'CD63', 'CD96', 'CDH13', 'CDK5', 'CDK6', 'CDKN2A', 'CEACAM6', 'CFLAR', 'CHADL', 'CIB1', 'CLASP1', 'CLASP2', 'CMA1', 'COL10A1', 'COL11A1', 'COL11A2', 'COL12A1', 'COL13A1', 'COL14A1', 'COL15A1', 'COL16A1', 'COL17A1', 'COL18A1', 'COL19A1', 'COL1A1', 'COL1A2', 'COL21A1', 'COL23A1', 'COL24A1', 'COL25A1', 'COL27A1', 'COL28A1', 'COL2A1', 'COL3A1', 'COL4A1', 'COL4A2', 'COL4A3', 'COL4A4', 'COL4A5', 'COL4A6', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2', 'COL6A3', 'COL6A5', 'COL6A6', 'COL7A1', 'COL8A1', 'COL8A2', 'COL9A1', 'COL9A2', 'COL9A3', 'COLGALT1', 'CORO1C', 'CORO2B', 'CPB2', 'CSF1', 'CST3', 'CTRB1', 'CTRB2', 'CTSG', 'CTSK', 'CTSL', 'CTSS', 'CTSV', 'CTTN', 'CX3CL1', 'DAG1', 'DAPK3', 'DCN', 'DDR1', 'DDR2', 'DEFB118', 'DISC1', 'DLC1', 'DMTN', 'DPP4', 'DUSP22', 'DUSP3', 'ECM2', 'EDA', 'EFEMP2', 'EFNA5', 'ELANE', 'ELN', 'EMILIN1', 'EMILIN2', 'EMILIN3', 'EMP2', 'ENAM', 'EPB41L5', 'EPDR1', 'EPHA1', 'EPHA3', 'ETS1', 'EXOC8', 'FAM107A', 'FAP', 'FBLN2', 'FBLN5', 'FBN1', 'FBN2', 'FERMT1', 'FERMT2', 'FERMT3', 'FGA', 'FGB', 'FGFR4', 'FGG', 'FKBP10', 'FLOT1', 'FMN1', 'FMOD', 'FN1', 'FREM1', 'FSCN1', 'FURIN', 'FUT1', 'GAS6', 'GFUS', 'GPM6B', 'GREM1', 'GSK3B', 'HAPLN1', 'HAPLN2', 'HAS1', 'HAS2', 'HAS3', 'HOXA7', 'HOXD3', 'HPN', 'HPSE', 'HRG', 'HSPG2', 'HTRA1', 'IHH', 'IL6', 'ILK', 'IQGAP1', 'ITGA1', 'ITGA10', 'ITGA11', 'ITGA2', 'ITGA2B', 'ITGA3', 'ITGA4', 'ITGA7', 'ITGA8', 'ITGAL', 'ITGAV', 'ITGB1', 'ITGB1BP1', 'ITGB2', 'ITGB3', 'ITGB4', 'ITGB5', 'ITGB6', 'ITGB7', 'ITGBL1', 'JAG1', 'JAM3', 'JUP', 'KDR', 'KIF9', 'KLK2', 'KLK4', 'KLK5', 'KLK7', 'KLKB1', 'L1CAM', 'LAMB1', 'LAMB2', 'LAMB3', 'LAMB4', 'LAMC1', 'LCP1', 'LDB1', 'LIMCH1', 'LIMS1', 'LOX', 'LRP1', 'LTBP3', 'LUM', 'LYPD3', 'LYPD5', 'LYVE1', 'MACF1', 'MADCAM1', 'MAP4K4', 'MELTF', 'MFAP4', 'MINK1', 'MIR192', 'MIR29B1', 'MIR29C', 'MIR92A1', 'MIR939', 'MIR98', 'MKLN1', 'MMP1', 'MMP10', 'MMP11', 'MMP12', 'MMP13', 'MMP14', 'MMP15', 'MMP16', 'MMP19', 'MMP2', 'MMP20', 'MMP3', 'MMP7', 'MMP8', 'MMP9', 'MSLN', 'MSLNL', 'MUC4', 'MYF5', 'MYH11', 'MYOC', 'NEXMIF', 'NF1', 'NF2', 'NID1', 'NID2', 'NOTCH1', 'NOXO1', 'NPNT', 'NRP1', 'NTN4', 'NTNG1', 'NTNG2', 'OGN', 'ONECUT1', 'ONECUT2', 'OTOA', 'PARVG', 'PDPK1', 'PDPN', 'PEAK1', 'PHLDB1', 'PHLDB2', 'PIK3CB', 'PIK3R1', 'PIP5K1A', 'PKD1', 'PKHD1', 'PLAU', 'PLEKHA2', 'PLET1', 'PLG', 'PLOD3', 'PODN', 'POSTN', 'PPFIA1', 'PPFIA2', 'PPM1F', 'PRELP', 'PRG2', 'PRG3', 'PRG4', 'PRKCZ', 'PRSS1', 'PRSS2', 'PTEN', 'PTK2', 'PTK2B', 'PTPRA', 'PTPRJ', 'PTPRK', 'PXDN', 'PXN', 'QSOX1', 'RAC1', 'RAMP2', 'RASA1', 'RB1', 'RCC2', 'RGCC', 'RHOA', 'RHOD', 'RIC1', 'RIN2', 'ROCK1', 'ROCK2', 'RRAS', 'RUNX1', 'S100A10', 'SCUBE1', 'SCUBE3', 'SDC4', 'SEMA3E', 'SERPINE1', 'SFRP1', 'SGCE', 'SH3PXD2B', 'SIGLEC1', 'SKAP1', 'SLC2A10', 'SLC9A1', 'SLK', 'SMAD3', 'SMPD3', 'SNED1', 'SORBS1', 'SOX9', 'SRC', 'SRF', 'STATH', 'STRC', 'STRCP1', 'TAOK2', 'TCF15', 'TECTA', 'TEK', 'TESK2', 'TGFB1', 'THBS1', 'THBS3', 'THSD1', 'THSD4', 'THY1', 'TIAM1', 'TIE1', 'TIMM10B', 'TIMP1', 'TIMP2', 'TLL1', 'TLL2', 'TMEM8B', 'TMPRSS6', 'TNFRSF1A', 'TNFRSF1B', 'TNN', 'TNXB', 'TPSAB1', 'TRIP6', 'TRPM7', 'TSC1', 'TUFT1', 'UTRN', 'VCAM1', 'VCAN', 'VCL', 'VEGFA', 'VTN', 'VWA2', 'WASHC1', 'WDPCP', 'WHAMM', 'WNT4', 'ZYX')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=3, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_piezo1.matrix.heatmap.pdf', width = 9, height = 6)

genes = c('ACVRL1', 'ADAMTS12', 'ADAMTS7', 'AMELX', 'ANXA6', 'AXIN2', 'BMP2', 'BMP4', 'BMP6', 'BMPR1A', 'BMPR1B', 'BMPR2', 'BPNT2', 'CCN2', 'CCN3', 'CCN4', 'CHADL', 'CHST11', 'CHSY1', 'COL11A1', 'COL27A1', 'COL2A1', 'COMP', 'CREB3L2', 'CYTL1', 'DDR2', 'ECM1', 'EFEMP1', 'EIF2AK3', 'EXT1', 'FGF18', 'FGFR3', 'GDF5', 'GDF6', 'GLG1', 'GLI2', 'GLI3', 'GPLD1', 'GREM1', 'HMGA2', 'HOXA11', 'IFT80', 'IHH', 'LNPK', 'LOXL2', 'LTBP3', 'LTF', 'MAF', 'MAPK14', 'MATN1', 'MBOAT2', 'MDK', 'MEF2C', 'MEF2D', 'MEX3C', 'MIR21', 'MMP14', 'MMP16', 'MUSTN1', 'NFIB', 'NKX3-2', 'NPPC', 'OSR1', 'OSR2', 'PKDCC', 'PTH', 'PTH1R', 'PTHLH', 'PTPN11', 'RARB', 'RARG', 'RFLNA', 'RFLNB', 'RUNX1', 'RUNX2', 'RUNX3', 'SCIN', 'SCX', 'SERPINH1', 'SFRP2', 'SHOX2', 'SIRT6', 'SIX2', 'SLC39A14', 'SMAD3', 'SMAD7', 'SMPD3', 'SNAI2', 'SNX19', 'SOX5', 'SOX6', 'SOX9', 'STC1', 'SULF1', 'SULF2', 'TGFB1', 'TGFBI', 'TGFBR1', 'TGFBR2', 'TRIP11', 'TRPS1', 'TSKU', 'WNT10B', 'WNT2B', 'WNT5B', 'WNT7A', 'WNT9A', 'ZBTB16', 'ZNF219')
foo <- assay(rld) #save assay as matrix
assaygenes <- row.names(foo) #extract rownames
idx <- assaygenes %in% genes #find target genes
newmat <- assay(rld)[idx,] #subset to target genes
newmat <- newmat - rowMeans(newmat)
newdf <- as.data.frame(colData(rld)[c("Condition", "Gene", "Sex", "Cohort")])
pheatmap(newmat, annotation_col=newdf, fontsize_row=4, fontsize_col=8)
dev.copy2pdf(file='cunn_flexcell_rnaseq_0321_ctl_piezo1.chondro.heatmap.pdf', width = 9, height = 6)



