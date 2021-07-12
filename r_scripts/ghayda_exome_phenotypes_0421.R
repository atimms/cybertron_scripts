##load libraries
library("pheatmap")
library("ggplot2")
library("FactoMineR")
library("factoextra")
library(ggpubr)

workingDir = "/archive/mirzaa_g/exomes/result_files_0321/phenotype_graphing_0421";
setwd(workingDir);

##heatmaps.....

##read in data
phenotype_data <- read.table('phenotypes_042821.txt', header=T, row.names=1, sep='\t')
##remove dx col
phenotype_data_edit <- subset (phenotype_data, select = -DxGroup1)
##make df for annotation and add rownames
pheno_df = data.frame("Dx" = phenotype_data$DxGroup1)
rownames(pheno_df) <- rownames(phenotype_data)
##graph and make pdf
pheatmap(phenotype_data_edit, annotation_row = pheno_df, fontsize_row=2)
dev.copy2pdf(file='phenotype_heatmap_042821.pdf', height = 12, width = 6 )
##cut data
out = pheatmap(phenotype_data_edit, annotation_row = pheno_df, fontsize_row=2,cutree_rows = 5)
dev.copy2pdf(file='phenotype_heatmap_cut5_042821.pdf', height = 12, width = 6 )
##get group names (need to double check if we're using this)
ped_groups = sort(cutree(out$tree_row, k=5))
head(ped_groups, 30)

##MCA analysis....

#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/114-mca-multiple-correspondence-analysis-in-r-essentials/

##read in data
phenotype_data <- read.table('phenotypes_adj_042921.txt', header=T, row.names=1, sep='\t')
##remove dx col and sumarize -- not needed
#phenotype_data_active <- subset(phenotype_data, select = -DxGroup1)
summary(phenotype_data)

##run mca
#res.mca <- MCA(phenotype_data_active, graph = FALSE)
res.mca <- MCA(phenotype_data, quali.sup = 1,  graph=FALSE)
res.mca
## visualize the percentages of inertia explained by each MCA dimensions
scree.plot = fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))
## biplot of individuals and variable categories
biplot = fviz_mca_biplot(res.mca, repel = F,  ggtheme = theme_minimal())
## extract the results for variable categories and view
var <- get_mca_var(res.mca)
var
head(var$coord)# Coordinates
head(var$cos2)# Cos2: quality on the factore map
head(var$contrib)# Contributions to the principal components
##correlation between variables and MCA principal dimensions
mca.var =fviz_mca_var(res.mca, choice = "mca.cor", 
             repel = TRUE, # Avoid text overlapping (slow)
             ggtheme = theme_minimal())

fviz_cos2(res.mca, choice = "var", axes = 1:2)
# Contributions of rows to different dimensions
cont1 = fviz_contrib(res.mca, choice = "var", axes = 1, top = 15)
cont2 = fviz_contrib(res.mca, choice = "var", axes = 2, top = 15)
cont3 = fviz_contrib(res.mca, choice = "var", axes = 3, top = 15)
##indivdual with contribution
ind_cont = fviz_mca_ind(res.mca, col.ind = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             ggtheme = theme_minimal())
##label individuals
ind_dx = fviz_mca_ind(res.mca, habillage = 1, label = "none")
ind_wm = fviz_mca_ind(res.mca, habillage = 2, addEllipses = TRUE, ellipse.type = "confidence", label = "none")
ind_ctx = fviz_mca_ind(res.mca, habillage = 3, addEllipses = TRUE, ellipse.type = "confidence", label = "none")
ind_ven = fviz_mca_ind(res.mca, habillage = 4, addEllipses = TRUE, ellipse.type = "confidence", label = "none")
ind_cbl = fviz_mca_ind(res.mca, habillage = 5, addEllipses = TRUE, ellipse.type = "confidence", label = "none")
ind_bs = fviz_mca_ind(res.mca, habillage = 6, addEllipses = TRUE, ellipse.type = "confidence", label = "none")
ind_bgth = fviz_mca_ind(res.mca, habillage = 7, addEllipses = TRUE, ellipse.type = "confidence", label = "none")
ind_cc = fviz_mca_ind(res.mca, habillage = 8, addEllipses = TRUE, ellipse.type = "confidence", label = "none")


##write pdfs
ggexport(plotlist = list(scree.plot, biplot, mca.var, cont1, cont2, cont3, ind_cont, ind_dx, ind_wm, ind_ctx, ind_ven, ind_cbl, ind_bs, ind_bgth, ind_cc), 
         filename = "MCA_analysis_042921.pdf")

