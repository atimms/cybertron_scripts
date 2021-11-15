# install.packages("dendsort")

##load libraries
library("pheatmap")
library("ggplot2")
library('dendsort')

workingDir = "/active/cherry_t/OrgManuscript_SingleCell_Data/mouse_retina_data_0521/human_mouse_comparisons_0621";
setwd(workingDir);


##graph human data
human_hm <- read.table('human_0621.ret_genes.heatmap_all.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm, cluster_rows = F, show_rownames = F)
dev.copy2pdf(file='human_0621.ret_genes.heatmap_all.no_cluster.pdf', width = 9, height = 6)
pheatmap(human_hm, show_rownames = F)
dev.copy2pdf(file='human_0621.ret_genes.heatmap_all.cluster.pdf', width = 9, height = 6)

human_hm <- read.table('human_0621.ret_genes.heatmap_markers.txt', row.names=1, header=T, sep='\t')
pheatmap(human_hm, cluster_rows = F, show_rownames = F)
dev.copy2pdf(file='human_0621.ret_genes.heatmap_markers.no_cluster.pdf', width = 9, height = 6)
pheatmap(human_hm, show_rownames = F)
dev.copy2pdf(file='human_0621.ret_genes.heatmap_markers.cluster.pdf', width = 9, height = 6)

##graph organoid data
organoid_hm <- read.table('organoid_0621.ret_genes.heatmap_all.txt', row.names=1, header=T, sep='\t')
pheatmap(organoid_hm, cluster_rows = F, show_rownames = F)
dev.copy2pdf(file='organoid_0621.ret_genes.heatmap_all.no_cluster.pdf', width = 9, height = 6)
pheatmap(organoid_hm, show_rownames = F)
dev.copy2pdf(file='organoid_0621.ret_genes.heatmap_all.cluster.pdf', width = 9, height = 6)

organoid_hm <- read.table('organoid_0621.ret_genes.heatmap_markers.txt', row.names=1, header=T, sep='\t')
pheatmap(organoid_hm, cluster_rows = F, show_rownames = F)
dev.copy2pdf(file='organoid_0621.ret_genes.heatmap_markers.no_cluster.pdf', width = 9, height = 6)
pheatmap(organoid_hm, show_rownames = F)
dev.copy2pdf(file='organoid_0621.ret_genes.heatmap_markers.cluster.pdf', width = 9, height = 6)

##graph mouse data
mouse_hm <- read.table('mouse_0621.ret_genes.heatmap_all.txt', row.names=1, header=T, sep='\t')
pheatmap(mouse_hm, cluster_rows = F, show_rownames = F)
dev.copy2pdf(file='mouse_0621.ret_genes.heatmap_all.no_cluster.pdf', width = 9, height = 6)
pheatmap(mouse_hm, show_rownames = F)
dev.copy2pdf(file='mouse_0621.ret_genes.heatmap_all.cluster.pdf', width = 9, height = 6)

mouse_hm <- read.table('mouse_0621.ret_genes.heatmap_markers.txt', row.names=1, header=T, sep='\t')
pheatmap(mouse_hm, cluster_rows = F, show_rownames = F)
dev.copy2pdf(file='mouse_0621.ret_genes.heatmap_markers.no_cluster.pdf', width = 9, height = 6)
pheatmap(mouse_hm, show_rownames = F)
dev.copy2pdf(file='mouse_0621.ret_genes.heatmap_markers.cluster.pdf', width = 9, height = 6)

##function for sort clusters.. https://slowkow.com/notes/pheatmap-tutorial/
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

##graph human/organoid data
ho_hm <- read.table('human_organoid_0621.ret_genes.heatmap_all.int.txt', row.names=1, header=T, sep='\t')
pheatmap(ho_hm, show_rownames = F)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_all.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(ho_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(ho_hm)))
pheatmap(ho_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_all.cluster2.pdf', width = 9, height = 6)

ho_hm <- read.table('human_organoid_0621.ret_genes.heatmap_markers.int.txt', row.names=1, header=T, sep='\t')
pheatmap(ho_hm, show_rownames = F)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_markers.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(ho_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(ho_hm)))
pheatmap(ho_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_markers.cluster2.pdf', width = 9, height = 6)

##graph human/mouse data
hm_hm <- read.table('human_mouse_0621.ret_genes.heatmap_all.int.txt', row.names=1, header=T, sep='\t')
pheatmap(hm_hm, show_rownames = F)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_all.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(hm_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(hm_hm)))
pheatmap(hm_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_all.cluster2.pdf', width = 9, height = 6)

hm_hm <- read.table('human_mouse_0621.ret_genes.heatmap_markers.int.txt', row.names=1, header=T, sep='\t')
pheatmap(hm_hm, show_rownames = F)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_markers.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(hm_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(hm_hm)))
pheatmap(hm_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_markers.cluster2.pdf', width = 9, height = 6)

##graph human/mouse data -- removing hu_Mature Ganglions
hm_hm <- read.table('human_mouse_0621.ret_genes.heatmap_all.int.txt', row.names=1, header=T, sep='\t')
hm_hm = subset(hm_hm, select=-c(hu_Mature.Ganglions))
pheatmap(hm_hm, show_rownames = F)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_all.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(hm_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(hm_hm)))
pheatmap(hm_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_all.no_RGC.cluster2.pdf', width = 9, height = 6)
hm_hm <- read.table('human_mouse_0621.ret_genes.heatmap_markers.int.txt', row.names=1, header=T, sep='\t')
hm_hm = subset(hm_hm, select=-c(hu_Mature.Ganglions))
pheatmap(hm_hm, show_rownames = F)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_markers.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(hm_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(hm_hm)))
pheatmap(hm_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_mouse_0621.ret_genes.heatmap_markers.no_RGC.cluster2.pdf', width = 9, height = 6)

##graph human/organoid data -- removing hu_Mature Ganglions
ho_hm <- read.table('human_organoid_0621.ret_genes.heatmap_all.int.txt', row.names=1, header=T, sep='\t')
ho_hm = subset(ho_hm, select=-c(hu_Mature.Ganglions))
pheatmap(ho_hm, show_rownames = F)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_all.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(ho_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(ho_hm)))
pheatmap(ho_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_all.no_RGC.cluster2.pdf', width = 9, height = 6)

ho_hm <- read.table('human_organoid_0621.ret_genes.heatmap_markers.int.txt', row.names=1, header=T, sep='\t')
ho_hm = subset(ho_hm, select=-c(hu_Mature.Ganglions))
pheatmap(ho_hm, show_rownames = F)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_markers.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(ho_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(ho_hm)))
pheatmap(ho_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_org_0621.ret_genes.heatmap_markers.no_RGC.cluster2.pdf', width = 9, height = 6)


##graph human/organoid data -- removing hu_Mature Ganglions... with new data 1021
ho_hm <- read.table('human_organoid_1021.ret_genes.heatmap_all.int.txt', row.names=1, header=T, sep='\t')
ho_hm = subset(ho_hm, select=-c(hu_Mature.Ganglions))
pheatmap(ho_hm, show_rownames = F)
dev.copy2pdf(file='human_org_1021.ret_genes.heatmap_all.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(ho_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(ho_hm)))
pheatmap(ho_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_org_1021.ret_genes.heatmap_all.no_RGC.cluster2.pdf', width = 9, height = 6)

ho_hm <- read.table('human_organoid_1021.ret_genes.heatmap_markers.int.txt', row.names=1, header=T, sep='\t')
ho_hm = subset(ho_hm, select=-c(hu_Mature.Ganglions))
pheatmap(ho_hm, show_rownames = F)
dev.copy2pdf(file='human_org_1021.ret_genes.heatmap_markers.no_RGC.cluster.pdf', width = 9, height = 6)
mat_cluster_cols <- hclust(dist(t(ho_hm)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(ho_hm)))
pheatmap(ho_hm, show_rownames = F, cluster_cols = mat_cluster_cols, cluster_rows = mat_cluster_rows)
dev.copy2pdf(file='human_org_1021.ret_genes.heatmap_markers.no_RGC.cluster2.pdf', width = 9, height = 6)

