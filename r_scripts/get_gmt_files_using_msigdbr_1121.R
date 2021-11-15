install.packages("msigdbr")
library(msigdbr)


kegg_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")

msigdbr_list = split(x = kegg_gene_sets$gene_symbol, f = kegg_gene_sets$gs_name)
write.table(msigdbr_list, file = 'Mm.kegg.1121.gmt', row.names=FALSE)

capture.output(summary(msigdbr_list), file = 'Mm.kegg.1121.gmt')

cat(capture.output(print(msigdbr_list), file="Mm.kegg.1121.gmt"))
getwd()

capture.output(msigdbr_list, file = "Mm.kegg.1121.gmt")

