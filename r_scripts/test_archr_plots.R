library(reshape2)
melted_hm <- melt(heatmap.matrix)
head(melted_hm)

rownames(heatmap.matrix)
colnames(heatmap.matrix)

p1 <- ggplot(data = melted_hm, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
p1

pheatmap(heatmap.matrix