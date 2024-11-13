# Trying to replicate the umap plot with a poplar genome

embeddings = read.csv("GCF_000001735.4_TAIR10.1._avg_embeddings.tsv", sep ="\t", head=F)
embeddings.umap = umap(embeddings)


plot(embeddings.umap $layout[,1], embeddings.umap $layout[,2], cex=.4, pch=19, xlim=c(-10,10), ylim=c(-10,10))



