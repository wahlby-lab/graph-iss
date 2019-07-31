library(Seurat)
library(dplyr)

exp_tab <- read.table('../data/results/expression.csv', sep=",",stringsAsFactors=F,header=T,row.names = 1)
norm_exp_tab <-  read.table('../data/results/residuals.csv', sep=",",stringsAsFactors=F,header=T,row.names = 1)
scaled_exp_tab <-  read.table('../data/results/scaled.csv', sep=",",stringsAsFactors=F,header=T,row.names = 1)


sample_df <- read.table('../data/results/sample_info.csv', sep=",",stringsAsFactors=F,header=T,row.names = 1)

sobj <- CreateSeuratObject(counts = exp_tab)
# Normalized data
sobj[["RNA"]]@data <- data.matrix(norm_exp_tab)
# Scaled data
sobj[["RNA"]]@scale.data <- data.matrix(scaled_exp_tab)

clustering.results <- data.frame(seurat_clusters = sample_df$top_cluster,row.names = colnames(x = sobj))
sobj <- AddMetaData(object = sobj, metadata = clustering.results)
Idents(object = sobj) <- colnames(x = clustering.results)[ncol(x = clustering.results)]
sobj[['seurat_clusters']] <- Idents(object = sobj)
n_clusters = length(unique(sample_df$top_cluster))
clusters = unique(sample_df$top_cluster)

# Find Cluster BioMarkers
marker.tables <- vector(mode = "list", length = n_clusters)
for (i in c(0:(n_clusters-1))) {
	markers <- FindMarkers(sobj, i, clusters[clusters!=i], only.pos = TRUE,  slot="scale.data", min.pct = 0.25, logfc.threshold = 0.25)
        df <- as.data.frame(markers)
	df$Gene <- rownames(df)
	df$cluster <- rep(i,dim(markers)[1])
        marker.tables[[i+1]] <- df
}
df <- do.call(rbind, marker.tables)


# Top 5 markers per cluster based on highest logFC
#sobj.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% View()
marker.table <- df %>% group_by(cluster) %>% top_n(n = 5, wt = avg_diff)
write.table(marker.table,'../data/results/markers.csv', sep=',',row.names = FALSE, quote=FALSE)


# DE between samples
sample_df[sample_df$s==1,'top_cluster'] <- sample_df[sample_df$s==1,'top_cluster']+n_clusters

clustering.results <- data.frame(seurat_clusters = sample_df$top_cluster,row.names = colnames(x = sobj))
sobj <- AddMetaData(object = sobj, metadata = clustering.results)
Idents(object = sobj) <- colnames(x = clustering.results)[ncol(x = clustering.results)]
sobj[['seurat_clusters']] <- Idents(object = sobj)

marker.tables <- vector(mode = "list", length = n_clusters)
for (i in c(0:(n_clusters-1))) {
	markers <- FindMarkers(sobj,i,i+n_clusters, slot="scale.data",  min.pct = 0.25, logfc.threshold = 0.25)
	df <- as.data.frame(markers)
	df$Gene <- rownames(df)
	df$cluster <- rep(i,dim(markers)[1])
	marker.tables[[i+1]] <- df
}

df <- do.call(rbind, marker.tables)
write.table(df,'../data/results/markers_samples.csv', sep=',', quote=FALSE)
