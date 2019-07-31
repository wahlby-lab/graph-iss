library(Seurat)
library(dplyr)

morpho_tab <-  read.table('../data/results/morpho_table.csv', sep=",",stringsAsFactors=F,header=T,row.names = 1)
norm_morpho_tab <- read.table('../data/results/normalized_morpho_table.csv', sep=",",stringsAsFactors=F,header=T,row.names = 1)
sample_df <-  read.table('../data/results/sample_info_morpho.csv', sep=",",stringsAsFactors=F,header=T,row.names = 1)

sobj <- CreateSeuratObject(counts = morpho_tab)

# Normalized data
sobj[["RNA"]]@data <- data.matrix(norm_morpho_tab)

all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes)

clustering.results <- data.frame(seurat_clusters = sample_df$top_cluster,row.names = colnames(x = sobj))
sobj <- AddMetaData(object = sobj, metadata = clustering.results)
Idents(object = sobj) <- colnames(x = clustering.results)[ncol(x = clustering.results)]
sobj[['seurat_clusters']] <- Idents(object = sobj)

# Find Cluster Morphological Markers
# Uses a logistic regression framework to determine differentially
# expressed genes. Constructs a logistic regression model predicting group
# membership based on each feature individually and compares this to a null
# model with a likelihood ratio test.
sobj.markers <- FindAllMarkers(sobj, test.use='LR', slot='data',  only.pos = TRUE)
marker.table <- sobj.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
write.table(marker.table,'../data/results/morpho_markers.csv', sep=',',row.names = FALSE, quote=FALSE)
