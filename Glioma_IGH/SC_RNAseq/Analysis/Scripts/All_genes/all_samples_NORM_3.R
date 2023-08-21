# Seurat object  
## Load required packages  

library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(clustree)


## Load seurat object (RAW)  

seurat_obj <- readRDS(file = "merged_seurat_FindNeighbors.rds")

res = 0.3
pdf("norm_3.pdf")
seurat_obj <- FindClusters(seurat_obj, resolution = res)
levels(seurat_obj$seurat_clusters)

DimPlot(seurat_obj, reduction = "pca", label = TRUE, group.by = "seurat_clusters") + ggtitle("PC1 vs PC2 with Clusters")

DimPlot(seurat_obj, reduction = "umap", pt.size = 1, label.size = 7, label = TRUE, group.by = "seurat_clusters")

DimPlot(seurat_obj, reduction = "tsne", pt.size = 1, label.size = 7, label = TRUE, group.by = "seurat_clusters")

VlnPlot(seurat_obj, features="nCount_RNA", group.by = "seurat_clusters", pt.size = 0)
RidgePlot(seurat_obj, features = "nCount_RNA", group.by = "seurat_clusters")

VlnPlot(seurat_obj, features="nFeature_RNA", group.by = "seurat_clusters", pt.size = 0)
RidgePlot(seurat_obj, features = "nFeature_RNA", group.by = "seurat_clusters")

VlnPlot(seurat_obj, features = "percent.mt", group.by = "seurat_clusters", pt.size = 0)
RidgePlot(seurat_obj, features = "percent.mt", group.by = "seurat_clusters")

VlnPlot(seurat_obj, features = "percent.ribosomal", group.by = "seurat_clusters", pt.size = 0)
RidgePlot(seurat_obj, features = "percent.ribosomal", group.by = "seurat_clusters")

VlnPlot(seurat_obj, features = "percent.largest_gene", group.by = "seurat_clusters", pt.size = 0)
RidgePlot(seurat_obj, features = "percent.largest_gene", group.by = "seurat_clusters")

seurat_obj[[]] %>%
  group_by(seurat_clusters, largest_gene) %>%
  count() %>%
  arrange(desc(n)) %>%
  group_by(seurat_clusters) %>%
  slice(1:2) %>%
  ungroup() %>%
  arrange(seurat_clusters, desc(n))

seurat_obj@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")

print("cluster markers")

lapply(
  levels(seurat_obj$seurat_clusters),
  function(x)FindConservedMarkers(seurat_obj, ident.1 = x, grouping.var='condition')
) -> cluster.markers

# This adds the cluster number to the results of FindConservedMarkers
sapply(0:(length(cluster.markers)-1),function(x) {
  cluster.markers[[x+1]]$gene <<- rownames(cluster.markers[[x+1]])
  cluster.markers[[x+1]]$cluster <<- x
})

print("population markers")
lapply(
  levels(seurat_obj$seurat_clusters),
  function(x)FindConservedMarkers(seurat_obj, ident.1 = x, grouping.var='orig.ident')
) -> population.markers

# This adds the cluster number to the results of FindConservedMarkers
sapply(0:(length(population.markers)-1),function(x) {
  population.markers[[x+1]]$gene <<- rownames(population.markers[[x+1]])
  population.markers[[x+1]]$cluster <<- x
})

dev.off()
write.table(cluster.markers, "cluster_markers.txt", append = FALSE, sep = " ")
write.table(population.markers, "population_markers.txt", append = FALSE, sep = " ")


saveRDS(seurat_obj, file = "merged_seurat_Markers.rds")

print("done")