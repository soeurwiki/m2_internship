# Seurat object  
## Load required packages  

library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(clustree)


## Load seurat object (RAW)  

seurat_obj <- readRDS(file = "Norm_merged_seurat.rds")

print(seurat_obj)
table(seurat_obj$orig.ident)

pc = 50
8482 -> saved.seed
set.seed(saved.seed)

pdf("norm_2.pdf") 

seurat_obj <- RunUMAP(seurat_obj, dims = 1:pc, seed.use = saved.seed, )

DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")
DimPlot(seurat_obj, reduction = "umap", group.by = "Phase", cols = c(G1='azure4',G2M='dodgerblue', S='tomato3') )

seurat_obj <- RunTSNE(seurat_obj, dims=1:pc,seed.use = saved.seed)

DimPlot(seurat_obj, reduction = "tsne", group.by = "orig.ident")
DimPlot(seurat_obj, reduction = "tsne", group.by = "Phase", cols = c(G1='azure4',G2M='dodgerblue', S='tomato3'))


seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pc)

saveRDS(seurat_obj, file ="merged_seurat_FindNeighbors.rds")

seurat_obj <- FindClusters(seurat_obj, resolution = 0.2)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.7)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.9)

seurat_obj <- FindClusters(seurat_obj, resolution = 1)

seurat_obj <- FindClusters(seurat_obj, resolution = 1.1)

seurat_obj <- FindClusters(seurat_obj, resolution = 1.2)

seurat_obj <- FindClusters(seurat_obj, resolution = 1.3)

seurat_obj <- FindClusters(seurat_obj, resolution = 1.4)

clustree(seurat_obj, prefix = "RNA_snn_res.")
dev.off()