library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)

library(ggplot2)
library("scales")
library(clustree)

metabo_obj <- readRDS(file = "./merged_seurat_metabo.rds")
metabo_obj
table(metabo_obj$orig.ident)

pc = 50
8482 -> saved.seed
set.seed(saved.seed)

pdf("metabo_2.pdf")

metabo_obj <- RunUMAP(metabo_obj, dims = 1:pc, seed.use = saved.seed, )

DimPlot(metabo_obj, reduction = "umap", group.by = "orig.ident")
DimPlot(metabo_obj, reduction = "umap", group.by = "Phase", cols = c(G1='azure4',G2M='dodgerblue', S='tomato3') )
metabo_obj <- RunTSNE(metabo_obj, dims=1:pc,seed.use = saved.seed)

DimPlot(metabo_obj, reduction = "tsne", group.by = "orig.ident")
DimPlot(metabo_obj, reduction = "tsne", group.by = "Phase", cols = c(G1='azure4',G2M='dodgerblue', S='tomato3'))

metabo_obj <- FindNeighbors(metabo_obj, dims = 1:pc)
saveRDS(metabo_obj, file ="merged_seurat_FindNeighbors_metabo.rds")


metabo_obj <- FindClusters(metabo_obj, resolution = 0.2)


metabo_obj <- FindClusters(metabo_obj, resolution = 0.3)


metabo_obj <- FindClusters(metabo_obj, resolution = 0.4)


metabo_obj <- FindClusters(metabo_obj, resolution = 0.5)


metabo_obj <- FindClusters(metabo_obj, resolution = 0.6)


metabo_obj <- FindClusters(metabo_obj, resolution = 0.7)


metabo_obj <- FindClusters(metabo_obj, resolution = 0.8)


metabo_obj <- FindClusters(metabo_obj, resolution = 0.9)


metabo_obj <- FindClusters(metabo_obj, resolution = 1)


metabo_obj <- FindClusters(metabo_obj, resolution = 1.1)


metabo_obj <- FindClusters(metabo_obj, resolution = 1.2)


metabo_obj <- FindClusters(metabo_obj, resolution = 1.3)


metabo_obj <- FindClusters(metabo_obj, resolution = 1.4)

clustree(metabo_obj, prefix = "RNA_snn_res.")
dev.off()

print("done")