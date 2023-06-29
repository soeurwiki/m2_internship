library(Seurat)

seurat_obj <- readRDS(file = "merged_seurat_Clusters_1.rds")

seurat_obj <- RunMCA(seurat_obj)

saveRDS(seurat_obj, file ="merged_seurat_Mca.rds")
