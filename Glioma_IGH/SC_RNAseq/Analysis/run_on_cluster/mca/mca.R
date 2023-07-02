library(Seurat)

seurat_obj <- readRDS(file = "merged_seurat_Markers_metabo.rds")
seurat_obj <- RunMCA(seurat_obj)
saveRDS(seurat_obj, file ="merged_seurat_Mca_metabo.rds")


seurat_obj <- readRDS(file = "merged_seurat_Markers.rds")
seurat_obj <- RunMCA(seurat_obj)
saveRDS(seurat_obj, file ="merged_seurat_Mca.rds")
