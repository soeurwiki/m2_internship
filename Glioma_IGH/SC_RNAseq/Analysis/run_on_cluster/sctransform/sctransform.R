library(sctransform)

seurat_obj <- readRDS(file = "Raw_merged_seurat.rds")

seurat_obj <- SCTransform( seurat_obj, method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE)

saveRDS(seurat_obj, file = "sctransform.rds")
