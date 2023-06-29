library(Seurat)

seurat_obj <- readRDS(file = "Norm_merged_seurat.rds")

seurat_obj <- JackStraw(seurat_obj, num.replicate = 100, dims = 40)
jack.data <- seurat_obj

seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:40)

JackStrawPlot(seurat_obj, dims = 1:40)


pdf("JackStrawPlot.pdf") 
JackStrawPlot(seurat_obj, dims = 1:40)
dev.off()

pdf("ElbowPlot.pdf") 
ElbowPlot(seurat_obj, ndims = 40)
dev.off()