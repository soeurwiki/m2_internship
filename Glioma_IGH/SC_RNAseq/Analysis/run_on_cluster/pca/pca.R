library(Seurat)

args <- commandArgs(trailingOnly=TRUE)
pc = 40

if (args[1] == "metabo") {
    seurat_obj <- readRDS(file = "merged_seurat_metabo.rds")
    pdf_jack = "JackStrawPlot_metabo.pdf"
    pdf_elbow = "ElbowPlot_metabo.pdf"
    pdf_heatmap = "DimHeatmap_metabo.pdf"
} else{
    seurat_obj <- readRDS(file = "Norm_merged_seurat.rds")
    pdf_jack = "JackStrawPlot.pdf"
    pdf_elbow = "ElbowPlot.pdf"
    pdf_heatmap = "DimHeatmap.pdf"
}

seurat_obj <- JackStraw(seurat_obj, num.replicate = 100, dims = pc)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:pc)

pdf(pdf_jack)
JackStrawPlot(seurat_obj, dims = 1:pc)
dev.off()

pdf(pdf_elbow)
ElbowPlot(seurat_obj, ndims = pc)
dev.off()

pdf(pdf_heatmap)
DimHeatmap(seurat_obj, dims = 1:20, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 21:pc, balanced = TRUE)
dev.off()