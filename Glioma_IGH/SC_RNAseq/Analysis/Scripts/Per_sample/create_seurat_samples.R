library(Seurat)

samples <- c("BT1", "BT2", "BT54", "BT88", "LGG85", "LGG275", "LGG336", "LGG349")

lapply(
  samples,
  function(x){
    ### sample GF
    condition = "prolif"
    path <- paste0("./samples/", x, "_", condition, "/")
    project <- paste0(x, "_", condition)

    mtx_obj <- Seurat::ReadMtx(mtx = paste0(path, "matrix.mtx.gz"),
                            cells = paste0(path, "barcodes.tsv.gz"),
                            features = paste0(path, "features.tsv.gz"))

    seurat_obj <- CreateSeuratObject(counts = mtx_obj,
                                            min.cells = 3, min.features = 200,
                                            project = project)
    seurat_obj@meta.data$condition <- condition

    saveRDS(seurat_obj, file = paste0("./results/mtx_", x, "_", condition, ".rds"))

    ### sample noGF
    condition = "diff"
    path <- paste0("./samples/", x, "_", condition, "/")
    project <- paste0(x, "_", condition)

    mtx_obj <- Seurat::ReadMtx(mtx = paste0(path,"matrix.mtx.gz"),
                            cells = paste0(path,"barcodes.tsv.gz"),
                            features = paste0(path,"features.tsv.gz"))

    seurat_obj <- CreateSeuratObject(counts = mtx_obj,
                                            min.cells = 3, min.features = 200,
                                            project = project)
    seurat_obj@meta.data$condition <- condition

    saveRDS(seurat_obj, file = paste0("./results/mtx_", x, "_", condition, ".rds"))
  }
)