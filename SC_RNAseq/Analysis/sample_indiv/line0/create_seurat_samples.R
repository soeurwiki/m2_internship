#   name : create_seurat_samples.R
#
#   Author (2023)  Safiya ATIA

library(Seurat)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

project <- args[1]
sample = str_split_1(project, '_')[1]
condition = str_split_1(project, '_')[2]

## !! Please verify sample names in cellranger directory, they must be sample_prolif / sample_diff !! ##


path <- paste0("../SC_cellranger/cellranger/count/sample-", project, "/outs/raw_feature_bc_matrix/")

mtx_obj <- Seurat::ReadMtx(mtx = paste0(path, "matrix.mtx.gz"),
cells = paste0(path, "barcodes.tsv.gz"),
features = paste0(path, "features.tsv.gz"))

seurat_obj <- CreateSeuratObject(counts = mtx_obj,
                                min.cells = 3, min.features = 200,
                                project = project)
seurat_obj@meta.data$condition <- condition

saveRDS(seurat_obj, file = paste0("./results/rds/mtx/mtx_", project, ".rds"))
print(paste0(" seurat object mtx_", project, ".rds created."))
