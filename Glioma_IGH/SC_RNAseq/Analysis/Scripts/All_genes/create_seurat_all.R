library(Seurat)

path = "./results/rds/"

seurat.BT1_prolif <- readRDS(file = paste0(path, "Filtered_BT1_prolif.rds"))
seurat.BT1_diff <- readRDS(file = paste0(path,"Filtered_BT1_diff.rds"))
seurat.BT2_prolif <- readRDS(file = paste0(path,"Filtered_BT2_prolif.rds"))
seurat.BT2_diff <- readRDS(file = paste0(path,"Filtered_BT2_diff.rds"))
seurat.BT54_prolif <- readRDS(file = paste0(path,"Filtered_BT54_prolif.rds"))
seurat.BT54_diff <- readRDS(file = paste0(path,"Filtered_BT54_diff.rds"))
seurat.BT88_prolif <- readRDS(file = paste0(path,"Filtered_BT88_prolif.rds"))
seurat.BT88_diff <- readRDS(file = paste0(path,"Filtered_BT88_diff.rds"))
seurat.LGG85_prolif <- readRDS(file = paste0(path,"Filtered_LGG85_prolif.rds"))
seurat.LGG85_diff <- readRDS(file = paste0(path,"Filtered_LGG85_diff.rds"))
seurat.LGG275_prolif <- readRDS(file = paste0(path,"Filtered_LGG275_prolif.rds"))
seurat.LGG275_diff <- readRDS(file = paste0(path,"Filtered_LGG275_diff.rds"))
seurat.LGG336_prolif <- readRDS(file = paste0(path,"Filtered_LGG336_prolif.rds"))
seurat.LGG336_diff <- readRDS(file = paste0(path,"Filtered_LGG336_diff.rds"))
seurat.LGG349_prolif <- readRDS(file = paste0(path,"Filtered_LGG349_prolif.rds"))
seurat.LGG349_diff <- readRDS(file = paste0(path,"Filtered_LGG349_diff.rds"))

## Merge all seurat objects  
seurat_obj <- merge(seurat.BT1_prolif, y = c(seurat.BT1_diff, seurat.BT2_prolif, seurat.BT2_diff ,seurat.BT54_prolif, seurat.BT54_diff, seurat.BT88_prolif, seurat.BT88_diff, seurat.LGG85_prolif, seurat.LGG85_diff, seurat.LGG275_prolif, seurat.LGG275_diff, seurat.LGG336_prolif, seurat.LGG336_diff, seurat.LGG349_prolif, seurat.LGG349_diff), add.cell.ids = c("BT1_prolif","BT1_diff", "BT2_prolif", "BT2_diff", "BT54_prolif", "BT54_diff", "BT88_prolif", "BT88_diff", "LGG85_prolif", "LGG85_diff", "LGG275_prolif", "LGG275_diff", "LGG336_prolif", "LGG336_diff", "LGG349_prolif", "LGG349_diff"))

saveRDS(seurat_obj, file = "./results/mtx_merged.rds")