library(Seurat)

seurat_obj <- readRDS(file = './results/rds/merged_seurat_Clusters.rds')

## Create seurat object (ALL GENES) for each samples
BT1 <- subset(x = seurat_obj, subset = orig.ident %in% c("BT1_prolif", "BT1_diff"))
saveRDS(BT1_prolif, file ='./results/rds/samples/norm_BT1.rds')

BT2 <- subset(x = seurat_obj, subset = orig.ident %in% c("BT2_prolif","BT2_diff") )
saveRDS(BT2_prolif, file ='./results/rds/samples/norm_BT2.rds')

BT54 <- subset(x = seurat_obj, subset = orig.ident %in% c("BT54_prolif", "BT54_diff"))
saveRDS(BT54_prolif, file ='./results/rds/samples/norm_BT54.rds')

BT88 <- subset(x = seurat_obj, subset = orig.ident %in% c("BT88_prolif", "BT88_diff"))
saveRDS(BT88_prolif, file ='./results/rds/samples/norm_BT88.rds')

LGG85 <- subset(x = seurat_obj, subset = orig.ident %in% c("LGG85_prolif","LGG85_diff"))
saveRDS(LGG85_prolif, file ='./results/rds/samples/norm_LGG85.rds')

LGG275 <- subset(x = seurat_obj, subset = orig.ident %in% c("LGG275_prolif","LGG275_diff"))
saveRDS(LGG275_prolif, file ='./results/rds/samples/norm_LGG275.rds')

LGG336 <- subset(x = seurat_obj, subset = orig.ident %in% c("LGG336_prolif", "LGG336_diff"))
saveRDS(LGG336_prolif, file ='./results/rds/samples/norm_LGG336.rds')

LGG349 <- subset(x = seurat_obj, subset = orig.ident %in% c("LGG349_prolif","LGG349_diff"))
saveRDS(LGG349_prolif, file ='./results/rds/samples/norm_LGG349.rds')