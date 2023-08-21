# Seurat object  
## Load required packages  

library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library("SCINA")


## Load seurat object (RAW)  

pdf("norm4.pdf")
seurat_obj <- readRDS(file = "merged_seurat_Mca.rds")

as.data.frame(seurat_obj@assays$RNA[,]) -> scina.data

mature_oligo <- c("SOX8", "SOX11", "CLDN11","MBP","SOX10","SOX4","MOG","MYT1", "CNP", "PLP1", "OLIG1", "OLIG2", "NKX2-2", "ERBB3", "UGT8", "SOX17", "GPR17", "TNR")
oligo_progenitors  <- c("PDGFRA", "CSPG4")
astro <- c("SLC1A3", "NFIA","SOX9","GFAP","APOE", "AQP4","ALDH1L1", "FABP7", "TNC")

markers_astro_oligo <- list("astro"=astro, "mature_oligo"=mature_oligo, "oligo_progenitors"=oligo_progenitors)

SCINA(
  scina.data,
  markers_astro_oligo, 
  max_iter = 100, 
  convergence_n = 10, 
  convergence_rate = 0.999, 
  sensitivity_cutoff = 0.9, 
  rm_overlap=TRUE, 
  allow_unknown=TRUE
) -> scina.results

seurat_obj$scina_labels <- scina.results$cell_labels

colors <- c(astro='seagreen4', oligo_progenitors = "#297fb8" ,mature_oligo='orangered3', unknown="lightgray")

DimPlot(seurat_obj,reduction = "umap", pt.size = 1, label = TRUE, group.by = "scina_labels", label.size = 5, cols = colors)
DimPlot(seurat_obj,reduction = "tsne", pt.size = 1, label = TRUE, group.by = "scina_labels", label.size = 5, cols = colors)


#panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

# restricting the analysis to brain specific gene signatues
#panglao_brain <- panglao %>% filter(organ == "Brain")
  ##             panglao %>%  filter(str_detect(species,"Hs"))
  ## To obtain gene signatures for all genes

# restricting to human specific genes
#panglao_brain <- panglao_brain %>%  filter(str_detect(species,"Hs"))

# converting dataframes into a list of vectors, which is the format needed as input for CellID
#panglao_brain <- panglao_brain %>%  
#  group_by(`cell type`) %>%  
#  summarise(geneset = list(`official gene symbol`))

#brain_signatures <- setNames(panglao_brain$geneset, panglao_brain$`cell type`)

enrichR_db = "KEGG_2021_Human"
nb_maxgenes = 200

lapply(
  levels(seurat_obj$seurat_clusters),
  function(x)DEenrichRPlot(seurat_obj,  ident.1 = x,  max.genes = nb_maxgenes, return.gene.list = FALSE, enrich.database =  enrichR_db)
)

dev.off()
print("done")