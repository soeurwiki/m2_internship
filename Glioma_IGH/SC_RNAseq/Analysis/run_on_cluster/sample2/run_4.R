args <- commandArgs(trailingOnly=TRUE)

if(args[1] == "mark"){
    print("mark")
    rmarkdown::render("all_samples_NORM_3.Rmd",
        output_dir = "./results/html/")
} else if(args[1] == "annot"){
    print("annot")
    rmarkdown::render("all_samples_NORM_4.Rmd",
        output_dir = "./results/html/")
} else if(args[1] == "mca"){
    print("mca")
    sample <- args[2]

    library(CellID)
    library(readr)
    library(dplyr)

    seurat_obj <- readRDS(file = paste0('./results/rds/', sample, 'merged_seurat_Markers.rds'))
    seurat_obj <- RunMCA(seurat_obj)
    saveRDS(seurat_obj, file = paste0('./results/rds/', sample, '_Mca.rds'))

    pdf("./results/PlotMC.pdf")
    DimPlotMC(seurat_obj, reduction = "mca", group.by = 'orig.ident', features = c("APOE"), as.text = T)
    dev.off()

    # download all cell-type gene signatures from panglaoDB
    panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

    # restricting the analysis to brain specific gene signatues
    panglao_brain <- panglao %>% filter(organ == "Brain")

    # restricting to human specific genes
    panglao_brain <- panglao_brain %>%  filter(grepl("Hs", species))

    # converting dataframes into a list of vectors, which is the format needed as input for CellID
    panglao_brain <- panglao_brain %>%  
                            group_by(`cell type`) %>%  
                            summarise(geneset = list(`official gene symbol`))

    brain_signatures <- setNames(panglao_brain$geneset, panglao_brain$`cell type`)

    # Performing per-cell hypergeometric tests against the gene signature collection
    HGT_brain_signatures <- RunCellHGT(seurat_obj, pathways = brain_signatures, dims = 1:50, n.features = 200)
    saveRDS(HGT_brain_signatures, file = paste0('./results/rds/', sample, '_HGT_brain_signatures.rds'))
}
