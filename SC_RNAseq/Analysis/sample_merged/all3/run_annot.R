#   name : run_annot.R
#
#   Author (2023)  Safiya ATIA

args <- commandArgs(trailingOnly=TRUE)

if(args[1] == "all_mark"){
    print("all mark")
    rmarkdown::render("all_samples_NORM_3.Rmd",
        output_dir = "./results/html/")

} else if(args[1] == "all_annot"){
    print("all annot")
    rmarkdown::render("all_samples_NORM_4.Rmd",
        output_dir = "./results/html/")

} else if(args[1] == "all_mca"){
    print("all mca")
    library(CellID)
    library(readr)
    library(dplyr)

    nb_maxgenes = 200
    pc = 40

    seurat_obj <- readRDS(file = './results/rds/all_samples_Markers.rds')
    seurat_obj <- RunMCA(seurat_obj)
    saveRDS(seurat_obj, file = './results/rds/all_samples_Mca.rds')

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
    HGT_brain_signatures <- RunCellHGT(seurat_obj, pathways = brain_signatures, dims = 1:pc, n.features = nb_maxgenes)
    saveRDS(HGT_brain_signatures, file = './results/rds/all_samples_HGT_brain_signatures.rds')
} else if(args == "metabo_mark"){
    print("metabo mark")
    rmarkdown::render("all_samples_metabo_3.Rmd",
        output_dir = "./results/html/")

} else if(args == "metabo_annot"){
    print("metabo annot")
    rmarkdown::render("all_samples_metabo_4.Rmd",
        output_dir = "./results/html/")

} else if(args == "metabo_mca"){
    print("metabo mca")
    library(CellID)
    library(readr)
    library(dplyr)

    nb_maxgenes = 200
    pc = 30

    seurat_obj <- readRDS(file = './results/rds/all_samples_Markers_metabo.rds')
    seurat_obj <- RunMCA(seurat_obj)
    saveRDS(seurat_obj, file = './results/rds/all_samples_Mca_metabo.rds')

    pdf("./results/PlotMC_metabo.pdf")
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
    HGT_brain_signatures <- RunCellHGT(seurat_obj, pathways = brain_signatures, dims = 1:pc, n.features = nb_maxgenes)
    saveRDS(HGT_brain_signatures, file = './results/rds/all_samples_HGT_brain_signatures_metabo.rds')
}

