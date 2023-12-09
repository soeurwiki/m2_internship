#   name : run_3.R
#
#   Author (2023)  Safiya ATIA

args <- commandArgs(trailingOnly=TRUE)

line <- args[2]

if(args[1] == "all_mark"){
    print("all mark")
    name <- paste0(line, "_NORM_3.html")
    rmarkdown::render("one_line_NORM_3.Rmd",
        output_file = name,
        params = list("line" = line, "res" = args[3]),
        output_dir = paste0("./results/html/", line))

} else if(args[1] == "all_annot"){
    print("all annot")
    name <- paste0(line, "_NORM_4.html")
    rmarkdown::render("one_line_NORM_4.Rmd",
        output_file = name,
        params = list("line" = line),
        output_dir = paste0("./results/html/", line))

} else if(args[1] == "all_mca"){
    print("all mca")
    library(CellID)
    library(readr)
    library(dplyr)

    nb_maxgenes = 200
    pc = 50

    seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Markers.rds'))
    seurat_obj <- RunMCA(seurat_obj)
    saveRDS(seurat_obj, file = paste0('./results/rds/samples/', line, '_Mca.rds'))

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
    saveRDS(HGT_brain_signatures, file = paste0('./results/rds/samples/', line, '_HGT_brain_signatures.rds'))

} else if(args[1] == "metabo_mark"){
    print("metabo mark")
    name <- paste0(line, "_metabo_3.html")
    rmarkdown::render("one_line_metabo_3.Rmd",
        output_file = name,
        params = list("line" = line, "res" = args[3]),
        output_dir = paste0("./results/html/", line))

} else if(args[1] == "metabo_annot"){
    print("metabo annot")
    name <- paste0(line, "_metabo_4.html")
    rmarkdown::render("one_line_metabo_4.Rmd",
        output_file = name,
        params = list("line" = line),
        output_dir = paste0("./results/html/", line))

} else if(args[1] == "metabo_mca"){
    print("metabo mca")
    library(CellID)
    library(readr)
    library(dplyr)

    nb_maxgenes = 200
    pc = 50

    seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Markers_metabo.rds'))
    seurat_obj <- RunMCA(seurat_obj)
    saveRDS(seurat_obj, file = paste0('./results/rds/samples/', line, '_Mca_metabo.rds'))

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
    saveRDS(HGT_brain_signatures, file = paste0('./results/rds/samples/', line, '_HGT_brain_signatures_metabo.rds'))
}
