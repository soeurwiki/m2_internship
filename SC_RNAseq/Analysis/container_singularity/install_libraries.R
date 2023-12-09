#   name : install_libraries.R
#
#   Author (2023)  Safiya ATIA

install.packages("BiocManager")
install.packages("devtools")
setRepositories(ind = c(1,2,3))
install.packages("remotes")

install.packages("tidyverse")

## ---------------- Rmarkdown --------------- ##
install.packages("rmarkdown")
install.packages("installr")
installr::install.pandoc()

## --------- Visualisations & plots --------- ##
install.packages("ggplot2")
install.packages("viridis")
install.packages("scales")
install.packages("ggpubr")

## ------------ Seurat functions ------------ ##
install.packages('Seurat')
BiocManager::install("glmGamPoi")
BiocManager::install("harmony")
install.packages("sctransform")
install.packages("clustree")
install.packages("future")

## DoubletFinder
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

## FindConservative Markers
BiocManager::install('multtest') 
install.packages('metap')

## --------- Annotations & Pathways --------- ##
BiocManager::install("scater")
BiocManager::install("destiny")
BiocManager::install("SingleR")

BiocManager::install('SCINA') 
BiocManager::install("SiPSiC")
install.packages('enrichR')