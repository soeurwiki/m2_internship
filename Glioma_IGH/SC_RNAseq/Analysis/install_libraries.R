if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("glmGamPoi")
install.packages("sctransform")
BiocManager::install("CelliD")


## convert Rmd to ipynb
devtools::install_github("mkearney/rmd2jupyter")
library(rmd2jupyter)
#rmd2jupyter(x = "Scripts/All_genes/one_sample_RAW.Rmd")