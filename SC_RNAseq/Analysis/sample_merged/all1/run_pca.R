#   name : run_pca.R
#
#   Author (2023)  Safiya ATIA

args <- commandArgs(trailingOnly=TRUE)

if( args == "all"){
        print("all")
        rmarkdown::render("all_samples_NORM_1.Rmd",
        output_dir = "./results/html/")

} else if( args == "metabo"){
        print("metabo")
        rmarkdown::render("all_samples_metabo_1.Rmd",
                output_dir = "./results/html/")
}
