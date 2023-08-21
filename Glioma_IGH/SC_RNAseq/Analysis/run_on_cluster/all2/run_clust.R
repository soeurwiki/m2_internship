args <- commandArgs(trailingOnly=TRUE)

if(args == "all"){
        print("all")
        rmarkdown::render("all_samples_NORM_2.Rmd",
                output_dir = "./results/html/")
} else if(args == "metabo"){
        print("metabo")
        rmarkdown::render("all_samples_metabo_2.Rmd",
                output_dir = "./results/html/")
}