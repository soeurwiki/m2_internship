library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
sample <- args[2]

 if(args[1] == 'first'){
    print("first")
    name <- paste0(sample, "_NORM_0.html")
    dir = str_split_1(sample, '_')[1]

    rmarkdown::render("one_sample_unmerged.Rmd",
        output_file = name,
        params = list("sample" = sample),
        output_dir = paste0("./results/html/", dir))

} else if(args[1] == "samp"){
    print("samp")
    name <- paste0(sample, "_NORM_1.html")
    dir = str_split_1(sample, '_')[1]

    rmarkdown::render("one_sample_NORM_1.Rmd",
        output_file = name,
        params = list("sample" = sample),
        output_dir = paste0("./results/html/", dir))

} else if(args[1] == "line"){
    print("line")
    name <- paste0("line_", sample, "_NORM_1.html")

    rmarkdown::render("one_line_NORM_1.Rmd",
        output_file = name,
        params = list("sample" = sample),
        output_dir = paste0("./results/html/", sample))

}