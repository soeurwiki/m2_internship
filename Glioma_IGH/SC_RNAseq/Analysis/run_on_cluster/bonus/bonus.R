#rmarkdown::render("clustree.Rmd",
#   output_dir = "./results/html/")

sample <- commandArgs(trailingOnly=TRUE)
name <- paste0(sample, "_plots_QC.html")
rmarkdown::render("plots_lines.Rmd",
    output_file = name,
    params = list("sample" = sample),
    output_dir = "./results/html/samples/")