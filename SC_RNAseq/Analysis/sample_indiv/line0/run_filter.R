#   name : run_filter.R
#
#   Author (2023)  Safiya ATIA

args <- commandArgs(trailingOnly=TRUE)
sample <- args[2]

if(args[1] == 'qc'){
        name <- paste0(sample,"_RAW_qc.html")
        rmarkdown::render("one_sample_RAW_1.Rmd",
                output_file = name,
                params = list("sample" = sample),
                output_dir = "./results/html/samples/")

} else if(args[1] == 'filter'){
        name <- paste0(sample, "_RAW_filtered.html")
        rmarkdown::render("one_sample_RAW_2.Rmd",
                output_file = name,
                params = list("sample" = sample, "pc" = args[3], "doublet_rate"=args[4]),
                output_dir = "./results/html/samples/")
}