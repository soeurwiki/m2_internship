args <- commandArgs(trailingOnly=TRUE)
name <- paste0(args[1],"_RAW_filtered.html")
rmarkdown::render("one_sample_RAW_2.Rmd",
        output_file = name,
        params = list("sample" = args[1], "pc" = args[2], "doublet_rate"=args[3]),
        output_dir = "./results/html/")