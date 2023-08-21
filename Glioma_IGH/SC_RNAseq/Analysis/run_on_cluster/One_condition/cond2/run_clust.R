args <- commandArgs(trailingOnly=TRUE)

samp <- args[2]

if( args[1] == "all"){
        print("all")
        name <- paste0(samp, "_NORM_2.html")

        rmarkdown::render("one_line_NORM_2.Rmd",
                output_file = name,
                params = list("samp" = samp),
                output_dir = paste0("./Lucille/results/html/", samp))
} else if( args[1] == "metabo"){
        print("metabo")
        name <- paste0(samp, "_metabo_2.html")

        rmarkdown::render("one_line_metabo_2.Rmd",
                output_file = name,
                params = list("samp" = samp),
                output_dir = paste0("./Lucille/results/html/", samp))
}
