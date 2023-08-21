args <- commandArgs(trailingOnly=TRUE)

line <- args[2]

if( args[1] == "all"){
        print("all")
        name <- paste0(line, "_NORM_2.html")

        rmarkdown::render("one_line_NORM_2.Rmd",
                output_file = name,
                params = list("line" = line),
                output_dir = paste0("./results/html/", line))
} else if( args[1] == "metabo"){
        print("metabo")
        name <- paste0(line, "_metabo_2.html")

        rmarkdown::render("one_line_metabo_2.Rmd",
                output_file = name,
                params = list("line" = line),
                output_dir = paste0("./results/html/", line))
}
