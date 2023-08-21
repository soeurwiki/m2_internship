args <- commandArgs(trailingOnly=TRUE)

line <- args[2]

if( args[1] == "all"){
        print("all")
        name <- paste0(line, "_NORM_1.html")

        rmarkdown::render("one_line_NORM_1.Rmd",
                output_file = name,
                params = list("line" = line),
                output_dir = paste0("./results/html/", line))

} else if( args[1] == "metabo"){
        print("metabo")
        name <- paste0(line, "_metabo_1.html")

        rmarkdown::render("one_line_metabo_1.Rmd",
        output_file = name,
        params = list("line" = line),
        output_dir = paste0("./results/html/", line))
}
