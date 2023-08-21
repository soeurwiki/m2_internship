sample.list <- commandArgs(trailingOnly=TRUE)

render_one <- function(sample) {
  rmarkdown::render(
    'one_sample_RAW_1.Rmd',
    output_file = paste0(sample,'_RAW.html'),
    params = list("sample" = sample),
    output_dir = "./results/html/"
  )
}

for (sample in sample.list) {
    print(paste0(" QC for ", sample, " ..."))
    render_one(sample)
    print(paste0(" seurat object QC_", sample, ".rds created."))
    print(" Rmarkdown created.")
}