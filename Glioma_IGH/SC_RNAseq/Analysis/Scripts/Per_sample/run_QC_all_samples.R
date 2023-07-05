sample.list <- commandArgs(trailingOnly=TRUE)

render_one <- function(sample) {
  rmarkdown::render(
    'one_sample_RAW.Rmd',
    output_file = paste0(sample,'_RAW.html'),
    params = list("sample" = sample),
    envir = parent.frame()
  )
}

for (sample in sample.list) {
    render_one(sample)
}