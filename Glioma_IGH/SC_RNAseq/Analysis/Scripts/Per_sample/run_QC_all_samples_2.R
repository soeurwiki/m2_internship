BT1_diff <- list("sample"= "BT1_diff", "pc"=20, "doublet_rate"= 0.031) 
BT1_prolif <- list("sample"= "BT1_prolif", "pc"=20, "doublet_rate"= 0.016) 
BT2_diff <- list("sample"= "BT2_diff", "pc"=20, "doublet_rate"= 0.023) 
BT2_prolif <- list("sample"= "BT2_prolif", "pc"=20, "doublet_rate"= 0.008) 
BT54_diff <- list("sample"= "BT54_diff", "pc"=20, "doublet_rate"= 0.016) 
BT54_prolif <- list("sample"= "BT54_prolif", "pc"=15, "doublet_rate"= 0.054) 
BT88_diff <- list("sample"= "BT88_diff", "pc"=20, "doublet_rate"= 0.023) 
BT88_prolif <- list("sample"= "BT88_prolif", "pc"=20, "doublet_rate"= 0.023)
LGG85_diff <- list("sample"= "LGG85_diff", "pc"=15, "doublet_rate"= 0.016) 
LGG85_prolif <- list("sample"= "LGG85_prolif", "pc"=20, "doublet_rate"= 0.008)
LGG275_diff <- list("sample"= "LGG275_diff", "pc"=15, "doublet_rate"= 0.016) 
LGG275_prolif <- list("sample"= "LGG275_prolif", "pc"=20, "doublet_rate"= 0.016)
LGG336_diff <- list("sample"= "LGG336_diff", "pc"=20, "doublet_rate"= 0.039)
LGG336_prolif <- list("sample"= "LGG336_prolif", "pc"=20, "doublet_rate"= 0.023)
LGG349_diff <- list("sample"= "LGG349_diff", "pc"=20, "doublet_rate"= 0.016)
LGG349_prolif <- list("sample"= "LGG349_prolif", "pc"=20, "doublet_rate"= 0.016) 


sample.list <- list(BT1_prolif, BT1_diff,
                 BT2_prolif, BT2_diff,
                 BT54_prolif, BT54_diff,
                 BT88_prolif, BT88_diff,
                 LGG85_prolif, LGG85_diff,
                 LGG275_prolif, LGG275_diff,
                 LGG336_prolif, LGG336_diff,
                 LGG349_prolif, LGG349_diff)


render_one <- function(sample, pc, doublet_rate) {
  rmarkdown::render(
    'one_sample_RAW_2.Rmd',
    output_file = paste0(sample,'_RAW_filtered.html'),
    params = list("sample" = sample, "pc" = pc, "doublet_rate"=doublet_rate),
    output_dir = "./results/html/"
  )
}

for (sample in sample.list) {
    print(paste0(" QC for ", sample$sample, " ..."))
    render_one(sample= sample$sample, sample$pc, doublet_rate = sample$doublet_rate )
    print(paste0(" seurat object Filtered_", sample, ".rds created."))
    print(" Rmarkdown created.")
}