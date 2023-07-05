# module load system/R-4.2.3_Miniconda3

Rscript create_seurat_samples.R

sample=$(ls ./results | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

Rscript run_QC_all_samples.R ${sample}