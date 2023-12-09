--- SC_RNAseq ---

# line0

This folder contains :

-   line_0.0.sh

-   line_0.1.sh

-   line_0.2.sh

-   one_sample_RAW_1.Rmd

-   one_sample_RAW_2.Rmd

-   run_filter.R

-   args.txt

#### args.txt :

*To fill the args.txt, use the table page 18:*

<https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf>

**/!\\ In order to be read properly, the args.txt must end with a** \n

(for 16 samples, there's 17 lines in args.txt)

#### Outputs :

-   mtx_project.rds
-   QC_project.rds
-   sample_RAW_qc.html
-   Filtered_project.rds
-   sample_RAW_filtered.html

(for examples : sample=BT1, project=BT1_diff, BT1_prolif )