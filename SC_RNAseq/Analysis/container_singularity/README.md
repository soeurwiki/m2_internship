--- SC_RNAseq ---

# container singularity

This folder contains :

-   singularity-r.def

-   install_libraries.R

-   test_libraries.R

-   test_libraries.Rmd

#### To build the image :

`sudo singularity build --writable-tmpfs r.sif singularity-r.def`

*It cannot be build on the genotoul cluster because you'll need admin rights.*

*You can run it on your personal computer and send `r.sif` on your genotoul space.*

#### Output :

-   r.sif

How to use it :

`singularity exec r.sif Rscript 'test_libraries.R`

`singularity exec r.sif Rscript -e  "rmarkdown::render('test_libraries.Rmd')"`

`singularity exec r.sif R`