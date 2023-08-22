# Scripts order

*Scripts are written to be run on the cluster genotoul*

## 1- Processing

-   Run `get_files.sh`
-   Run `run_pipeline.sh`

``` bash
source ./SC_RNAseq/Processing/get_files.sh
source ./SC_RNAseq/Processing/run_pipeline.sh
```

## 2- Analysis

### a) Construction of singularity image for R
``` bash
sudo singularity build --writable-tmpfs r.sif ./SC_RNAseq/Analysis/container_singularity/singularity-r.def
```

### b) Filtering samples
``` bash
sbatch --wrap=" ./SC_RNAseq/Analysis/sample_indiv/line0/line_0.0.sh"
```
-   `line_0.0.sh` 
    -   `create_seurat_samples.R`

``` bash
sbatch --wrap=" ./SC_RNAseq/Analysis/sample_indiv/line0/line_0.1.sh"
```
-   `line_0.1.sh` 
    -   `run_filter.R`
        -   `one_sample_RAW_1.Rmd`: Filter one sample based on nCount, nFeature_RNA and percent Mitochondrion.
--- Modify args.txt (pc + doublet_rate) ---

``` bash
sbatch --wrap=" ./SC_RNAseq/Analysis/sample_indiv/line0/line_0.2.sh"
```
-   `line_0.2.sh` 
    -   `run_filter.R`
        -   `one_sample_RAW_2.Rmd`

### c) Merging samples + analysis
*For an analysis line specific, replace `all1/all_1.sh` by `line1/line_1.sh`*

``` bash
sbatch --wrap=" ./SC_RNAseq/Analysis/sample_merged/all1/all_1.sh"
```
-   `all_1.sh` 
    -   `run_pca.R`
        -   `all_samples_NORM_1.Rmd`: run PCA
        -   `all_samples_metabo_1.Rmd`
--- Choose number of PC for all_genes & metabo ---

``` bash
sbatch --wrap=" ./SC_RNAseq/Analysis/sample_merged/all2/all_2.sh"
```
-   `all_2.sh` 
    -   `run_clust.R`
        -   `all_samples_NORM_2.Rmd`: run UMAP, tSNE and clustree (to visualize resolution for clusters)
        -   `all_samples_metabo_2.Rmd`
--- Choose resolution for all_genes & metabo ---

``` bash
sbatch --wrap=" ./SC_RNAseq/Analysis/sample_merged/all3/all_3.sh"
```
-   `all_3.sh` 
    -   `run_annot.R`
        -   `all_samples_NORM_3.Rmd`: Define clusters and Markers
        -   `all_samples_metabo_3.Rmd`
        -   `all_samples_NORM_4.Rmd`: Cell-type annotation and Ontology
        -   `all_samples_metabo_4.Rmd`


------------------------------------------------------------------------