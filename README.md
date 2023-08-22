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

**Construction of singularity image for R**
``` bash
sudo singularity build --writable-tmpfs r.sif ./SC_RNAseq/Analysis/container_singularity/singularity-r.def
```

### Filtering samples
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

--- Choisir le nombre de pc ---
``` bash
sbatch --wrap=" ./SC_RNAseq/Analysis/sample_indiv/line0/line_0.2.sh"
```

### a) Per sample

#### --- All genes ---

-   Run `.sh`

    ``` bash
    sbatch --wrap=" source ./.sh"
    ```

    -   run .R

        -   `one_sample_RAW_1.Rmd` : Default filtering : [nCount=500, nFeature=220, Mitochondrion=15]{.underline}

-   Run `.sh`

    ``` bash
    sbatch --wrap=" source ./.sh"
    ```

    -   run .R

        -   `args.txt` To **modify** ( let a line break at the end of the file ! )

        -   `one_sample_RAW_2.Rmd`

-   Run `line_1.sh`

    ``` bash
    sbatch --wrap=" source ./line_1.sh"
    ```

    -   `run_filter.R`

        -   `one_line_NORM_1.Rmd` : Default [pc=60]{.underline}

-   Run `line_2.sh`

````         
``` bash
sbatch --wrap=" source ./line_2.sh"
```

-   `run_clust.R`

    -   `one_line_NORM_2.Rmd` : **Choose** pc *(it will be common for all samples, with SCT the most PC the better)*
````

-   Run `line_3.sh`

    ``` bash
    sbatch --wrap=" source ./line_3.sh"
    ```

    -   `run_3.R`

        -   `res.txt` To **modify** ( let a line break at the end of the file ! )

        -   `one_line_NORM_3.Rmd`

        -   *(run mca)*

        -   `one_line_NORM_4.Rmd`

#### --- Metabo only ---

### b) All samples

#### --- All genes ---

#### --- Metabo only ---

# Scripts content

-   `one_sample_RAW_1.Rmd` : Filter one sample based on nCount, nFeature_RNA and percent Mitochondrion.
    -   *QC_sample.rds* (RAW)
    -   *Doublet_sample.rds* (NORM)
-   `one_sample_RAW_2.Rmd`
    -   *Filtered_sample.rds* (RAW)
-   `one_line_NORM_1.Rmd`
    -   ***Norm_line.rds***
    -   *Norm_sample.rds*
-   `one_line_NORM_2.Rmd`
    -   *line_FindNeighbors.rds*
-   `one_line_NORM_3.Rmd`
    -   ***line_cluster.markers.rds*** (in html directory ! )
    -   *line_Markers.rds*
-   ( run_3.R )
    -   *line_Mca.rds*
    -   ***line_HGT_brain_signatures.rds***
-   `one_line_NORM_4.Rmd`
    -   ***line_Clusters.rds***

------------------------------------------------------------------------

### All **genes**

-   Run `run_QC_all_samples.sh`

    -   `create_seurat_samples.R` : Change the samples names if needed.

    -   `run_QC_all_samples.R`

        -   `one_sample_RAW_1.Rmd` : Choose `x1`, `x2` and `x3` for the filtering. (default: 500, 220, 10)

``` bash
source ./Analysis/Scripts/Per_sample/run_QC_all_samples.sh
```

sbatch -c 64 --mem-per-cpu=2048 --wrap="Rscript 'all_samples_metabo_2.R'"

``` bash
module load system/R-4.2.3_Miniconda3
Rscript -e ./Analysis/Scripts/Per_sample/run_QC_all_samples_2.R
```

-   Run `create_seurat_all.R`
    -   Be careful of the **path** of all samples (and all the names if script used for different data)

``` bash
module load system/singularity-3.7.3
singularity pull docker://
sbatch -c 16 --mem-per-cpu=2048 --wrap="singularity exec seurat_latest.sif Rscript 'create_seurat_all.R'"
```

-   Run `all_samples_RAW.Rmd`

``` bash

mkdir -p ./Analysis/Scripts/All_genes/results/html
Rscript -e  "rmarkdown::render('./Analysis/Scripts/All_genes/all_samples_RAW.Rmd', output_dir = './Analysis/Scripts/All_genes/results/html/')"
```

-   Run `sctransform.sh`

    -   Careful to **which object** you're working on (metabo or all)

-   Run `all_samples_NORM_1.Rmd`

-   Run `pca.sh`

    -   Careful to **which object** you're working on (metabo or all)
    -   See if the number of components is enough (default: 40)

-   Run `all_samples_NORM_2.Rmd`

    -   **Modify** the number of pc (**/!\\** as sctransform is used, you can use more pc than needed, even recommended)
    -   **Modify** the number of thread available (default: 8)

-   Run `all_samples_NORM_3.Rmd`

    -   **Modify** the resolution to define clusters. (default: 0.5)

-   Run `mca.sh`

    -   Careful to **which object** you're working on (metabo or all)

-   Run `all_samples_NORM_4.Rmd`

    -   Choose the Panglao database (default: humain brain)
    -   Choose the enrichR database (default: KEGG_2021_Human)
    -   Choose the maximum nb of genes to use as input to enrichR (default: 200)

### Metabo genes

-   Create *all_metabo.tsv*
-   Run `all_samples_metabo_1.Rmd`
-   Run `pca.sh`
    -   Careful to **which object** you're working on (metabo or all)
    -   See if the number of components is enough (default: 40)
-   Run `all_samples_metabo_2.Rmd`
    -   **Modify** the number of pc (**/!\\** as sctransform is used, you can use more pc than needed, even recommended)
    -   **Modify** the number of thread available (default: 8)
-   Run `all_samples_metabo_3.Rmd` *run it on notebook format to adapt cluster_colors*
    -   **Modify** `cluster_colors` (metabo)
    -   **Modify** the resolution to define clusters. (default: 0.5)
-   Run `mca.sh`
    -   Careful to **which object** you're working on (metabo or all)
-   Create *cluster_markers_signatures.csv*
-   Run `all_samples_metabo_4.Rmd`
    -   **Modify** `cluster_colors` (metabo)
    -   Choose the Panglao database (default: humain brain)
    -   Choose the enrichR database (default: KEGG_2021_Human)
    -   Choose the maximum nb of genes to use as input to enrichR (default: 200)

# Scripts content

-   `create_seurat_samples.R`
    -   mtx\_*project*.rds
-   `run_QC_all_samples.sh` (run `run_QC_all_samples.sh` and `one_sample_RAW.Rmd` for each sample): Filter one sample based on nFeature_RNA and percent Mitochondrion.
    -   Filtered\_*project*.rds
-   `all_samples_RAW.Rmd` : Combined all raw filtered samples (+ remove zero count genes)
    -   Raw_merged_seurat.rds
-   `sctransform.sh` : normalize and scale data
    -   sctransform.rds
-   `all_samples_NORM_1.Rmd` : run PCA
    -   Norm_merged_seurat.rds
-   `pca.sh` : Plot JackStraw, Elbow and Heatmap to determine how many components are relevant.
    -   JackStrawPlot.pdf
    -   ElbowPlot.pdf
    -   DimHeatmap.pdf
-   `all_samples_NORM_2.Rmd` : run UMAP, tSNE and clustree (to visualize resolution for clusters)
    -   merged_seurat_FindNeighbors.rds
-   `all_samples_NORM_3.Rmd` : Define clusters and Markers
    -   merged_seurat_Markers.rds
-   `mca.R` : run MCA
    -   merged_seurat_Mca.rds
-   `all_samples_NORM_4.Rmd` : Cell-type annotation and Ontology
    -   merged_seurat_Clusters.rds