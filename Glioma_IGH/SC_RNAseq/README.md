# Scripts order

## Processing

*On the cluster genotoul*

Run `get_files.sh` \* Run `run_pipeline.sh`

```         
source ./Processing/get_files.sh
source ./Processing/run_pipeline.sh
```

## Analysis

### All genes

-   Run `one_sample_RAW.Rmd` for **each** sample.
    -   **Modify** `x1`, `x2` and `x3` for the filtering.
-   Run `all_sample_RAW.Rmd`
    -   Be careful of the **path** of all samples (and all the names if script used for different data)
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
-   Run `all_samples_NORM_3.Rmd`
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

-   `one_sample_RAW.Rmd` : Filter one sample based on nFeature_RNA and percent Mitochondrion.
    -   QC\_*project*.rds
-   `all_sample_RAW.Rmd` : Combined all raw filtered samples (+ remove zero count genes)
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