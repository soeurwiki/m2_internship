## ----- Create & organize workspace -----  ##
start=$(date +%s)
module load system/singularity-3.7.3

echo -e "\n                ---------------------------------------------------------------------------------\n"
echo -e "                                                 Defining clusters                               "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 64 --mem-per-cpu=2048 --output="./out/norm2.out" --wrap="singularity exec r.sif Rscript run_clust.R "

## merged_seurat_FindNeighbors.rds
## all_samples_NORM_2.html

while [[ ! -f ./results/html/all_samples_NORM_2.html || ! -f ./results/rds/merged_seurat_FindNeighbors.rds ]];do
    sleep 1; done


echo -e "\n                --------------------                Done !                 --------------------\n"



echo -e " Next step : Choose resolution & run rmd_4.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                 Defining markers                                "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 64 --mem-per-cpu=2048 --output="./out/norm3.out" --wrap="singularity exec r.sif Rscript run_4.R mark"

## merged_seurat_Markers.rds
## all_samples_NORM_3.html

while [[ ! -f ./results/html/all_samples_NORM_3.html || ! -f ./results/rds/merged_seurat_Markers.rds ]]; do
    sleep 1; done


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                       MCA                                       "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 64 --mem-per-cpu=2048 --output="./out/mca.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_4.R mca "

## merged_seurat_Mca.rds

while [[ ! -f ./results/rds/merged_seurat_Mca.rds || ! -f ./results/rds/HGT_brain_signatures.rds ]]; do
    sleep 1; done


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                    Annotation                                   "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 64 --mem-per-cpu=2048 --output="./out/norm4.out" --wrap="singularity exec r.sif Rscript run_4.R annot"

## merged_seurat_Clusters.rds
## all_samples_NORM_4.html

while [[ ! -f ./results/html/all_samples_NORM_4.html || ! -f ./results/rds/merged_seurat_Clusters.rds ]]; do
    sleep 1; done

echo -e "\n                --------------------                Done !                 --------------------\n"


echo -e " Next step : run create_objects.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"

mkdir -p ./results/rds/samples

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                      Creation of seurat object [metabo]                         "
echo -e "\n                ---------------------------------------------------------------------------------"
    
sbatch -c 64 --mem-per-cpu=2048 --output=./out/"create_metabo.out" --wrap="singularity exec r.sif Rscript create_sample_norm.R metabo "

## merged_seurat_metabo1.rds

while [[ ! -f ./results/rds/merged_seurat_metabo1.rds ]]; do
    sleep 1; done

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                       Creation of samples object [norm]                         "
echo -e "\n                ---------------------------------------------------------------------------------"
    
sbatch -c 64 --mem-per-cpu=2048 --output=./out/"norm_samp.out" --wrap="singularity exec r.sif Rscript create_sample_norm.R samples "

## norm_project.rds

while [[ ! -f ./results/rds/samples/norm_BT1_prolif.rds ]]; do
    sleep 1; done

while [[ ! $(ls ./results/rds/samples/norm_* | wc -l) -eq 16  ]];do
    sleep 1; done

echo -e "\n                --------------------                Done !                 --------------------\n"

end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"

echo -e " Next step : run metabo scripts ! "