## ----- Create & organize workspace -----  ##
start=$(date +%s)
#module load system/singularity-3.7.3
module load containers/singularity/3.9.9 

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                 Defining clusters                               "
echo -e "\n                ---------------------------------------------------------------------------------"


sbatch -c 32 --mem-per-cpu=2048 --output="./out/all2.out" --wrap="singularity exec r.sif Rscript run_clust.R all"
sbatch -c 32 --mem-per-cpu=2048 --output="./out/all_metabo2.out" --wrap="singularity exec r.sif Rscript run_clust.R metabo"

## all_samples_FindNeighbors.rds
## all_samples_NORM_2.html

## all_samples_FindNeighbors_metabo.rds
## all_samples_metabo_2.html

while [[ ! -f ./results/html/all_samples_NORM_2.html || ! -f ./results/html/all_samples_metabo_2.html ]];do
    sleep 1; done

echo -e "\n                --------------------                Done !                 --------------------\n"


echo -e " Next step : Choose resolution & run all_3.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"
