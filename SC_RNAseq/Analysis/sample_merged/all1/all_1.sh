start=$(date +%s)
#module load system/singularity-3.7.3
module load containers/singularity/3.9.9 

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                        PCA                                      "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 32 --mem-per-cpu=2048 --output="./out/all1.out" --wrap="singularity exec r.sif Rscript run_pca.R all"

## Norm_all_samples.rds
## all_samples_NORM_1.html

while [[ ! -f ./results/html/all_samples_NORM_1.html ]];do
        sleep 1; done


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                   PCA [metabo]                                  "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 32 --mem-per-cpu=2048 --output="./out/all_metabo1.out" --wrap="singularity exec r.sif Rscript run_pca.R metabo"

## all_metabo.rds
## all_samples_metabo_1.html

while [[ ! -f ./results/html/all_samples_metabo_1.html ]];do
        sleep 1; done


echo -e "\n                --------------------                 Done !                  --------------------\n"


echo -e " Next step : Choose pc & run all_2.sh "
echo -e " Next step : Choose pc & run all_5.sh "

end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"