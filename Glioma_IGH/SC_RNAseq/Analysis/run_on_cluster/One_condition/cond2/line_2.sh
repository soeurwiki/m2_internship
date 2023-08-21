## ----- Create & organize workspace -----  ##
start=$(date +%s)
#module load system/singularity-3.7.3
module load containers/singularity/3.9.9
 

samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                 Defining clusters                               "
echo -e "\n                ---------------------------------------------------------------------------------"


for samp in ${samples}; do
    echo -e "Running $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp2.out" --wrap="singularity exec r.sif Rscript run_clust.R all $samp"
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo2.out" --wrap="singularity exec r.sif Rscript run_clust.R metabo $samp"
done

## samp_FindNeighbors.rds
## one_samp_NORM_2.html

## samp_FindNeighbors_metabo.rds
## samp_metabo_2.html

for samp in ${samples}; do
    while [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_2.html || ! -f ./Lucille/results/html/${samp}/${samp}_metabo_2.html ]];do
        sleep 1; done
done

#######################
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )
for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_2.html ]];then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp2.out" --wrap="singularity exec r.sif Rscript run_clust.R all $samp"; fi
done
for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_metabo_2.html ]];then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo2.out" --wrap="singularity exec r.sif Rscript run_clust.R metabo $samp"; fi
done
#######################


echo -e "\n                --------------------                Done !                 --------------------\n"



echo -e " Next step : Choose resolution & run samp_3.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"
