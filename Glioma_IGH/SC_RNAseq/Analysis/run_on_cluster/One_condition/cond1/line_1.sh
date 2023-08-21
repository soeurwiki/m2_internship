start=$(date +%s)
#module load system/singularity-3.7.3
module load containers/singularity/3.9.9 

samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                               Normalisation & PCA                               "
echo -e "\n                ---------------------------------------------------------------------------------"

for samp in ${samples}; do
    mkdir -p ./Lucille/results/rds/samples
    mkdir -p ./Lucille/results/html/${samp}
    mkdir -p ./Lucille/out/${samp}

    echo -e "Running $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp1.out" --wrap="singularity exec r.sif Rscript run_norm.R all $samp"
done


## Norm_samp_diff.rds
## Norm_samp_prolif.rds
## Norm_samp.rds
## one_samp_NORM_1.html

for samp in ${samples}; do
    while [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_1.html ]];do
        sleep 1; done
done

#######################
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )
for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_1.html ]];then
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp1.out" --wrap="singularity exec r.sif Rscript run_norm.R all $samp"
    fi
done
#######################


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                   PCA [metabo]                                  "
echo -e "\n                ---------------------------------------------------------------------------------"

for samp in ${samples}; do
    echo -e "Running $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo1.out" --wrap="singularity exec r.sif Rscript run_norm.R metabo $samp"
done

## samp_metabo.rds
## one_sample_metabo_1.html

for samp in ${samples}; do
    while [[ ! -f ./Lucille/results/html/${samp}/${samp}_metabo_1.html ]];do
        sleep 1; done
done

#######################
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )
for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_metabo_1.html ]];then
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo1.out" --wrap="singularity exec r.sif Rscript run_norm.R metabo $samp"
    fi
done
#######################

echo -e "\n                --------------------                 Done !                  --------------------\n"


echo -e " Next step : Choose pc & run samp_2.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"