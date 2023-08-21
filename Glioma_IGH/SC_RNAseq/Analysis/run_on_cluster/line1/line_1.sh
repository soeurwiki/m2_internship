start=$(date +%s)
module load system/singularity-3.7.3

samples=$(ls ./results/rds/mtx/ | grep '^mtx' | grep 'prolif.rds$' | sed 's/mtx_//g' | sed 's/_prolif.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                               Normalisation & PCA                               "
echo -e "\n                ---------------------------------------------------------------------------------"

for line in ${samples}; do
    mkdir -p ./results/rds/samples
    mkdir -p ./results/html/${line}
    mkdir -p ./out/${line}

    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/line1.out" --wrap="singularity exec r.sif Rscript run_norm.R all $line"
done


## Norm_line_diff.rds
## Norm_line_prolif.rds
## Norm_line.rds
## one_line_NORM_1.html

for line in ${samples}; do
    while [[ ! -f ./results/html/${line}/${line}_NORM_1.html ]];do
        sleep 1; done
done


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                   PCA [metabo]                                  "
echo -e "\n                ---------------------------------------------------------------------------------"

for line in ${samples}; do
    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/metabo1.out" --wrap="singularity exec r.sif Rscript run_norm.R metabo $line"
done

## line_metabo.rds
## one_sample_metabo_1.html

for line in ${samples}; do
    while [[ ! -f ./results/html/${line}/${line}_metabo_1.html ]];do
        sleep 1; done
done

echo -e "\n                --------------------                 Done !                  --------------------\n"


echo -e " Next step : Choose pc & run line_2.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"