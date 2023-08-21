## ----- Create & organize workspace -----  ##
start=$(date +%s)
module load system/singularity-3.7.3

samples=$(ls ./results/rds/mtx/ | grep '^mtx' | grep 'prolif.rds$' | sed 's/mtx_//g' | sed 's/_prolif.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                 Defining clusters                               "
echo -e "\n                ---------------------------------------------------------------------------------"


for line in ${samples}; do
    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/line2.out" --wrap="singularity exec r.sif Rscript run_clust.R all $line"
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/metabo2.out" --wrap="singularity exec r.sif Rscript run_clust.R metabo $line"
done

## line_FindNeighbors.rds
## one_line_NORM_2.html

## line_FindNeighbors_metabo.rds
## line_metabo_2.html

for line in ${samples}; do
    while [[ ! -f ./results/html/${line}/${line}_NORM_2.html || ! -f ./results/html/${line}/${line}_metabo_2.html ]];do
        sleep 1; done
done

echo -e "\n                --------------------                Done !                 --------------------\n"



echo -e " Next step : Choose resolution & run line_3.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"
