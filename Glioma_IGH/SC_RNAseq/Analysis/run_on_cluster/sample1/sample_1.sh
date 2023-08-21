start=$(date +%s)
module load system/singularity-3.7.3

samples=$(ls ./results/rds/mtx/ | grep '^mtx' | grep 'prolif.rds$' | sed 's/mtx_//g' | sed 's/_prolif.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                      PCA                                        "
echo -e "\n                ---------------------------------------------------------------------------------"

for samp in ${samples}; do
    mkdir -p ./results/html/$samp
    mkdir -p ./results/rds/samples/$samp
    mkdir -p ./out/$samp
    #rm ./out/$samp/*

    ## diff
    diff="${samp}_diff"
    echo -e "Running $diff..."
    sbatch -c 4 --mem-per-cpu=2048 --output="./out/${samp}/${diff}_first.out" --wrap="singularity exec r.sif Rscript run_clust.R first $diff"
    sbatch -c 4 --mem-per-cpu=2048 --output="./out/${samp}/${diff}_cond1.out" --wrap="singularity exec r.sif Rscript run_clust.R samp $diff"

    ## prolif
    prolif="${samp}_prolif"
    echo -e "Running $prolif..."
    sbatch -c 4 --mem-per-cpu=2048 --output="./out/${samp}/${prolif}_first.out" --wrap="singularity exec r.sif Rscript run_clust.R first $prolif"
    sbatch -c 4 --mem-per-cpu=2048 --output="./out/${samp}/${prolif}_cond1.out" --wrap="singularity exec r.sif Rscript run_clust.R samp $prolif"
    
    echo -e "Running $samp ..."
    ## both
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/$samp/${samp}_line1.out" --wrap="singularity exec r.sif Rscript run_clust.R line $samp"
done

## sample_FindNeighbors.rds
## line_sample_FindNeighbors.rds
## sample_NORM_0.html
## sample_NORM_1.html
## line_sample_NORM_1.html

for samp in ${samples}; do
    while [[ ! -f ./results/html/$samp/$samp_diff_NORM_0.html || ! -f ./results/html/$samp/$samp_prolif_NORM_0.html || \
             ! -f ./results/html/$samp/$samp_diff_NORM_1.html || ! -f ./results/html/$samp/$samp_prolif_NORM_1.html || \
             ! -f ./results/html/$samp/line_$samp_diff_NORM_1.html || ! -f ./results/html/$samp/line_$samp_prolif_NORM_1.html ]];do
        sleep 1; done
done

echo -e "\n                --------------------                 Done !                  --------------------\n"

echo -e " Next step : Choose pc & run sample_2.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"
