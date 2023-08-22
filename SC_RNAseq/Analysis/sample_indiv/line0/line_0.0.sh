start=$(date +%s)
module load system/singularity-3.7.3
samples=("BT1" "BT2" "BT54" "BT88" "LGG85" "LGG275" "LGG336" "LGG349")

mkdir -p ./out/mtx
mkdir -p ./results/html/samples
mkdir -p ./results/rds/mtx

echo " Please check the .out files during each step to see if the script was not halted at one moment ... "


echo -e "\n                ---------------------------------------------------------------------------------\n"
echo -e "                                            Creation of seurat objects                           "
echo -e "\n                ---------------------------------------------------------------------------------"
    
for line in ${samples}; do
    diff="$line_diff"
    prolif="$line_prolif"

    echo -e " The mtx rds files for $line are not present in  ./results/rds/mtx/. \n
                Creation of the corresponding rds ... \n "

    sbatch -c 4 --output=./out/mtx/"$diff.out" --wrap="singularity exec r.sif Rscript 'create_seurat_samples.R $diff'"
    sbatch -c 4 --output=./out/mtx/"$prolif.out" --wrap="singularity exec r.sif Rscript 'create_seurat_samples.R $prolif'"
done

## mtx_project.rds

for line in ${samples}; do
    while [[ ! -f ./results/rds/mtx/mtx_${line}_prolif.rds && ! -f ./results/rds/mtx/mtx_${line}_diff.rds ]];do
        sleep 1; done
done

echo -e "\n                --------------------                Done !                 --------------------\n"

end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"

echo -e " Next step : run rmd_2.sh "