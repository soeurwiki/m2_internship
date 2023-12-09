#   name : line_0.0.sh
#
#   Author (2023)  Safiya ATIA

module load containers/singularity/3.9.9

start=$(date +%s)

## samples names
samples=("BT1" "BT2" "BT54" "BT88" "LGG85" "LGG275" "LGG336" "LGG349")

## ----- Create & organize workspace -----  ##
mkdir -p ./out/mtx
mkdir -p ./results/html/samples
mkdir -p ./results/rds/mtx

echo " Please check the .out files during each step to see if the script was not halted at one moment ... "

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                            Creation of seurat objects                           "
echo -e "\n                ---------------------------------------------------------------------------------"
    
for line in ${samples}; do
    diff="$line_diff"
    prolif="$line_prolif"

    sbatch -c 4 --output=./out/mtx/"$diff.out" --wrap="singularity exec r.sif Rscript 'create_seurat_samples.R $diff'"
    sbatch -c 4 --output=./out/mtx/"$prolif.out" --wrap="singularity exec r.sif Rscript 'create_seurat_samples.R $prolif'"
done

# Output: 
## mtx_project.rds

for line in ${samples}; do
    while [[ ! -f ./results/rds/mtx/mtx_${line}_prolif.rds && ! -f ./results/rds/mtx/mtx_${line}_diff.rds ]];do
        sleep 1; done
done

echo -e "\n                --------------------                Done !                 --------------------\n"

end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"

echo -e " Next step : run rmd_2.sh "