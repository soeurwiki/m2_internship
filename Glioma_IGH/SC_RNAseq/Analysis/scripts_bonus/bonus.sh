## ----- Create & organize workspace -----  ##
module load system/singularity-3.7.3

sbatch -c 64 --mem-per-cpu=3072 --output="./out/plots.out" --wrap="singularity exec r.sif Rscript bonus.R"

while [[ ! -f ./results/html/plots.html ]];do
    sleep 1; done
echo -e "\n                --------------------                Done !                 --------------------\n"
