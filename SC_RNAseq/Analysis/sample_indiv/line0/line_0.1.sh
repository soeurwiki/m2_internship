#   name : line_0.1.sh
#
#   Author (2023)  Safiya ATIA

module load containers/singularity/3.9.9

## ----- Create & organize workspace -----  ##
start=$(date +%s)

mkdir -p ./out
mkdir -p ./results/html/samples
mkdir -p ./results/rds/mtx
mkdir -p ./results/rds/samples

echo " Please check the .out files during each step to see if the script was not halted at one moment ... "


samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                   Quality Check                                 "
echo -e "\n                ---------------------------------------------------------------------------------"

mkdir -p ./results/rds/samples

for samp in ${samples}; do
    mkdir -p ./out/${samp}

    echo -e " QC for $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${samp}/QC_${samp}.out" --wrap="singularity exec r.sif Rscript run_filter.R qc $samp "

done

#  Outputs:
## QC_project.rds
## sample_RAW_qc.html

for samp in ${samples}; do
    while [[ ! -f ./results/html/samples/${samp}_RAW_qc.html ]];do
        sleep 1; done
    echo -e " QC for $samp done."
done

########### [Copy&paste] if code is stuck in loop ###########
#############################################################
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

## check samples
for samp in ${samples}; do
    if [[ ! -f ./results/html/samples/${samp}_RAW_qc.html ]];then
        echo "missing $samp"; fi
done

## run missing samples
for samp in ${samples}; do
    if [[ ! -f ./results/html/samples/${samp}_RAW_qc.html ]]; then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./out/${samp}/QC_${samp}.out" --wrap="singularity exec r.sif Rscript run_qc.R $samp "
    fi
done
#############################################################
#############################################################

echo -e "\n                --------------------                Done !                 --------------------\n"

end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"

echo -e " Next step : modify args.txt (pc + doublet_rate) & run rmd_2.sh "