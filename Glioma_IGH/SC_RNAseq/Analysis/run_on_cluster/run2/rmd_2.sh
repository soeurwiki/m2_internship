start=$(date +%s)
module load containers/singularity/3.9.9 

# https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                     Filtering                                   "
echo -e "\n                ---------------------------------------------------------------------------------"

while read samples; do 
    samp=$( echo ${samples} | head -n1 | cut -d " " -f1)
    echo " Filtering of "$samp" ..."
    
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/Filter_$samp.out" --wrap="singularity exec r.sif Rscript run_filter.R $samples "

done < args.txt  


## Filtered_project.rds
## sample_RAW_filtered.html

for samp in ${samples}; do
    while [[ ! -f ./results/html/samples/${samp}_RAW_filtered.html ]];do
        sleep 1; done
    echo -e " Filtering for $samp done."
done

########### [Copy&paste] if code is stuck in loop ###########
#############################################################
module load containers/singularity/3.9.9 
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

## check samples
for samp in ${samples}; do
    if [[ ! -f ./results/html/samples/${samp}_RAW_filtered.html ]];then
        echo "missing $samp"; fi
done

## run missing samples
for samp in ${samples}; do
    if [[ ! -f ./results/html/samples/${samp}_RAW_filtered.html ]]; then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./out/Filter_$samp.out" --wrap="singularity exec r.sif Rscript run_filter.R $samples "
    fi
done
#############################################################
#############################################################

echo " Filtering of all samples done !"

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                  Merging samples                                "
echo -e "\n                ---------------------------------------------------------------------------------"


sbatch -c 64 --mem-per-cpu=2048 --output="./out/RAW_all.out" --wrap="singularity exec r.sif Rscript run_merge.R "

while [[ ! -f ./results/html/all_samples_RAW.html ]];do
    sleep 1; done

echo " Seurat object merging all samples created !"

## Raw_merged_seurat.rds
## all_samples_RAW.html


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                   Normalisation                                 "
echo -e "\n                ---------------------------------------------------------------------------------"


sbatch -c 64 --mem-per-cpu=2048 --output="./out/norm1.out" --wrap="singularity exec r.sif Rscript run_norm.R "

## Norm_merged_seurat.rds
## all_samples_NORM_1.html

while [[ ! -f ./results/html/all_samples_NORM_1.html ]];do
    sleep 1; done

echo " Normalisation was done !"


echo -e "\n                --------------------                Done !                 --------------------\n"

echo -e " Next step : Choose pc & run rmd_3.sh "
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"