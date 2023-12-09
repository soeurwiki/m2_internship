#   name : line_0.2.sh
#
#   Author (2023)  Safiya ATIA

## ----- Create & organize workspace -----  ##
start=$(date +%s)

mkdir -p ./out
mkdir -p ./results/html/samples
mkdir -p ./results/rds/mtx
mkdir -p ./results/rds/samples

echo " Please check the .out files during each step to see if the script was not halted at one moment ... "

# https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                     Filtering                                   "
echo -e "\n                ---------------------------------------------------------------------------------"

while read args; do 
    samp=$( echo ${args} | head -n1 | cut -d " " -f1)
    echo -e "\n Filtering of $samp ..."
    
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${samp}/Filter_${samp}.out" --wrap="singularity exec r.sif Rscript run_filter.R filter $args"

done < args.txt  

#  Outputs:
## Filtered_project.rds
## sample_RAW_filtered.html

for samp in ${samples}; do
    while [[ ! -f ./results/html/samples/${samp}_RAW_filtered.html ]];do
        sleep 1; done
    echo -e " Filtering of $samp done."
done

echo -e " Filtering of all samples done !"

echo -e "\n                --------------------                 Done !                  --------------------\n"

end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"