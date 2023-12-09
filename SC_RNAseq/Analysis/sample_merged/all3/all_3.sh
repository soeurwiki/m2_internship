#   name : all_3.sh
#
#   Author (2023)  Safiya ATIA

module load containers/singularity/3.9.9 

start=$(date +%s)

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                 Defining markers                                "
echo -e "\n                ---------------------------------------------------------------------------------"


sbatch -c 32 --mem-per-cpu=2048 --output="./out/all3.out" --wrap="singularity exec r.sif Rscript run_annot.R all_mark"
sbatch -c 32 --mem-per-cpu=2048 --output="./out/all_metabo3.out" --wrap="singularity exec r.sif Rscript run_annot.R metabo_mark"

#  Outputs:

## all_samples_Markers.rds
## all_samples_NORM_3.html

## all_samples_Markers_metabo.rds
## all_samples_metabo_3.html

while [[ ! -f ./results/html/all_samples_NORM_3.html || ! -f ./results/html/all_samples_metabo_3.html ]];do
    sleep 1; done


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                       MCA                                       "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 32 --mem-per-cpu=2048 --output="./out/all_mca.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_annot.R all_mca"
sbatch -c 32 --mem-per-cpu=2048 --output="./out/all_mca_metabo.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_annot.R metabo_mca"

#  Outputs:
## all_samples_Mca.rds
## all_samples_Mca_metabo.rds

while [[ ! -f ./results/rds/all_samples_Mca.rds || ! -f ./results/rds/samples/all_samples_Mca_metabo.rds ]];do
    sleep 1; done


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                    Annotation                                   "
echo -e "\n                ---------------------------------------------------------------------------------"

sbatch -c 64 --mem-per-cpu=2048 --output="./out/all4.out" --wrap="singularity exec r.sif Rscript run_annot.R all_annot"

#  Outputs:
## all_samples_Clusters.rds
## all_samples_NORM_4.html

while [[ ! -f ./results/html/all_samples_NORM_4.html ]];do
    sleep 1; done

sbatch -c 32 --mem-per-cpu=2048 --output="./out/all_metabo4.out" --wrap="singularity exec r.sif Rscript run_annot.R metabo_annot"

#  Outputs:
## all_samples_Clusters_metabo.rds
## all_samples_metabo_4.html

while [[ ! -f ./results/html/all_samples_metabo_4.html ]];do
    sleep 1; done

echo -e "\n                --------------------                Done !                 --------------------\n"


end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"