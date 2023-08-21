## ----- Create & organize workspace -----  ##
start=$(date +%s)
module load system/singularity-3.7.3

samples=$(ls ./results/rds/mtx/ | grep '^mtx' | grep 'prolif.rds$' | sed 's/mtx_//g' | sed 's/_prolif.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                 Defining markers                                "
echo -e "\n                ---------------------------------------------------------------------------------"


while read args; do 
    line=$( echo ${args} | head -n1 | cut -d " " -f1)
    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/line3.out" --wrap="singularity exec r.sif Rscript run_3.R all_mark $args"
done < res.txt

while read args; do 
    line=$( echo ${args} | head -n1 | cut -d " " -f1)
    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/metabo3.out" --wrap="singularity exec r.sif Rscript run_3.R metabo_mark $args"
done < res_metabo.txt

## line_Markers.rds
## one_line_NORM_3.html

## line_Markers_metabo.rds
## one_line_metabo_3.html

for line in ${samples}; do
    while [[ ! -f ./results/html/${line}/${line}_NORM_3.html || ! -f ./results/html/${line}/${line}_metabo_3.html ]];do
        sleep 1; done
done


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                       MCA                                       "
echo -e "\n                ---------------------------------------------------------------------------------"

for line in ${samples}; do
    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/mca.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_3.R all_mca $line"
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/mca_metabo.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_3.R metabo_mca $line"
done

## line_Mca.rds
## line_Mca_metabo.rds

for line in ${samples}; do
    while [[ ! -f ./results/rds/samples/${line}_Mca.rds || ! -f ./results/rds/samples/${line}_Mca_metabo.rds ]];do
        sleep 1; done
done

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                    Annotation                                   "
echo -e "\n                ---------------------------------------------------------------------------------"

for line in ${samples}; do
    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/line4.out" --wrap="singularity exec r.sif Rscript run_3.R all_annot $line"
done

## line_Clusters.rds
## one_line_NORM_4.html

for line in ${samples}; do
    while [[ ! -f ./results/html/${line}/${line}_NORM_4.html ]];do
        sleep 1; done
done

for line in ${samples}; do
    echo -e "Running $line ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./out/${line}/metabo4.out" --wrap="singularity exec r.sif Rscript run_3.R metabo_annot $line"
done

## line_Clusters_metabo.rds
## one_line_metabo_4.html

for line in ${samples}; do
    while [[ ! -f ./results/html/${line}/${line}_metabo_4.html ]];do
        sleep 1; done
done

done
echo -e "\n                --------------------                Done !                 --------------------\n"


end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"