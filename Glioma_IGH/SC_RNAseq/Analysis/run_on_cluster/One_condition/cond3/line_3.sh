## ----- Create & organize workspace -----  ##
start=$(date +%s)
#module load system/singularity-3.7.3
module load containers/singularity/3.9.9 

samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                 Defining markers                                "
echo -e "\n                ---------------------------------------------------------------------------------"


while read args; do 
    samp=$( echo ${args} | head -n1 | cut -d " " -f1)
    echo -e "Running $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp3.out" --wrap="singularity exec r.sif Rscript run_3.R all_mark $args"
done < res.txt

while read args; do 
    samp=$( echo ${args} | head -n1 | cut -d " " -f1)
    echo -e "Running $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo3.out" --wrap="singularity exec r.sif Rscript run_3.R metabo_mark $args"
done < res_metabo.txt

## samp_Markers.rds
## one_samp_NORM_3.html

## samp_Markers_metabo.rds
## one_samp_metabo_3.html

for samp in ${samples}; do
    while [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_3.html || ! -f ./Lucille/results/html/${samp}/${samp}_metabo_3.html ]];do
        sleep 1; done
done

#######################
module load containers/singularity/3.9.9 
while read args; do 
    samp=$( echo ${args} | head -n1 | cut -d " " -f1)
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_3.html ]];then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp3.out" --wrap="singularity exec r.sif Rscript run_3.R all_mark $args"
    fi
done < res.txt


while read args; do 
    samp=$( echo ${args} | head -n1 | cut -d " " -f1)
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_metabo_3.html ]];then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo3.out" --wrap="singularity exec r.sif Rscript run_3.R metabo_mark $args"
    fi
done < res_metabo.txt
#######################

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                       MCA                                       "
echo -e "\n                ---------------------------------------------------------------------------------"

for samp in ${samples}; do
    echo -e "Running $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/mca.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_3.R all_mca $samp"
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/mca_metabo.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_3.R metabo_mca $samp"
done

## samp_Mca.rds
## samp_Mca_metabo.rds

for samp in ${samples}; do
    while [[ ! -f ./Lucille/results/rds/samples/${samp}_Mca.rds || ! -f ./Lucille/results/rds/samples/${samp}_Mca_metabo.rds ]];do
        sleep 1; done
done

########### [Copy&paste] if code is stuck in loop ###########
#############################################################
module load containers/singularity/3.9.9 
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

## run missing samples
for samp in ${samples}; do
    if [[ -f ./Lucille/results/rds/samples/${samp}_Markers.rds && ! -f ./Lucille/results/rds/samples/${samp}_Mca.rds ]]; then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/mca.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_3.R all_mca $samp"
    fi
done
for samp in ${samples}; do
    if [[ -f ./Lucille/results/rds/samples/${samp}_Markers_metabo.rds && ! -f ./Lucille/results/rds/samples/${samp}_Mca_metabo.rds ]]; then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/mca_metabo.out" --wrap="singularity exec cellid_0.1.0.sif Rscript run_3.R metabo_mca $samp"
    fi
done

## check samples
for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/rds/samples/${samp}_Mca.rds && ! -f ./Lucille/results/rds/samples/${samp}_Mca_metabo.rds ]];then
        echo "missing $samp"; fi
done
#############################################################
#############################################################


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                    Annotation                                   "
echo -e "\n                ---------------------------------------------------------------------------------"

for samp in ${samples}; do
    echo -e "Running $samp ..."
    sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp4.out" --wrap="singularity exec r.sif Rscript run_3.R all_annot $samp"
done

## samp_Clusters.rds
## one_samp_NORM_4.html

for samp in ${samples}; do
    while [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_4.html ]];do
        sleep 1; done
done

########### [Copy&paste] if code is stuck in loop ###########
#############################################################
module load containers/singularity/3.9.9 
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

## run missing samples
for samp in ${samples}; do
    if [[ -f ./Lucille/results/rds/samples/${samp}_Mca.rds && ! -f ./Lucille/results/html/${samp}/${samp}_NORM_4.html ]]; then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/samp4.out" --wrap="singularity exec r.sif Rscript run_3.R all_annot $samp"
    fi
done

## check samples
for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_NORM_4.html ]];then
        echo "missing $samp"; fi
done
#############################################################
#############################################################

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                                Annotation (metabo)                              "
echo -e "\n                ---------------------------------------------------------------------------------"

for samp in ${samples}; do
    #if [[ -f ./Lucille/results/html/${samp}/${samp}_NORM_4.html ]];then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo4.out" --wrap="singularity exec r.sif Rscript run_3.R metabo_annot $samp"
    #fi
done

## samp_Clusters_metabo.rds
## one_samp_metabo_4.html

for samp in ${samples}; do
    while [[ ! -f ./Lucille/results/html/${samp}/${samp}_metabo_4.html ]];do
        sleep 1; done
done

########### [Copy&paste] if code is stuck in loop ###########
#############################################################
module load containers/singularity/3.9.9 
samples=$(ls ./results/rds/mtx/ | grep '^mtx' | sed 's/mtx_//g' | sed 's/.rds//g' )

## check samples
for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_metabo_4.html ]];then
        echo "missing $samp"; fi
done

for samp in ${samples}; do
    if [[ ! -f ./Lucille/results/html/${samp}/${samp}_metabo_4.html ]];then
        echo -e "Running $samp ..."
        sbatch -c 8 --mem-per-cpu=2048 --output="./Lucille/out/${samp}/metabo4.out" --wrap="singularity exec r.sif Rscript run_3.R metabo_annot $samp"
    fi
done
#############################################################
#############################################################


echo "  ____                     _ ";
echo " |  _ \  ___  _ __   ___  | |";
echo " | | | |/ _ \| '_ \ / _ \ | |";
echo " | |_| | (_) | | | |  __/ |_|";
echo " |____/ \___/|_| |_|\___| (_)";
echo "                             ";


end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"