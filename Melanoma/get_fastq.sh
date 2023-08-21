module load bioinfo/sratoolkit.3.0.0
start=$(date +%s)

mkdir -p MELANOMA/sc/
mkdir -p MELANOMA/bulk/

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                       Retrieving SC RNA-seq files                               "
echo -e "\n                ---------------------------------------------------------------------------------"


samp='SRR742'
num=4152

nb=30
i=0

while [[ ! $num -eq 4826 ]]; do
    sample="$samp$num"
    echo -e " Retrieving $sample ..."
    sbatch --wrap=" fasterq-dump $sample -O ./MELANOMA/sc/"
    num=$((num + 1)) 
    i=$((i + 1)) 

    if [[ $i -eq $nb ]]; then
        echo -e "*********************************************************"
        while [[ ! -f ./MELANOMA/sc/SRR7424152.fastq ]]; do
            sleep 1; done

        while [[ ! $(ls ./MELANOMA/sc/* | wc -l) -eq $nb  ]]; do
            sleep 1; done

        echo " gzip files ...."
        gzip ./MELANOMA/sc/*.fastq
        echo " done."

        rm slurm*

        nb=$((nb + 30)) 
        echo $nb
    fi
done

while [[ ! $(ls ./MELANOMA/sc/* | wc -l) -eq 674 ]]; do
    sleep 1; done

gzip ./MELANOMA/sc/*.fastq
rm slurm*


echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                       Retrieving Bulk RNA-seq files                             "
echo -e "\n                ---------------------------------------------------------------------------------"

num=5010
while [[ ! $num -eq 5017 ]]; do
    sample="$samp$num"
    echo -e " Retrieving $sample ..."
    sbatch --wrap=" fasterq-dump $sample -O ./MELANOMA/bulk/ "
    num=$((num + 1)) 
done

while [[ ! -f ./MELANOMA/bulk/SRR7425010.fastq ]]; do
    sleep 1; done

while [[ ! $(ls ./MELANOMA/bulk/* | wc -l) -eq 7 ]]; do
    sleep 1; done

gzip ./MELANOMA/bulk/*.fastq
rm slurm*

echo -e "\n                ---------------------                Done !                 ---------------------\n"
end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"