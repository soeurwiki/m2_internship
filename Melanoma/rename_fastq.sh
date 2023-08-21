start=$(date +%s)

echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                          Renaming SC RNA-seq files                              "
echo -e "\n                ---------------------------------------------------------------------------------"


samp='SRR742'
num=4152

echo -e "\n                -----------------               Before treatment                -----------------"

i=0
while [[ ! $num -eq 4324 ]]; do
    sample="$samp$num"
    if [[ -f ./MELANOMA/sc/$sample.fastq.gz ]]; then
        echo -e " Renaming $sample ..."
        mv ./MELANOMA/sc/$sample.fastq.gz ./MELANOMA/sc/bef_$sample.fastq.gz
        i=$((i + 1))  
    else
        echo -e " !!!!! $sample doesn't exist !!!!!"
    fi
    num=$((num + 1)) 
done   
echo "--------  $i samples  --------"

echo -e "\n                -----------------                4d on treatment                 -----------------"

i=0
while [[ ! $num -eq 4479 ]]; do
    sample="$samp$num"
    if [[ -f ./MELANOMA/sc/$sample.fastq.gz ]]; then
        echo -e " Renaming $sample ..."
        mv ./MELANOMA/sc/$sample.fastq.gz ./MELANOMA/sc/4d_$sample.fastq.gz
        i=$((i + 1))  
    else
        echo -e " !!!!! $sample doesn't exist !!!!!"
    fi
    num=$((num + 1)) 
done
echo "--------  $i samples  --------"

echo -e "\n                -----------------                28d on treatment                 -----------------"

i=0
while [[ ! $num -eq 4678 ]]; do
    sample="$samp$num"
    if [[ -f ./MELANOMA/sc/$sample.fastq.gz ]]; then
        echo -e " Renaming $sample ..."
        mv ./MELANOMA/sc/$sample.fastq.gz ./MELANOMA/sc/28d_$sample.fastq.gz
        i=$((i + 1))  
    else
        echo -e " !!!!! $sample doesn't exist !!!!!"
    fi
    num=$((num + 1))
done
echo "--------  $i samples  --------"

echo -e "\n                -----------------                57d on treatment                 -----------------"

i=0
while [[ ! $num -eq 4826 ]]; do
    sample="$samp$num"
    if [[ -f ./MELANOMA/sc/$sample.fastq.gz ]]; then
        echo -e " Renaming $sample ..."
        mv ./MELANOMA/sc/$sample.fastq.gz ./MELANOMA/sc/57d_$sample.fastq.gz
        i=$((i + 1))  
    else
        echo -e " !!!!! $sample doesn't exist !!!!!"
    fi
    num=$((num + 1)) 
done
echo "--------  $i samples  --------"
 
echo -e "\n                ---------------------------------------------------------------------------------"
echo -e "\n                                         Renaming bulk RNA-seq files                             "
echo -e "\n                ---------------------------------------------------------------------------------"

i=0
num=5010
while [[ ! $num -eq 5017 ]]; do
    sample="$samp$num"
    if [[ -f ./MELANOMA/bulk/$sample.fastq.gz ]]; then
        echo -e " Renaming $sample ..."

            ## 28d on treatment ##
        if [[ $num -eq 5012 || $num -eq 5013]]; then
            mv ./MELANOMA/bulk/$sample.fastq.gz ./MELANOMA/bulk/28d_$sample.fastq.gz
        else
            ## 10d on treatment ##
            mv ./MELANOMA/bulk/$sample.fastq.gz ./MELANOMA/bulk/10d_$sample.fastq.gz
        fi
        i=$((i + 1))  
    else
        echo -e " !!!!! $sample doesn't exist !!!!!"
    fi
    num=$((num + 1)) 
done   
echo "--------  $i samples  --------"

echo -e "\n                --------------------                Done !                 --------------------\n"

end=$(date +%s)
echo "Elapsed Time in total : $(($end-$start)) seconds"
