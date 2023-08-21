#!/bin/bash

##  sudo apt install sshpass

read -p " Genotoul username: " username
read -s -p " Genotoul Password: " password
echo -e "\n connexion to genotoul ..."
path="/work/${username}/SC_cellranger/cellranger/count"
print='yes'

mkdir -p samples/

for sample in $(sshpass -p ${password} ssh ${username}@genologin.toulouse.inra.fr ls ${path} | grep -v .yml$ | sed 's/sample-//g')
do
    if [ ${print} = 'yes' ];then
        echo -e "\n ---------- connexion to genotoul established ----------\n"
        print='no'
    fi
    mkdir -p samples/${sample}/
    echo -e " retrieving ${sample} ..."
    sshpass -p ${password} scp -r ${username}@genologin.toulouse.inra.fr:${path}/sample-${sample}/outs/raw_feature_bc_matrix/* ./samples/${sample}/
    echo -e " ${sample} done\n"
done
if [ ${print} = 'no' ];then
    echo -e " ---------- All samples were retrieved in the samples directory ! ----------\n"
fi