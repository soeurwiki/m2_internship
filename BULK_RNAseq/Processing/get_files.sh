#   name : get_files.sh
#
#   Author (2023)  Safiya ATIA

## Create & organize workspace
mkdir ./Reference
mkdir ./Annotation

## Téléchargement du génome de référence
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
mv ./Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz ./Reference
gunzip ./Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz


## Téléchargement de l'annotation
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
mv ./Homo_sapiens.GRCh38.109.gtf.gz ./Annotation
gunzip ./Annotation/Homo_sapiens.GRCh38.109.gtf.gz