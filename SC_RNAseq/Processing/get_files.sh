#   name : get_files.sh
#
#   Author (2023)  Safiya ATIA
#
#   This code is to be runned on your genotoul space
#   before processing your data

## Create & organize workspace
mkdir ./Reference
mkdir ./Annotation

## Download reference genome
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
mv ./Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz ./Reference
gunzip ./Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz


## Download reference genome
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
mv ./Homo_sapiens.GRCh38.109.gtf.gz ./Annotation
gunzip ./Annotation/Homo_sapiens.GRCh38.109.gtf.gz