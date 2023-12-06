#!/bin/bash

#   name : run_pipeline.sh
#
#   Author (2023)  Safiya ATIA
#
#   This code is to be runned on your genotoul space
#   to process your single cell data
#
#   (Genobioinfo cluster)

module load devel/java/17.0.6
module load bioinfo/Nextflow/22.12.0-edge

module load containers/singularity/3.9.9


sbatch --output=sc_pipeline.out -c 8 --wrap="nextflow run nf-core/scrnaseq -profile genotoul \
        --input SC_Glioma.csv \
        --outdir SC_results/ \
        --fasta ./Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        --gtf ./Annotation/Homo_sapiens.GRCh38.109.gtf \
        --protocol 10XV3 --aligner cellranger -resume"