#!/bin/bash

module load system/jdk-12.0.2
module load bioinfo/Nextflow-v22.12.0-edge

module load system/singularity-3.7.3

sh /usr/local/bioinfo/src/NextflowWorkflows/create_nfx_dirs.sh

sbatch --output=sc_pipeline.out -c 8 --wrap="nextflow run nf-core/scrnaseq  --input Gliome_SC.csv --outdir SC_results/ \
--fasta ./Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
--gtf ./Annotation/Homo_sapiens.GRCh38.109.gtf --protocol 10XV3 \
--aligner cellranger  -profile genotoul -resume"