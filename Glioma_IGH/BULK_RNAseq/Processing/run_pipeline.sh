module load system/jdk-12.0.2
module load bioinfo/Nextflow-v22.12.0-edge
module load system/singularity-3.7.3
sh /usr/local/bioinfo/src/NextflowWorkflows/create_nfx_dirs.sh

sbatch --output=bulk_pipeline.out -c 8 --wrap=" nextflow run nf-core/rnaseq -profile genotoul \
        --input Bulk_Glioma.csv \
        --outdir ./BULK_nextflow/ \
        --fasta ./Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        --gtf ./Annotation/Homo_sapiens.GRCh38.109.gtf -resume"