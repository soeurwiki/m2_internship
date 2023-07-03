module load system/R-4.2.3_Miniconda3
module load system/singularity-3.7.3

sample=$(ls /work/project/LeCam_U1194/AIRELLE/RNASEQ/nextflow_full/star_salmon/featurecounts | grep .featureCounts.txt$ | sed 's/.featureCounts.txt//g')

python count_mtx.py  ${sample}
