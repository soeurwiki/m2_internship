import pandas as pd
import sys

args = sys.argv
print(len(args))
print(args[0])
print(args[82])

samples = {}
for i in range(1,len(args)) :
        path = '/work/project/LeCam_U1194/AIRELLE/RNASEQ/nextflow_full/star_salmon/'+ args[i] +'/quant.genes.sf'
        samples[args[i]] = pd.read_csv(path,sep='\t', usecols = ["Name", "NumReads"], index_col=0)
        samples[args[i]] = samples[args[i]].rename(columns={'NumReads': args[i]})

merged_df = samples[args[1]]
print('---------------------------------\n')
for i in samples :
    if(i != args[1]) : ## sample different from first one
        merged_df = pd.merge(merged_df, samples[i], on=['Name'], how='outer')
    else :
        print("out")

merged_df.to_csv("mtx_RL.tsv",sep='\t',header=True, index=True)