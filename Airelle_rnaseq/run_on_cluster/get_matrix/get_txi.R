library(tximport)

args = commandArgs(trailingOnly=TRUE) # all samples name
files <- file.path("RL_data", args, "quant.sf")
names(files) <-  args

# tx2gene : must be dataframe with col1 = Id, col2 = names
tx2gene <- read.csv("RL_data/salmon_tx2gene.tsv", sep="\t")
tx2gene  <- tx2gene [, c(1, 2)]
txi <- tximport(files, type="salmon", tx2gene = tx2gene,countsFromAbundance="lengthScaledTPM")
saveRDS(txi, "txi.rds")