# module load containers/singularity/3.9.9 
# sbatch -c 8 --mem-per-cpu=2048 --wrap=" singularity exec r.sif Rscript get_top_markers_metabo.R"

library(tidyverse)

t=25
markers <- readRDS(file = "./results/html/all_samples_cluster.markers_metabo.rds")


prolif <- markers[which( markers$diff_p_val == markers$max_pval), ]
diff <- markers[which( markers$prolif_p_val == markers$max_pval), ]

table(rownames(prolif) == rownames(diff))

print(dim(prolif))
print(dim(diff))


prolif <- prolif[, c('prolif_p_val', 'prolif_avg_log2FC', 'prolif_p_val_adj', 'gene', 'cluster' )]
diff <- diff[, c('diff_p_val', 'diff_avg_log2FC', 'diff_p_val_adj', 'gene', 'cluster' )]

## remove genes that are both in prolif & diff (because of same values)
diff <- diff[! rownames(diff) %in% rownames(prolif),]
table(rownames(prolif) == rownames(diff))

diff_markers <- diff[order(diff$diff_p_val_adj, decreasing = FALSE),]
prolif_markers <- prolif[order(prolif$prolif_p_val_adj, decreasing = FALSE),]

print(dim(prolif_markers))
print(dim(diff_markers))

## keep best markers
diff_markers %>%
    group_by(cluster) %>%
    top_n(n = t, wt = diff_avg_log2FC) -> top_diff

top_diff <- top_diff[, c('gene', 'cluster' )]

prolif_markers %>%
    group_by(cluster) %>%
    top_n(n = t, wt = prolif_avg_log2FC) -> top_prolif

top_prolif <- top_prolif[, c('gene', 'cluster' )]

# create dataframe of top
markers <- rbind(top_diff, top_prolif)
print(head(markers))

markers_0 <- markers[which(markers$cluster==0),]$gene
markers_1 <- markers[which(markers$cluster==1),]$gene
markers_2 <- markers[which(markers$cluster==2),]$gene
markers_3 <- markers[which(markers$cluster==3),]$gene
markers_4 <- markers[which(markers$cluster==4),]$gene
markers_5 <- markers[which(markers$cluster==5),]$gene
markers_6 <- markers[which(markers$cluster==6),]$gene
markers_7 <- markers[which(markers$cluster==7),]$gene
markers_8 <- markers[which(markers$cluster==8),]$gene
markers_9 <- markers[which(markers$cluster==9),]$gene
markers_10 <- markers[which(markers$cluster==10),]$gene
rm(markers)

max_ln <- max(c(length(markers_0), length(markers_1), length(markers_0),length(markers_3),
                length(markers_4),length(markers_5),length(markers_6),length(markers_7),
                length(markers_8),length(markers_9),length(markers_10)))

dft <- data.frame(cl0 = c(markers_0, rep(NA, max_ln - length(markers_0))),
                 cl1 = c(markers_1, rep(NA, max_ln - length(markers_1))),
                 cl2 = c(markers_2, rep(NA, max_ln - length(markers_2))),
                 cl3 = c(markers_3, rep(NA, max_ln - length(markers_3))),
                 cl4 = c(markers_4, rep(NA, max_ln - length(markers_4))),
                 cl5 = c(markers_5, rep(NA, max_ln - length(markers_5))),
                 cl6 = c(markers_6, rep(NA, max_ln - length(markers_6))),
                 cl7 = c(markers_7, rep(NA, max_ln - length(markers_7))),
                 cl8 = c(markers_8, rep(NA, max_ln - length(markers_8))),
                 cl9 = c(markers_9, rep(NA, max_ln - length(markers_9))),
                 cl10 = c(markers_10, rep(NA, max_ln - length(markers_10)))
)

write.csv(dft, "./results/html/top_markers_metabo.csv", row.names=FALSE)

print('Finished.')