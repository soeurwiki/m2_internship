# module load containers/singularity/3.9.9 
# sbatch -c 8 --mem-per-cpu=2048 --wrap=" singularity exec r.sif Rscript markers_to_csv.R"

library(tidyverse)

## all_genes
markers <- readRDS(file = "./results/html/all_samples_cluster.markers.rds")

markers_0 <- rownames(markers[which(markers$cluster==0),])
markers_1 <- rownames(markers[which(markers$cluster==1),])
markers_2 <- rownames(markers[which(markers$cluster==2),])
markers_3 <- rownames(markers[which(markers$cluster==3),])
markers_4 <- rownames(markers[which(markers$cluster==4),])
markers_5 <- rownames(markers[which(markers$cluster==5),])
markers_6 <- rownames(markers[which(markers$cluster==6),])
markers_7 <- rownames(markers[which(markers$cluster==7),])
markers_8 <- rownames(markers[which(markers$cluster==8),])
markers_9 <- rownames(markers[which(markers$cluster==9),])
markers_10 <- rownames(markers[which(markers$cluster==10),])
markers_11 <- rownames(markers[which(markers$cluster==11),])
markers_12 <- rownames(markers[which(markers$cluster==12),])
markers_13 <- rownames(markers[which(markers$cluster==13),])
markers_14 <- rownames(markers[which(markers$cluster==14),])
markers_15 <- rownames(markers[which(markers$cluster==15),])
markers_16 <- rownames(markers[which(markers$cluster==16),])
markers_17 <- rownames(markers[which(markers$cluster==17),])
rm(markers)

max_ln <- max(c(length(markers_0), length(markers_1), length(markers_0),length(markers_3),
                length(markers_4),length(markers_5),length(markers_6),length(markers_7),
                length(markers_8),length(markers_9),length(markers_10),length(markers_11),
                length(markers_12),length(markers_13),length(markers_14),length(markers_15),
                length(markers_16),length(markers_17)))

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
                 cl10 = c(markers_10, rep(NA, max_ln - length(markers_10))),
                 cl11 = c(markers_11, rep(NA, max_ln - length(markers_11))),
                 cl12 = c(markers_12, rep(NA, max_ln - length(markers_12))),
                 cl13 = c(markers_13, rep(NA, max_ln - length(markers_13))),
                 cl14 = c(markers_14, rep(NA, max_ln - length(markers_14))),
                 cl15 = c(markers_15, rep(NA, max_ln - length(markers_15))),
                 cl16 = c(markers_16, rep(NA, max_ln - length(markers_16))),
                 cl17 = c(markers_17, rep(NA, max_ln - length(markers_17)))
)

write.csv(dft, "./results/html/markers.csv", row.names=FALSE)

## metabo
markers <- readRDS(file = "./results/html/all_samples_cluster.markers_metabo.rds")

markers_0 <- rownames(markers[which(markers$cluster==0),])
markers_1 <- rownames(markers[which(markers$cluster==1),])
markers_2 <- rownames(markers[which(markers$cluster==2),])
markers_3 <- rownames(markers[which(markers$cluster==3),])
markers_4 <- rownames(markers[which(markers$cluster==4),])
markers_5 <- rownames(markers[which(markers$cluster==5),])
markers_6 <- rownames(markers[which(markers$cluster==6),])
markers_7 <- rownames(markers[which(markers$cluster==7),])
markers_8 <- rownames(markers[which(markers$cluster==8),])
markers_9 <- rownames(markers[which(markers$cluster==9),])
markers_10 <- rownames(markers[which(markers$cluster==10),])
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

write.csv(dft, "./results/html/markers_metabo.csv", row.names=FALSE)

print("Finished.")