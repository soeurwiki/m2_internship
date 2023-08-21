library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)

library(ggplot2)
library("scales")

seurat_obj <- readRDS(file = "./rds_obj/Norm_merged_seurat.rds")

pdf("metabo_1.pdf")

as.tibble(
  seurat_obj@assays$SCT@data[,1:100]
) %>%
  pivot_longer(
    cols=everything(),
    names_to="cell",
    values_to="expression"
  ) %>%
  ggplot(aes(x=expression, group=cell)) +
  geom_density() +
  coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))



genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)
head(genes)


features <- rownames(seurat_obj@assays$SCT@counts)
length(features)

# keep genes names seeen un genes$Gene_name
matched_features <- features[features %in% genes$Gene_name]
length(matched_features)

counts <- GetAssayData(seurat_obj, assay = "SCT")
length(rownames(counts) )

counts <- counts[(which(rownames(counts) %in% matched_features)),]
metabo_obj <- subset(seurat_obj, features = rownames(counts))
metabo_obj


apply(seurat_obj@assays$SCT@data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression



FeatureScatter(metabo_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by ="orig.ident")
FeatureScatter(metabo_obj,feature1 = "nCount_RNA", feature2 = "percent.largest_gene",group.by ="orig.ident")

FeatureScatter(metabo_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by ="orig.ident") + geom_smooth(method = 'lm')

metadata <- metabo_obj@meta.data

VlnPlot(metabo_obj,features="nCount_RNA", group.by = "orig.ident", pt.size = 0)
RidgePlot(metabo_obj, features = "nCount_RNA", group.by = "orig.ident")

metadata %>% 
  ggplot(aes( x=nCount_RNA,y=orig.ident, fill= orig.ident)) + 
	geom_boxplot() + 
	theme_classic()+
  scale_x_log10()

nCount_metabo <- as.data.frame(rowSums(metabo_obj[["nCount_RNA"]]))
nCount_metabo$Genes <- "metabo"
colnames(nCount_metabo)[1] <- "nCount"

nCount_all <- as.data.frame(rowSums(seurat_obj[["nCount_RNA"]]))
nCount_all$Genes <- "all"
colnames(nCount_all)[1] <- "nCount"

other <- nCount_all
other$Genes <- "others"
dim(other)
other <- other[!(rownames(other) %in% rownames(metabo_obj$RNA)),]

expr <- rbind(nCount_all, nCount_metabo, other)
head(expr)
ggplot(expr, aes(x=Genes,y=nCount)) + geom_boxplot()

VlnPlot(metabo_obj,features="nFeature_RNA", group.by = "orig.ident", pt.size = 0)
RidgePlot(metabo_obj, features = "nFeature_RNA", group.by = "orig.ident")

metadata %>% 
  ggplot(aes( x=nFeature_RNA,y=orig.ident, fill= orig.ident)) + 
	geom_boxplot() + 
	theme_classic()+
    scale_x_log10()
    scale_x_log10()
    
  scale_x_log10()

ggplot(metabo_obj[[]]) + 
   geom_histogram(aes(nFeature_RNA), 
                  color = "#558bdc", fill= "#173664",
                  binwidth = 50) + 
  ggtitle("Distribution of nFeature_RNA ") + NoLegend()


expression_metabo <- as.data.frame(rowSums(metabo_obj@assays$SCT@data))
expression_metabo$Genes <- "metabo"
colnames(expression_metabo)[1] <- "expression"

expression_all <- as.data.frame(rowSums(seurat_obj@assays$SCT@data))
expression_all$Genes <- "all"
colnames(expression_all)[1] <- "expression"

other <- expression_all
other$Genes <- "others"
dim(other)
other <- other[!(rownames(other) %in% rownames(metabo_obj$RNA)),]

expr <- rbind(expression_all, expression_metabo, other)
head(expr)
ggplot(expr, aes(x=Genes,y=log(expression+1))) + geom_boxplot()


as.tibble(
  metabo_obj@assays$SCT@data[,1:100]
) %>%
  pivot_longer(
    cols=everything(),
    names_to="cell",
    values_to="expression"
  ) %>%
  ggplot(aes(x=expression, group=cell)) +
  geom_density() +
  coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))

expression_metabo <- as.data.frame(rowSums(metabo_obj@assays$SCT@data))
expression_metabo$Genes <- "metabo"
colnames(expression_metabo)[1] <- "expression"

expression_all <- as.data.frame(rowSums(seurat_obj@assays$SCT@data))
expression_all$Genes <- "all"
colnames(expression_all)[1] <- "expression"

other <- expression_all
other$Genes <- "others"
dim(other)
other <- other[!(rownames(other) %in% rownames(metabo_obj$RNA)),]

expr <- rbind(expression_all, expression_metabo, other)
head(expr)
ggplot(expr, aes(x=Genes,y=log(expression+1))) + geom_boxplot()


VlnPlot(metabo_obj, features = "percent.mt", group.by = "orig.ident", pt.size = 0)
RidgePlot(metabo_obj, features = "percent.mt", group.by = "orig.ident")

ggplot(metabo_obj[[]]) + 
   geom_histogram(aes(percent.mt), 
                  color = "#558bdc", fill= "#173664",
                  binwidth = 0.5) + 
  ggtitle("Distribution of Percentage Mitochondrion") + NoLegend()

ggplot(metabo_obj[[]]) + 
   geom_histogram(aes(percent.ribosomal), 
                  color = "#558bdc", fill= "#173664",
                  binwidth = 0.7) + 
  ggtitle("Distribution of Percentage Ribosomal") + NoLegend()

VlnPlot(metabo_obj, features = "percent.largest_gene", group.by = "orig.ident", pt.size = 0)
RidgePlot(metabo_obj, features = "percent.largest_gene", group.by = "orig.ident")

ggplot(metabo_obj[[]]) + 
   geom_histogram(aes(percent.largest_gene), 
                  color = "#558bdc", fill= "#173664",
                  binwidth = 0.7) + 
  ggtitle("Distribution of Percentage Largest Gene") + NoLegend()

metabo_obj@meta.data %>%
  group_by(orig.ident,Phase) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=orig.ident,y=percent, fill=Phase)) +
  scale_fill_manual( values=c(G1='azure4',G2M='dodgerblue',S='tomato3' )) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per sample")

as_tibble(metabo_obj[[]]) %>%
  ggplot(aes(Phase, fill=Phase)) +
  scale_fill_manual( values=c(G1='azure4',G2M='dodgerblue',S='tomato3' ))+ geom_bar()


as_tibble(metabo_obj[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  scale_color_manual( values=c(G1='azure4',G2M='dodgerblue',S='tomato3' )) +
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))

# identify the 10 most highly variable genes
top10 <- head(VariableFeatures(metabo_obj), 10)

# plot the variable features 
plot <- VariableFeaturePlot(metabo_obj)
LabelPoints(plot = plot, points = head(top10,10), repel = TRUE)


metabo_obj <- RunPCA(metabo_obj, features = VariableFeatures(object = metabo_obj))
VizDimLoadings(metabo_obj, dims = 1:2, reduction = "pca")

mat <- Seurat::GetAssayData(metabo_obj, assay = "SCT", slot = "scale.data")
pca <- metabo_obj[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues

varExplained = eigValues / total_variance

pc1 = percent(varExplained[1], accuracy = 0.01)
pc2 = percent(varExplained[2], accuracy = 0.01)

DimPlot(metabo_obj, reduction = "pca", group.by = "orig.ident") + xlab(paste0("PC_1 : ", pc1)) + ylab(paste0("PC_2 :", pc2 ))
DimPlot(metabo_obj, reduction = "pca", group.by = "Phase", cols = c(G1='azure4',G2M='dodgerblue',S='tomato3'))
DimPlot(metabo_obj, reduction = "pca", group.by = "largest_gene", label = TRUE, label.size = 3) + NoLegend()

saveRDS(metabo_obj, file ="merged_seurat_metabo.rds")

ElbowPlot(metabo_obj, ndims = 40)

dev.off()
print("done")
