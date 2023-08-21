# Seurat object  
## Load required packages  

library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library("scales")


## Load seurat object (RAW)  

seurat_obj <- readRDS(file = "./rds_obj/Raw_merged_seurat.rds")

print(seurat_obj)
table(seurat_obj$orig.ident)

metadata <- seurat_obj@meta.data

pdf("raw.pdf")
VlnPlot(seurat_obj,features="nCount_RNA", group.by = "orig.ident", pt.size = 0)
RidgePlot(seurat_obj, features = "nCount_RNA", group.by = "orig.ident")

metadata %>% 
  ggplot(aes( x=nCount_RNA,y=orig.ident, fill= orig.ident)) + 
	geom_boxplot() + 
	theme_classic()+
  scale_x_log10()

VlnPlot(seurat_obj,features="nFeature_RNA", group.by = "orig.ident", pt.size = 0)
RidgePlot(seurat_obj, features = "nFeature_RNA", group.by = "orig.ident")

metadata %>% 
  ggplot(aes( x=nFeature_RNA,y=orig.ident, fill= orig.ident)) + 
	geom_boxplot() + 
	theme_classic()+
  scale_x_log10()

ggplot(seurat_obj[[]]) + 
   geom_histogram(aes(nFeature_RNA, 
                  color = "#173664"),
                  binwidth = 50) + 
  ggtitle("Distribution of nFeature_RNA ") + NoLegend()



apply(seurat_obj@assays$RNA@data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression

as.tibble(
  seurat_obj@assays$RNA@data[,1:100]
) %>%
  pivot_longer(
    cols=everything(),
    names_to="cell",
    values_to="expression"
  ) %>%
  ggplot(aes(x=expression, group=cell)) +
  geom_density() +
  coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))


seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj,features="percent.mt", group.by = "orig.ident", pt.size = 0)
RidgePlot(seurat_obj, features = "percent.mt", group.by = "orig.ident")


ggplot(seurat_obj[[]]) + 
   geom_histogram(aes(percent.mt, 
                  color = "#173664"),
                  binwidth = 0.5) + 
  ggtitle("Distribution of Percentage Mitochondrion") + NoLegend()


PercentageFeatureSet(seurat_obj,pattern="^RP[LS]") -> seurat_obj[["percent.ribosomal"]] 
VlnPlot(seurat_obj,features="percent.ribosomal", group.by = "orig.ident", pt.size = 0)
RidgePlot(seurat_obj, features = "percent.ribosomal", group.by = "orig.ident")


ggplot(seurat_obj[[]]) + 
   geom_histogram(aes(percent.ribosomal, 
                  color = "#173664"),
                  binwidth = 0.5) + 
  ggtitle("Distribution of Percentage Ribosomal") + NoLegend()

apply(
  seurat_obj@assays$RNA@counts,
  2,
  max
) -> largest_count

apply(
  seurat_obj@assays$RNA@counts,
  2,
  which.max
) -> largest_index

rownames(seurat_obj)[largest_index] -> seurat_obj$largest_gene

100 * largest_count / seurat_obj$nCount_RNA -> seurat_obj$percent.largest_gene

VlnPlot(seurat_obj,features="percent.largest_gene", group.by = "orig.ident", pt.size = 0)
RidgePlot(seurat_obj, features = "percent.largest_gene", group.by = "orig.ident")



ggplot(seurat_obj[[]]) + 
   geom_histogram(aes(percent.largest_gene, 
                  color = "#173664"),
                  binwidth = 0.7) + 
  ggtitle("Distribution of Percentage Largest Gene") + NoLegend()



FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by ="orig.ident")
FeatureScatter(seurat_obj,feature1 = "nCount_RNA", feature2 = "percent.largest_gene",group.by ="orig.ident")

FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by ="orig.ident") + geom_smooth(method = 'lm')

seurat_obj <- CellCycleScoring(seurat_obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)

seurat_obj@meta.data %>%
  group_by(orig.ident,Phase) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=orig.ident,y=percent, fill=Phase)) +
  scale_fill_manual( values=c(G1='azure4',G2M='dodgerblue',S='tomato3' )) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per sample")

as_tibble(seurat_obj[[]]) %>%
  ggplot(aes(Phase, fill=Phase)) +
  scale_fill_manual( values=c(G1='azure4',G2M='dodgerblue',S='tomato3' ))+ geom_bar()


as_tibble(seurat_obj[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  scale_color_manual( values=c(G1='azure4',G2M='dodgerblue',S='tomato3' )) +
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))

dev.off()
print("done")