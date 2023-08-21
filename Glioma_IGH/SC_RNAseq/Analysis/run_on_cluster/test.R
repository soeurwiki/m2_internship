
# Seurat objects  
## Load required packages  
library(Seurat)
library(tidyverse)
library(ggplot2)
library("scales")


## Load seurat object (ALL GENES)

line <- "LGG275"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))

line <- "LGG336"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))

line <- "BT1"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))

line <- "LGG85"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))

line <- "BT2"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))

line <- "BT54"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))

line <- "LGG349"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))

line <- "BT88"

seurat_obj <- readRDS(file = paste0('./results/rds/samples/', line, '_Clusters.rds'))

print(seurat_obj)

genes <- read.csv("all_metabo.tsv", sep = '\t', header = TRUE)

counts <- GetAssayData(seurat_obj, assay = "integrated")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "RNA")
print(length(rownames(counts)))

counts <- GetAssayData(seurat_obj, assay = "SCT")
print(length(rownames(counts)))