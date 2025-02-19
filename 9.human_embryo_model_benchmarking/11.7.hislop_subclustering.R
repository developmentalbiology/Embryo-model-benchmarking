library(readxl)
library(devtools)
library(RColorBrewer) 
library(reshape2)
library(plotly)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(tidyverse)

setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250119_embryo_model_benchmarking_metrics")
outputDir = getwd()

# load data
obj_new <- readRDS("human_embryo_models_filtered_v2.rds")

################################subset hislop cluster_10 to check whether has cardiac mesoderm######################
Hislop_cl10 <- subset(obj_new[["Hislop"]], idents = "10" )

Hislop_cl10 <- subset(Hislop_cl10, subset = human_ref_transferred_final_anno_uncert < threshold)

Hislop_cl10 <- FindVariableFeatures(Hislop_cl10, selection.method = "vst", nfeatures = 2000)
Hislop_cl10 <- ScaleData(Hislop_cl10)
Hislop_cl10 <- RunPCA(Hislop_cl10, features = VariableFeatures(object = Hislop_cl10))
Hislop_cl10 <- FindNeighbors(Hislop_cl10, dims = 1:10) 
Hislop_cl10 <- FindClusters(Hislop_cl10, resolution = 0.5)
Hislop_cl10 <- RunUMAP(Hislop_cl10, dims = 1:10)

p <- DimPlot(Hislop_cl10, reduction = "umap")
ggsave("Hislop_cl10_umap.pdf", plot = p, width = 4, height = 4)


genes_to_plot2 = c("COL6A2", "PDGFRA",
                   "MYL7", "TNNI1", "TNNT2","NKX2-5",   # cardiomyocyte
                   "MEOX2", "MEOX1", "ALDH1A2", "FOXC2"    # somite
)


p<- FeaturePlot(Hislop_cl10, features = genes_to_plot2, reduction = "umap")
ggsave("Hislop_cl10_features.pdf", plot = p, width = 12, height = 7)


p<- DimPlot(Hislop_cl10, reduction = "umap", group.by = "human_ref_transferred_final_anno_unfiltered")
ggsave("Hislop_cl10_umap_human_ref_transferred.pdf", plot = p, width = 7, height = 4)



