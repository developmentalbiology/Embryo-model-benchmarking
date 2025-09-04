
.libPaths(c("/home/liuxiaodongLab/fanxueying/R_LIBS",
             "/home/liuxiaodongLab/fanxueying/miniconda3/envs/r_env/lib/R/library",
             "/soft/devtools/R/R-4.4.3_installation/lib64/R/library"
))


library(Seurat)
library(tidyverse)
library(spacexr)
library(Matrix)
library(doParallel)
library(parallel)
library(quadprog)

setwd("/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250706_monkey_RCTD_reanno")
outputDir = getwd()

# load data
cs8 <- readRDS("/storage2/liuxiaodongLab/liuyifang/Projects/JiabinChen/Embryo_JiabinChen/2024-11-25_341_CS910-with-updated-Spateo-aligned-coordinate/data/sobj_CS910_spateo_341.rds")
ref <- readRDS("monkey_clustering.rds")

ref$reanno <- as.character(ref$reanno)

ref <- UpdateSeuratObject(ref)
Idents(ref) <- "reanno"

#debug, remove the cluster which only has one cell 
cell_counts <- table(Idents(ref))
valid_clusters <- names(cell_counts[cell_counts > 1])

ref <- subset(ref, idents = valid_clusters)

#Rename the Problematic Cell Types
ref$reanno <- gsub("/", "_", ref$reanno)

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$reanno)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)


# set up query with the RCTD function SpatialRNA
counts <- cs8@assays[["RNA"]]@layers$counts
rownames(counts) <- rownames(cs8)
colnames(counts) <- colnames(cs8)

#get xy coordinates
coords <- cs8@reductions$Spatial45@cell.embeddings
coords <- as.data.frame(coords)
colnames(coords) <- c("x","y")
query <- SpatialRNA(coords, counts, colSums(counts))

#run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 20, CELL_MIN_INSTANCE = 3, UMI_min = 0, counts_MIN = 0, UMI_min_sigma = 0)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

saveRDS(RCTD, file = "RCTD_cs910_monkey_reanno_20250706.rds")

