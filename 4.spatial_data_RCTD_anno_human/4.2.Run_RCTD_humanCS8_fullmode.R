
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

setwd("/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250724_human_refined_anno_RCTD_v3")
outputDir = getwd()

# load data
cs8 <- readRDS("/storage2/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/R_runing/cs8_human_embryo.rds")
ref <- readRDS("human_clustering_20250724_v3.rds")

ref$'reanno' <- as.character(ref$'reanno')

ref <- UpdateSeuratObject(ref)
Idents(ref) <- "reanno"
#Idents(ref) <- "lineage"

#debug, remove the cluster which only has one cell 
cell_counts <- table(Idents(ref))
valid_clusters <- names(cell_counts[cell_counts > 1])

ref <- subset(ref, idents = valid_clusters)

#Rename the Problematic Cell Types
#ref$leiden_3 <- gsub("/", "_", ref$leiden_3)

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$reanno)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)


# set up query with the RCTD function SpatialRNA
counts <- cs8@assays[["RNA"]]@counts
#get xy coordinates
coords <- cs8@meta.data[,c("x","y")]
coords <- as.data.frame(coords)
query <- SpatialRNA(coords, counts, colSums(counts))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 20, CELL_MIN_INSTANCE = 3, UMI_min = 0, counts_MIN = 0, UMI_min_sigma = 0)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

saveRDS(RCTD, file = "RCTD_human_cs8_leiden_3_20250724_v3.rds")

