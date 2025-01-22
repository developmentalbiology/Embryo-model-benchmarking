
.libPaths(c("/home/liuxiaodongLab/fanxueying/R_LIBS",
            "/home/liuxiaodongLab/liuyifang/miniconda3/envs/R-4.4.0/lib/R/library"
))


library(Seurat)
library(tidyverse)
library(spacexr)
library(Matrix)
library(doParallel)
library(parallel)
library(quadprog)

setwd("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/RCTD")
outputDir = getwd()

# load data
cs8 <- readRDS("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/R runing/cs8_human_embryo.rds")
ref <- readRDS("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/2025Jan/human_20250108.Rds")


ref <- UpdateSeuratObject(ref)
Idents(ref) <- "final_anno"
#Idents(ref) <- "lineage"

#debug, remove the cluster which only has one cell 
cell_counts <- table(Idents(ref))
valid_clusters <- names(cell_counts[cell_counts > 1])

ref <- subset(ref, idents = valid_clusters)


class_df <- ref@meta.data[,c("final_lineage","final_anno")]
class_df <- unique(class_df)
class_df <- data.frame(class_df$final_lineage, row.names = class_df$final_anno)
colnames(class_df)[1] = "class"
print(class_df)


#Rename the Problematic Cell Types
ref$final_anno <- gsub("/", "_", ref$final_anno)

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$final_anno)
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

RCTD <- create.RCTD(query, reference, max_cores = 1, CELL_MIN_INSTANCE = 3, UMI_min = 0, counts_MIN = 0, UMI_min_sigma = 0)
#RCTD <- create.RCTD(query, reference,class_df = class_df, max_cores = 8, CELL_MIN_INSTANCE = 3, UMI_min = 0, counts_MIN = 0, UMI_min_sigma = 0)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

saveRDS(RCTD, file = "RCTD_cs8_human_final_anno_fullmode_celltype_20250108.rds")

