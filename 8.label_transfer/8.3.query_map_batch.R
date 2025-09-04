
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readr)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_querymap_batch_v3")
outputDir = getwd()

# Define function
process_and_export_querymap <- function(ref_path, query_object, output_csv) {
  # Read RDS file for reference
  refer <- readRDS(ref_path)
  
  # Process the query dataset (Seurat object)
  query <- query_object
  
  # Process reference dataset
  refer <- NormalizeData(refer)
  refer <- FindVariableFeatures(refer)
  refer <- ScaleData(refer)
  refer <- RunPCA(refer)
  refer <- FindNeighbors(refer, dims = 1:30)
  refer <- FindClusters(refer)
  refer <- RunUMAP(refer, dims = 1:30, return.model = TRUE)
  DimPlot(refer, group.by = c("reanno", "orig.ident"))
  
  # Process query dataset
  query <- NormalizeData(query)
  query <- FindVariableFeatures(query)
  query <- ScaleData(query)
  query <- RunPCA(query)
  query <- FindNeighbors(query, dims = 1:30)
  query <- FindClusters(query)
  query <- RunUMAP(query, dims = 1:30)
  
  # Transfer celltype
  anchors <- FindTransferAnchors(reference = refer, query = query, dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = refer$reanno, dims = 1:30)
  query <- AddMetaData(query, metadata = predictions)
  
  # Check if 'X' and 'predicted.id' are in the metadata before extracting
  meta <- query@meta.data
  meta$X <- rownames(meta)
  meta <- meta[, c("X", "predicted.id")]
  colnames(meta) <- c("X", "querymap_human_ref_reanno")
  
  # Transfer lineage
  refer$lineage <- factor(refer$lineage)
  predictions <- TransferData(anchorset = anchors, refdata = refer$lineage, dims = 1:30)
  query <- AddMetaData(query, metadata = predictions)
  
  # Export metadata for lineage
  meta2 <- query@meta.data
  meta2$X <- rownames(meta2)
  meta2 <- meta2[, c("X", "predicted.id")]
  colnames(meta2) <- c("X", "querymap_human_ref_lineage")
  
  # Merge celltype and lineage metadata
  meta <- merge(meta, meta2, by = "X", all = TRUE)
  
  # Save the data as a CSV file
  write.csv(meta, file = output_csv, row.names = FALSE)
  print(paste("CSV saved:", output_csv))
}

# Define reference path (replace with actual path)
ref_path <- "D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/data/human_clustering_20250724_v3.rds"

# Load Seurat objects from query path (replace with actual path)
obj_new <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/R runing/human_embryo_models_filtered_v2.rds")

# Loop through each Seurat object and process it
for (obj_name in names(obj_new)) {
  query_object <- obj_new[[obj_name]]
  output_csv <- paste0(obj_name, "_human_ref_querymap.csv")
  
  # Run the function for each query object
  process_and_export_querymap(ref_path, query_object, output_csv)
}

