
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readr)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250205_querymap_embryomodel_batches")
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
  DimPlot(refer, group.by = c("final_anno", "orig.ident"))
  
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
  predictions <- TransferData(anchorset = anchors, refdata = refer$final_anno, dims = 1:30)
  query <- AddMetaData(query, metadata = predictions)
  
  # Check if 'X' and 'predicted.id' are in the metadata before extracting
  meta <- query@meta.data
  meta$X <- rownames(meta)
  meta <- meta[, c("X", "predicted.id")]
  colnames(meta) <- c("X", "querymap_human_ref_final_anno")
  
  # Transfer lineage
  refer$final_lineage <- factor(refer$final_lineage)
  predictions <- TransferData(anchorset = anchors, refdata = refer$final_lineage, dims = 1:30)
  query <- AddMetaData(query, metadata = predictions)
  
  # Export metadata for lineage
  meta2 <- query@meta.data
  meta2$X <- rownames(meta2)
  meta2 <- meta2[, c("X", "predicted.id")]
  colnames(meta2) <- c("X", "querymap_human_ref_final_lineage")
  
  # Merge celltype and lineage metadata
  meta <- merge(meta, meta2, by = "X", all = TRUE)
  
  # Save the data as a CSV file
  write.csv(meta, file = output_csv, row.names = FALSE)
  print(paste("CSV saved:", output_csv))
}

# Define reference path (replace with actual path)
ref_path <- "D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250107_reference_comparasion/human_20250108.Rds"

# Load Seurat objects from query path (replace with actual path)
obj_new <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/R runing/human_embryo_models_filtered_v2.rds")

# Loop through each Seurat object and process it
for (obj_name in names(obj_new)) {
  query_object <- obj_new[[obj_name]]
  output_csv <- paste0(obj_name, "_human_ref_querymap.csv")
  
  # Run the function for each query object
  process_and_export_querymap(ref_path, query_object, output_csv)
}
