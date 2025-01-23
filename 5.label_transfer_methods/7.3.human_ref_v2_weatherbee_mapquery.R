
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readr)



##read rds and hvg
refer <- readRDS("human_20250108.Rds")
query <- readRDS("processed_weatherbee.rds")

ref_hvg <- read.csv("highly_variable_genes2000_human_ref.csv", header = TRUE)

query_hvg <- read.csv("highly_variable_genes2000_weatherbee.csv", header = TRUE)



####ref
refer <- NormalizeData(refer)
refer <- FindVariableFeatures(refer)
refer <- ScaleData(refer)
refer <- RunPCA(refer)
refer <- FindNeighbors(refer, dims = 1:30)
refer <- FindClusters(refer)
refer <- RunUMAP(refer, dims = 1:30,return.model = TRUE)
DimPlot(refer, group.by = c("final_anno", "orig.ident"))

####query
query <- NormalizeData(query)
query <- FindVariableFeatures(query)
query <- ScaleData(query)
query <- RunPCA(query)
query <- FindNeighbors(query, dims = 1:30)
query <- FindClusters(query)
query <- RunUMAP(query, dims = 1:30)
#DimPlot(query, group.by = c("reanno", "orig.ident"))

##transfer celltype
anchors <- FindTransferAnchors(reference = refer, query = query, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = refer$final_anno, dims = 1:30)
query <- AddMetaData(query, metadata = predictions)

######export metadata
meta <- query@meta.data
meta <- meta[,c("X","predicted.id")]
colnames(meta) <- c("X","querymap_human_ref_final_anno")


##transfer lineage
refer$final_lineage <- factor(refer$final_lineage)
predictions <- TransferData(anchorset = anchors, refdata = refer$final_lineage, dims = 1:30)
query <- AddMetaData(query, metadata = predictions)

######export metadata
meta2 <- query@meta.data
meta2 <- meta2[,c("X","predicted.id")]
colnames(meta2) <- c("X","querymap_human_ref_final_lineage")

meta <- merge(meta, meta2, by="X", all=TRUE)


# Save the data as a CSV file
write.csv(meta, file = "human_ref_weatherbee_querymap.csv", row.names = FALSE)


meta2 <- meta[match( rownames(query@meta.data),meta$X),]

query$final_anno <- meta2$querymap_human_ref_final_anno
query$final_lineage <- meta2$querymap_human_ref_final_lineage

DimPlot(query, group.by = "final_anno")


