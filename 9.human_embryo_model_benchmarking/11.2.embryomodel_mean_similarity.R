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
library(pheatmap)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250119_embryo_model_benchmarking_metrics")
outputDir = getwd()

# Define the lineages
lineages <- c("Endoderm","epi","ExE_endo","Exe_meso", "hemogenic","mesoderm", "neural_ecto", "non_neuro_ecto","Notochord",  "PGC",  "TE_TrB", "Gastru")

#load human reference data
human <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250107_reference_comparasion/human_20250108.Rds")

# load similarity tables
human_sim <- read.csv("correlation_human_ref_celltype.csv", row.names = 1)

# calculate cell_type mean within each lineage
data <- human_sim

correlation_human_ref_celltype <- do.call(cbind, lapply(lineages, function(lineage) {
  
  sub <- unique(human$final_anno[which(human$final_lineage == lineage)])
  name <- paste0("RNA.", sub)
  name <- gsub("_",".",name)
  name <- gsub(" ",".",name)
  
  # Ensure the rows exist in 'data'
  valid_rows <- name[name %in% rownames(data)]
  
  if (length(valid_rows) > 0) {
    # Subset data using valid row names and compute column means
    df <- data[valid_rows, , drop = FALSE]
    df <- as.data.frame(colMeans(df, na.rm = TRUE))
    colnames(df) <- lineage
  } else {
    # Handle case where no valid rows are found
    df <- data.frame(matrix(NA, nrow = ncol(data), ncol = 1))
    colnames(df) <- lineage
  }
  
  return(df)
  
  
}))

###########export heatmap###########################################

# Define a custom color palette with grey for zero values
custom_colors <- colorRampPalette(c("#003f5c", "#2f4b7c", "#bc5090", "#ef5675", "#ffa600"))(100)  # Create a palette for non-zero values
custom_colors <- c("light grey", custom_colors)  # Add grey for the zero value

# Define the row order
dataset_order <- c("Ai_model", "Hislop", "Liu", "Oldak", "Pedroza", "Weatherbee", "Rowan", "zheng_2019", "zheng_2022")

data=correlation_human_ref_celltype
# Reorder the matrix rows
data <- t(data)
#data <- data[,dataset_order, drop = FALSE]

# Open a PDF device to save the plot
pdf(paste0("correlation_human_ref_celltype", "_heatmap.pdf"), width = 7, height = 5)  # Save as PDF, adjust size if needed

pheatmap(data,
         scale = "none",  # Optional: scale rows to have mean 0 and SD 1
         color = custom_colors,  # Use the custom color palette
         cluster_rows = FALSE,  # Add row dendrogram
         cluster_cols = TRUE,  # Add column dendrogram
         display_numbers = FALSE,  # Optional: show expression values on heatmap
         fontsize_row = 10,  # Adjust font size for row names
         fontsize_col = 10,  # Adjust font size for column names
         main = "correlation_human_ref_celltype" # Use the extracted name as title
)

# Close the PDF device
dev.off()

######################## load lineage similarity tables#############################################
human_sim <- read.csv("correlation_human_ref_lineage.csv", row.names = 1)

rownames(human_sim) <- gsub("RNA.","",rownames(human_sim))
rownames(human_sim) <- gsub("\\.","_",rownames(human_sim))

##########human ref###################
data=human_sim
# Ensure 'lineages' is a factor to define the desired order
lineages <- factor(lineages, levels = lineages)

# Check for missing lineages in the row names of the data
missing_lineages <- setdiff(levels(lineages), rownames(data))

# Add rows for missing lineages with all values set to 0
if (length(missing_lineages) > 0) {
  missing_rows <- matrix(0, nrow = length(missing_lineages), ncol = ncol(data))
  rownames(missing_rows) <- missing_lineages
  colnames(missing_rows) <- colnames(data)
  
  # Bind the missing rows to the existing data
  data <- rbind(data, missing_rows)
}

# Reorder rows of the data based on the 'lineages' order
data <- data[as.character(lineages), ,drop = FALSE]

# Open a PDF device to save the plot
pdf(paste0("correlation_human_ref_lineage", "_heatmap.pdf"), width = 7, height = 5)  # Save as PDF, adjust size if needed

pheatmap(data,
         scale = "none",  # Optional: scale rows to have mean 0 and SD 1
         color = custom_colors,  # Use the custom color palette
         cluster_rows = FALSE,  # Add row dendrogram
         cluster_cols = TRUE,  # Add column dendrogram
         display_numbers = FALSE,  # Optional: show expression values on heatmap
         fontsize_row = 10,  # Adjust font size for row names
         fontsize_col = 10,  # Adjust font size for column names
         main = "correlation_human_ref_lineage" # Use the extracted name as title
)

# Close the PDF device
dev.off()
