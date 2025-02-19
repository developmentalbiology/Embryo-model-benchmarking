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

# Define the ordered labels

humanref_ordered_labels = c('TE','CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'STBs_3', 'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4',
                            'Epi_1', 'Epi_2', 'Epi_3','Epi_4',
                            'Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm',  'Exe.endothelium', 
                            'Amnion', 'Amniotic_epi',  'Ectoderm_1', 'Ectoderm_2',
                            'Neural tube', 'Neural crest',
                            'Primitive.streak', 'Nascent mesoderm','PGC',
                            'Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
                            'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4',
                            'Lateral plate mesoderm_5', 'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm',
                            'Cardiac myocyte', 
                            'Notochord', 'DE', 'Gut',
                            'Hypoblast', 'AVE', 'VE/YE', 'YS.Endoderm_1', 'YS.Endoderm_2', 
                            'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor'
                            
)


# Define the ordered labels
lineage_ordered_labels = c("Endoderm","epi","ExE_endo","Exe_meso", "hemogenic","mesoderm", "neural_ecto", "non_neuro_ecto","Notochord",  "PGC",  "TE_TrB", "Gastru")


metrics <- list()

########1.calculate percentage of high-confidence cells for benchmarking#####################################
# check weather all certain cells belong to certain lineage

# Loop over each dataset in obj_new
for (table_name in names(obj_new)) {
  
  dataset <- obj_new[[table_name]]
  
  # Subset rownames based on the threshold
  cert_celltype <- rownames(dataset)[which(dataset$human_ref_transferred_final_anno_uncert < threshold)]
  cert_lineage <- rownames(dataset)[which(dataset$human_ref_transferred_final_lineage_uncert < threshold)]
  
  # Calculate the differences
  celltype_diff <- setdiff(cert_celltype, cert_lineage)
  lineage_diff <- setdiff(cert_lineage, cert_celltype)
  
  # Print the results for this table_name
  cat("\nFor table:", table_name, "\n")
  cat("Row names in cert_celltype but not in cert_lineage:", length(celltype_diff), "elements\n")
  cat("Row names in cert_lineage but not in cert_celltype:", length(lineage_diff), "elements\n")
}

##################################

# Define threshold
threshold = 0.2

# Initialize lists to store results
certain_lineage <- list()
certain_celltype_lineage <- list()

# Loop over each dataset in obj_new
for (table_name in names(obj_new)) {
  
  dataset <- obj_new[[table_name]]
  metadata <- dataset@meta.data
  
  ################### Subset certain lineage #######################
  
  avg_lineage <- do.call(rbind, lapply(lineage_ordered_labels, function(lineage) {
    
    # Subset data for the current lineage
    data <- metadata[metadata$human_ref_transferred_final_lineage_unfiltered == lineage, ]
    
    # Calculate the percentage of cells with 'human_ref_transferred_final_lineage_uncert' less than the threshold
    data_lineage <- length(which(data$human_ref_transferred_final_lineage_uncert < threshold)) / nrow(data)
    
    # Return the lineage and the calculated percentage as a dataframe row
    return(data.frame(lineage = lineage, percentage = data_lineage))
    
  }))
  
  # Store the results in the list with table_name as the name of the column
  # Ensure lineage column is not duplicated
  avg_lineage$lineage <- as.character(avg_lineage$lineage)
  colnames(avg_lineage)[2] <- table_name  # Set the column name to the table_name
  
  certain_lineage[[table_name]] <- avg_lineage
  
  
  ################# Subset certain cell type ##################
  
  avg_celltype_lineage <- do.call(rbind, lapply(lineage_ordered_labels, function(lineage) {
    
    # Subset data for the current lineage
    data <- metadata[metadata$human_ref_transferred_final_lineage_unfiltered == lineage, ]
    
    # Calculate the percentage of cells with 'human_ref_transferred_final_anno_uncert' less than the threshold
    data_celltype <- length(which(data$human_ref_transferred_final_anno_uncert < threshold)) / nrow(data)
    
    # Return the lineage and the calculated percentage as a dataframe row
    return(data.frame(lineage = lineage, percentage = data_celltype))
    
  }))
  
  # Store the results in the list with table_name as the name of the column
  # Ensure lineage column is not duplicated
  avg_celltype_lineage$lineage <- as.character(avg_celltype_lineage$lineage)
  colnames(avg_celltype_lineage)[2] <- table_name  # Set the column name to the table_name
  
  certain_celltype_lineage[[table_name]] <- avg_celltype_lineage
  
}

# Combine results into final dataframes, ensuring no repeated "lineage" column
final_certain_lineage <- Reduce(function(x, y) merge(x, y, by = "lineage", all = TRUE), certain_lineage)
final_certain_celltype_lineage <- Reduce(function(x, y) merge(x, y, by = "lineage", all = TRUE), certain_celltype_lineage)

metrics [["perct.certain_lineage"]] <- final_certain_lineage
metrics [["perct.certain_celltype_lineage"]] <- final_certain_celltype_lineage

######################2.calculate lineage coverage for benchmarking#####################################

coverage_lineage <- list()
coverage_celltype_lineage <- list()

for (table_name in names(obj_new)) {
  
  dataset <- obj_new[[table_name]]
  
  # subset certain lineage
  
  dataset <- subset(dataset,subset = human_ref_transferred_final_lineage_uncert < threshold)
  metadata <- dataset@meta.data
  
  ################### calculate lineage coverage #######################
  # Create a dataframe with NA values of the same length as lineage_ordered_labels
  cov_lineage <- data.frame(
    lineage = lineage_ordered_labels,
    coverage = rep(NA, length(lineage_ordered_labels)),
    stringsAsFactors = FALSE
  )
  # Get unique lineages
  unique_lineages <- unique(metadata$human_ref_transferred_final_lineage_unfiltered)
  
  # Assign the value 1 to the coverage column for each unique lineage
  cov_lineage$coverage[which(cov_lineage$lineage %in% unique_lineages)] <- 1
  cov_lineage$coverage[which(is.na(cov_lineage$coverage))] <- 0
 
  colnames(cov_lineage)[2] <- table_name  # Set the column name to the table_name
  coverage_lineage[[table_name]] <- cov_lineage
    
  ################### calculate cell type coverage #######################
  
  dataset <- subset(dataset, subset = human_ref_transferred_final_anno_uncert < threshold)
  metadata <- dataset@meta.data
  
  cov_celltype_lineage <- do.call(rbind, lapply(lineage_ordered_labels, function(lineage) {
    
    # get human cell type label 
    lineage_celltype <- unique(human$final_anno[which(human$final_lineage==lineage)])
    
    # Subset data for the current lineage
    data <- metadata[metadata$human_ref_transferred_final_anno_unfiltered %in% lineage_celltype, ]
    
    # Calculate the percentage of cells with 'human_ref_transferred_final_anno_uncert' less than the threshold
    
    data_celltype <- length(unique(data$human_ref_transferred_final_anno_unfiltered))/length(lineage_celltype)
    
    # Return the lineage and the calculated percentage as a dataframe row
    return(data.frame(lineage = lineage, coverage = data_celltype))
    
  }))
  
  colnames(cov_celltype_lineage)[2] <- table_name  # Set the column name to the table_name
  coverage_celltype_lineage[[table_name]] <- cov_celltype_lineage
    
  }
  
# Combine results into final dataframes, ensuring no repeated "lineage" column
final_coverage_lineage <- Reduce(function(x, y) merge(x, y, by = "lineage", all = TRUE), coverage_lineage)
final_coverage_celltype_lineage <- Reduce(function(x, y) merge(x, y, by = "lineage", all = TRUE), coverage_celltype_lineage)
  
metrics [["coverage_lineage"]] <- final_coverage_lineage
metrics [["coverage_celltype_lineage"]] <- final_coverage_celltype_lineage
  

######################3.include lineage similarity for benchmarking#####################################

# load lineage similarity tables
human_sim <- read.csv("correlation_human_ref_lineage.csv", row.names = 1)

rownames(human_sim) <- gsub("RNA.","",rownames(human_sim))
rownames(human_sim) <- gsub("\\.","_",rownames(human_sim))

data=human_sim
# Ensure 'lineages' is a factor to define the desired order
lineage_ordered_labels <- factor(lineage_ordered_labels, levels = lineage_ordered_labels)

# Check for missing lineages in the row names of the data
missing_lineages <- setdiff(levels(lineage_ordered_labels),rownames(data))

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

data <- data.frame(
  lineage=rownames(data),
  data
)

metrics [["similarity_lineage"]] <- data


# load cell type similarity tables
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

correlation_human_ref_celltype <- t(correlation_human_ref_celltype)

correlation_human_ref_celltype <- data.frame(
  lineage=rownames(correlation_human_ref_celltype),
  correlation_human_ref_celltype
)

metrics [["similarity_celltype"]] <- correlation_human_ref_celltype


##########################organize metrice table#########################################
# Initialize a list to hold the reorganized dataframes
reorganized_metrics <- list()

# List of metric names
metric_names <- names(metrics)
# Loop through each metric and reorganize the corresponding dataframe
for (metric_name in metric_names) {
  
  # Get the current metric dataframe (make sure the name matches the dataframe in 'metrics')
  current_df <- metrics[[metric_name]]
  
  # Ensure 'lineage' column is the row name
  rownames(current_df) <- current_df$lineage  # Set lineage as row names
  current_df <- current_df[, -1]  # Remove the 'lineage' column from the dataframe
  
  # Transpose the dataframe so that 'lineage' becomes rows and the metric names become columns
  transposed_df <- t(current_df)
  
  # Store the transposed dataframe in the list
  reorganized_metrics[[metric_name]] <- transposed_df
}

# Generate the reorganized dataframe for each lineage
# Create an empty list to store the final dataframes
metric <- list()

# Loop through each lineage in lineage_ordered_labels
for (lineage in lineage_ordered_labels) {
  
  # Combine columns for each dataframe in reorganized_metrics
  df <- do.call(cbind, lapply(1:length(reorganized_metrics), function(i) {
    
    # Get the dataframe from the list
    temp_df <- reorganized_metrics[[i]]
    
    # Select the column corresponding to the lineage
    temp_df <- temp_df[, lineage, drop = FALSE]
    
    # Ensure the row names match before combining
    temp_df <- temp_df[rownames(temp_df) == rownames(reorganized_metrics[[1]]), , drop = FALSE]
    
    # Rename the column to the metric name (from reorganized_metrics)
    colnames(temp_df) <- names(reorganized_metrics)[i]
    
    return(temp_df)
  }))
  
  # Ensure the row names are consistent across all columns
  rownames(df) <- rownames(reorganized_metrics[[1]])
  
  # Store the combined dataframe in the metric list
  metric[[lineage]] <- df
}


#################export tables###########################
# Define the output folder name
output_folder <- "output_folder"
# Create the folder if it doesn't already exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}
# Loop through each lineage in the metric list
for (i in 1:length(metric)) {
  # Get the matrix for the current lineage
  mat <- metric[[i]]
  
  # Convert matrix to a data frame
  df <- as.data.frame(mat)
  
  # Define the file name based on the lineage
  file_name <- paste0(names(metric)[i], "_lineage.csv")
  
  # Define the full path to save the file in the output folder
  file_path <- file.path(output_folder, file_name)

  # Write the dataframe to a CSV file
  write.csv(df, file = file_path, row.names = TRUE)
}
