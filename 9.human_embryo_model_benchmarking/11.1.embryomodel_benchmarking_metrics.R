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

########1.calculate percentage of high-confidence cells for benchmarking#####################################
# define threshold
threshold=0.2

certain_lineage <- list()
certain_celltype <- list()

for (table_name in names(obj_new)) {
  
  dataset <- obj_new[[table_name]]
  
  certain_lineage[[table_name]] <- length(which(dataset$human_ref_transferred_final_lineage_uncert < threshold))/length(dataset$human_ref_transferred_final_lineage_uncert)
  certain_celltype[[table_name]] <- length(which(dataset$human_ref_transferred_final_anno_uncert<threshold))/length(dataset$human_ref_transferred_final_anno_uncert)
  
}

# Combine all results into a single data frame
certain_lineage <- do.call(rbind, certain_lineage)
certain_celltype <- do.call(rbind, certain_celltype)

metrics <- cbind(certain_lineage,certain_celltype)
colnames(metrics) <- c("certain_lineage","certain_celltype")

#######2.lineage and cell_type coverge###############

cov_lineage <- list()
cov_celltype <- list()

for (table_name in names(obj_new)) {
  
  dataset <- obj_new[[table_name]]
  
  #subset certain lineage
  certain.lineage <- subset(dataset,subset = human_ref_transferred_final_lineage_uncert < threshold)
  cov_lineage[[table_name]] <- length(unique(certain.lineage$human_ref_transferred_final_lineage_unfiltered))/length(lineage_ordered_labels)
 
  #subset certain cell type
  certain.celltype <- subset(dataset,subset = human_ref_transferred_final_anno_uncert < threshold)
  cov_celltype[[table_name]] <- length(unique(certain.celltype$human_ref_transferred_final_anno_unfiltered))/length(humanref_ordered_labels)
  
}

# Combine all results into a single data frame
cov_lineage <- do.call(rbind, cov_lineage)
cov_celltype <- do.call(rbind, cov_celltype)

cov_merge <- cbind(cov_lineage,cov_celltype)
colnames(cov_merge) <- c("cov_lineage","cov_celltype")

metrics <- cbind(metrics, cov_merge)

#######plot percentage#################
# Generate heatmap with proper color breaks
library(pheatmap)

# Define a custom color palette with grey for zero values
custom_colors <- colorRampPalette(c("#003f5c", "#2f4b7c", "#bc5090", "#ef5675", "#ffa600"))(100)  # Create a palette for non-zero values

# Define the row order
dataset_order <- c("Ai_model", "Hislop", "Liu", "Oldak", "Pedroza", "Weatherbee", "Rowan", "zheng_2019", "zheng_2022")

# Reorder the matrix rows
cov_merge <- cov_merge[dataset_order, , drop = FALSE]
# Melt the dataframe to reorganize it
cov_melted <- melt(cov_merge, id.vars = "Model", 
                   variable.name = "Coverage_Type", 
                   value.name = "Coverage_Value")
Var2 <- as.factor(cov_melted$Var2)

# Create the boxplot
p <- ggplot(cov_melted, aes(x = Var1, y = Coverage_Value, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Var2, scales = "free_y", ncol = 1) +  # Wrap by Coverage_Type
  theme_minimal(base_size = 14) +
  labs(
    title = "Coverage by Lineage and Cell Type",
    x = "Model",
    y = "Coverage Value"
  ) +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none"  # Remove legend for cleaner visualization
  )+
  scale_fill_manual(
    values = c("cov_lineage" = "#f26b21", "cov_celltype" = "#68b8be"),  # Custom colors
  )

p

# Save the plot as a PDF
ggsave("embryo_model_coverage.pdf", plot = p, device = "pdf", width = 4, height = 4)

############3.calculate human lineage similarity#################################

#load human reference data
human <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250107_reference_comparasion/human_20250108.Rds")
#human <- NormalizeData(human,normalization.method = "LogNormalize", scale.factor = 10000) #already normalized
#error check
unique(human$final_anno)

###calculate lineage mean of each gene
human$final_lineage <- as.factor(human$final_lineage)
Idents(human) <- human$final_lineage
# Calculate average expression using AggregateExpression
avgexp <- AggregateExpression(human, assays = "RNA", slot = "data")
avgexp <- avgexp %>% as.data.frame()

similarity <- list()
similarity[["human_ref_lineage"]] <- avgexp

###calculate cell type mean of each gene

human$final_anno <- as.factor(human$final_anno)
Idents(human) <- human$final_anno
# Calculate average expression using AggregateExpression
avgexp <- AggregateExpression(human, assays = "RNA", slot = "data")
avgexp <- avgexp %>% as.data.frame()

similarity[["human_ref_celltype"]] <- avgexp

###calculate lineage and cell type mean of each gene for each dataset###########

for (table_name in names(obj_new)) {
  
  dataset <- obj_new[[table_name]]
  
  #subset certain lineage
  dataset <- subset(dataset,subset = human_ref_transferred_final_lineage_uncert < threshold)
  
  ##calculate humanref_lineage mean
  dataset$human_ref_transferred_final_lineage_unfiltered <- as.factor(dataset$human_ref_transferred_final_lineage_unfiltered)
  Idents(dataset) <- dataset$human_ref_transferred_final_lineage_unfiltered
  # Calculate average expression using AggregateExpression
  avgexp <- AggregateExpression(dataset, assays = "RNA", slot = "data")
  avgexp <- avgexp %>% as.data.frame()
  
  similarity[[paste0(table_name,"_human_ref_lineage")]] <- avgexp
  
  ######################################################################
  
  dataset <- obj_new[[table_name]]
  
  #subset certain cell type
  dataset <- subset(dataset,subset = human_ref_transferred_final_anno_uncert < threshold)
  
  ##calculate humanref_celltype mean
  dataset$human_ref_transferred_final_anno_unfiltered <- as.factor(dataset$human_ref_transferred_final_anno_unfiltered)
  Idents(dataset) <- dataset$human_ref_transferred_final_anno_unfiltered
  # Calculate average expression using AggregateExpression
  avgexp <- AggregateExpression(dataset, assays = "RNA", slot = "data")
  avgexp <- avgexp %>% as.data.frame()
  
  similarity[[paste0(table_name,"_human_ref_celltype")]] <- avgexp
  
}

################calculate correlation matrix####################################
#seperate matrics and organize lineage lists
lineage <- names(similarity)[grep("human_ref_lineage", names(similarity))]
human_ref_lineage<- similarity[[lineage[1]]] # reference
query_human_ref_lineage<- similarity[lineage[c(2:10)]]   # query list

#seperate matrics and organize lists
celltype <- names(similarity)[grep("human_ref_celltype", names(similarity))]
human_ref_celltype<- similarity[[celltype[1]]]  # reference
query_human_ref_celltype<- similarity[celltype[c(2:10)]]   # query list


###############calculate 1.human_ref_lineage similarity############

#define ref and query list
ref= human_ref_lineage
query= query_human_ref_lineage

correlation_human_ref_lineage <- list()
# Loop through all matrices in the list (starting from the second one)
for (i in 1:length(query)) {
  
  current_matrix <- query[[i]]
  # Find common row names
  common_genes <- intersect(rownames(ref), rownames(current_matrix))
  
  # Subset both matrices to keep only common genes
  ref_subset <- ref[common_genes, , drop = FALSE]
  curr_subset <- current_matrix[common_genes, , drop = FALSE]
  
  # Check if the subsets are valid
  if (ncol(ref_subset) > 0 && ncol(curr_subset) > 0) {
    # Calculate pairwise correlation using only the common groups
    corr_matrix <- cor(ref_subset, curr_subset, method = "pearson", use = "pairwise.complete.obs")
    
  } else {
    warning("No common groups found for correlation.")
  }
  
  # Store the result in the list
  correlation_human_ref_lineage[[paste0("Correlation_with_", names(query)[i])]] <- corr_matrix
}


###############calculate 2.human_ref_celltype similarity############

#define ref and query list
ref= human_ref_celltype
query= query_human_ref_celltype

correlation_human_ref_celltype <- list()
# Loop through all matrices in the list (starting from the second one)
for (i in 1:length(query)) {
  
  current_matrix <- query[[i]]
  # Find common row names
  common_genes <- intersect(rownames(ref), rownames(current_matrix))
  
  # Subset both matrices to keep only common genes
  ref_subset <- ref[common_genes, , drop = FALSE]
  curr_subset <- current_matrix[common_genes, , drop = FALSE]
  
  # Check if the subsets are valid
  if (ncol(ref_subset) > 0 && ncol(curr_subset) > 0) {
    # Calculate pairwise correlation using only the common groups
    corr_matrix <- cor(ref_subset, curr_subset, method = "pearson", use = "pairwise.complete.obs")
    
  } else {
    warning("No common groups found for correlation.")
  }
  
  # Store the result in the list
  correlation_human_ref_celltype[[paste0("Correlation_with_", names(query)[i])]] <- corr_matrix
}


#########################calculate 3.mean similarity################################

#1.correlation_human_ref_lineage

query = correlation_human_ref_lineage

datalist <-list()

for (x in 1:length(query)) {
  
  #extract name
  data <- query[[x]]
  name <- names(query)[x]
  name <- sub("Correlation_with_", "", name)
  name <- sub("_human_ref_lineage", "", name)
  
  #generate df
  df <- do.call(cbind, lapply(1:ncol(data), function(i) {
    
    #only keep values of the same group comparison
    label <- colnames(data)[i]
    df <- as.data.frame(data[label,i])
    colnames(df) <- label
    return(df)
    
  }))
  
  df <- as.data.frame(t(df))
  colnames(df) <- name 
  
  print(paste(name,"can be excuted!"))
  datalist[[x]] <- df
  
}

#generate one dataframe from the list
# Define the function
combine_dataframes <- function(df_list) {
  # Get all unique row names from all data frames
  all_rownames <- unique(unlist(lapply(df_list, rownames)))
  
  # Initialize an empty list to store modified data frames
  combined_list <- list()
  
  # Loop through each data frame in the list
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    
    # Convert row names to a column
    df_with_rownames <- df %>% 
      rownames_to_column(var = "rowname") 
    
    # Create a new data frame with all row names and fill missing values with 0
    full_df <- data.frame(rowname = all_rownames)
    full_df <- full_df %>%
      left_join(df_with_rownames, by = "rowname") %>%
      replace(is.na(.), 0)  # Replace NAs with 0
    
    # Assign the data to the combined list (excluding the rowname column)
    combined_list[[i]] <- full_df[-1]  # Remove the rowname column
  }
  
  # Combine all data frames with cbind
  final_combined_df <- do.call(cbind, combined_list)
  
  # Set the row names to the combined row names
  rownames(final_combined_df) <- all_rownames
  
  return(final_combined_df)
}

correlation_human_ref_lineage <- combine_dataframes(datalist)


#2.correlation_human_ref_celltype

query = correlation_human_ref_celltype

datalist <-list()

for (x in 1:length(query)) {
  
  #extract name
  data <- query[[x]]
  name <- names(query)[x]
  name <- sub("Correlation_with_", "", name)
  name <- sub("_human_ref_celltype", "", name)
  
  #generate df
  df <- do.call(cbind, lapply(1:ncol(data), function(i) {
    
    #only keep values of the same group comparison
    label <- colnames(data)[i]
    df <- as.data.frame(data[label,i])
    colnames(df) <- label
    return(df)
    
  }))
  
  df <- as.data.frame(t(df))
  colnames(df) <- name 
  
  print(paste(name,"can be excuted!"))
  
  datalist[[x]] <- df
  
}

correlation_human_ref_celltype <- combine_dataframes(datalist)

########################generate heatmaps base on these correlation metrics####################

custom_colors <- c("light grey", custom_colors)  # Add grey for the zero value

## 1. correlation_human_ref_lineage
data.list <- list(
  correlation_human_ref_lineage = correlation_human_ref_lineage,
  correlation_human_ref_celltype = correlation_human_ref_celltype
)

for (i in seq_along(data.list)) {
  
  data.name <- names(data.list)[i]
  data <- data.list[[i]]  # Reference the current dataframe
  
  # Open a PDF device to save the plot
  pdf(paste0(data.name, "_heatmap.pdf"), width = 7, height = 5)  # Save as PDF, adjust size if needed
  
  pheatmap(data,
           scale = "none",  # Optional: scale rows to have mean 0 and SD 1
           color = custom_colors,  # Use the custom color palette
           cluster_rows = TRUE,  # Add row dendrogram
           cluster_cols = TRUE,  # Add column dendrogram
           display_numbers = FALSE,  # Optional: show expression values on heatmap
           fontsize_row = 10,  # Adjust font size for row names
           fontsize_col = 10,  # Adjust font size for column names
           main = data.name # Use the extracted name as title
  )
  
  # Close the PDF device
  dev.flush()  #
  dev.off()
  
  # print(p)  # Print the plot
}

# export similarity tables
write.csv(correlation_human_ref_celltype, file = "correlation_human_ref_celltype.csv", row.names = TRUE)
write.csv(correlation_human_ref_lineage, file = "correlation_human_ref_lineage.csv", row.names = TRUE)

#######################calculate mean similarity for benchmarking#####################################

mean_sim <- do.call(cbind, lapply(seq_along(data.list), function(i) {
  
  data.name <- names(data.list)[i]
  data <- data.list[[i]]  
  
  # Function to calculate column means excluding zeros
  col_means_exclude_zeros <- function(data) {
    apply(data, 2, function(x) mean(x[x != 0], na.rm = TRUE))
  }
  
  # Calculate column means excluding zeros
  means <- col_means_exclude_zeros(data)
  means <- as.data.frame(means)
  colnames(means)<- data.name
  
  return(means)
  
  
}))

## add mean_sim into metrics
metrics <- merge(metrics, mean_sim, by = "row.names", all = TRUE)
rownames(metrics) <- metrics$Row.names
metrics <- metrics[,2:7, drop=FALSE]


# Save the metrics dataframe as a CSV file
write.csv(metrics, file = "embryomodel_metrics.csv", row.names = TRUE)


