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


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3")
outputDir = getwd()

#load garfield result
garfield <- read.csv("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3/embryo_model_garfield/corrected_processed_Weatherbee.h5ad_garfield_query.obs.csv")
unique(garfield$cell_assignment)

weatherbee_sum <- garfield[,c("X","cell_assignment", "transferred_reanno_unfiltered","transferred_lineage_unfiltered"),drop=FALSE]

#load scArches result
scarches <- read.csv("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3/embryo_model_scArches/human_Weatherbee_scArches.csv")
scarches_sum <- scarches[,c("X", "scArches_final_anno_pre" ,"scArches_final_lineage_pre"),drop=FALSE]

#load scPoli result
scpoli <- read.csv("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3/embryo_model_scPoli/lineage/corrected_processed_Weatherbee.h5ad_scPoli_query.csv")
scpoli_2 <- read.csv("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3/embryo_model_scPoli/reanno/corrected_processed_Weatherbee.h5ad_scPoli_query.csv")

# Merge by rowname and manually add unique columns
scpoli_2 <- scpoli_2[match(scpoli$X, scpoli_2$X),]
scpoli_2 <- scpoli_2[,c("reanno_pred", "reanno_uncert")]
scpoli <- cbind(scpoli, scpoli_2)

scpoli_sum <- scpoli[,c("X","reanno_pred", "lineage_pred" ),drop=FALSE]

## fix scpoli cell names
scpoli_sum$X <- gsub("-1$", "", scpoli_sum$X)

## check whether cell names match
setdiff(scpoli_sum$X, scarches_sum$X)   

#load scGPT result
scgpt <- read.csv("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3/embryo_model_scGPT/human_Weatherbee_scgpt.csv")
scgpt_sum <- scgpt[,c("X", "scgpt_final_anno_pre", "scgpt_final_lineage_pre"),drop=FALSE]

#load querymap result
querymap <- read.csv("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3/embryo_model_querymap/Weatherbee_human_ref_querymap.csv")

#merge table
merged_all <- merge(weatherbee_sum, scarches_sum, by = 'X', all = TRUE)
merged_all <- merge(merged_all, scpoli_sum, by = 'X', all = TRUE)
merged_all <- merge(merged_all, scgpt_sum, by = 'X', all = TRUE)
merged_all <- merge(merged_all, querymap, by = 'X', all = TRUE)


#only focus on selected lineage to calculate prediction metrics
merged_all <- merged_all[which(merged_all$cell_assignment%in% c("HYPO/VE","L-EPI","AM-3","MESO-1","MESO-2","EXMC" )  ),] 


# Define a function to rename values in a column
rename_values <- function(column) {
  column[column %in% c("AVE","VE", "Hypoblast_1","Hypoblast_2", "ExE_endo", "DE","Endoderm","YS.endoderm")] <- "HYPO/VE"
  column[column %in% c( "Epi_1","Epi_2", "Epi_3", "Ectoderm","epi")] <- "L-EPI"
  column[column %in% c( 'Amniontic.ectoderm',"Amniotic_ecto")] <- "AM-3"
  column[column %in% c('Primitive.streak','NMP','Neuromesodermal.progenitor')] <- "MESO-1"
  column[column %in% c('Axial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 'Myocyte.progenitor', 'Lateral.plate.mesoderm_1',
                       'Lateral.plate.mesoderm_2','Cardiac.mesoderm','Connecting.stalk')] <- "MESO-2"
  column[column %in% c('Amniotic.mesoderm','Exe.meso.progenitor', 'YS.mesoderm_1', 'YS.mesoderm_2','YS.mesoderm_1', 'YS.mesoderm_2' )] <- "EXMC"
  
  
  return(column)
}



# Apply the renaming function to the relevant columns
columns_to_rename <- colnames(merged_all)[3:12]

for (col in columns_to_rename) {
  merged_all[[col]] <- rename_values(merged_all[[col]])
}


# align "meso_Exe.meso" in lineage columns with the "cell_assignment" 
lineage_cols <- grep("lineage", colnames(merged_all))

# Loop over each row
for (i in seq_len(nrow(merged_all))) {
  # Check if the row meets the condition
  if (any(merged_all[i, lineage_cols] == "meso_Exe.meso") &&
      merged_all[i, "cell_assignment"] %in% c("MESO-2", "EXMC")) {
    
    # Replace "meso_Exe.meso" with cell_assignment value in that row
    merged_all[i, lineage_cols] <- ifelse(
      merged_all[i, lineage_cols] == "meso_Exe.meso",
      merged_all[i, "cell_assignment"],
      merged_all[i, lineage_cols]
    )
  }
}


# Load required libraries
library(caret)

# Initialize an empty list to store results
results_list <- list()

# Loop through columns 3 to ncol(merged_all)
for (i in 3:ncol(merged_all)) {
  # Extract predicted and true labels
  preds <- merged_all[[2]]
  gt <- merged_all[[i]]
  
  # Convert to factors for caret functions
  all_levels <- unique(c(preds, gt))
  preds <- factor(preds, levels = all_levels)
  gt <- factor(gt, levels = all_levels)
  
  # Calculate confusion matrix
  conf_matrix <- confusionMatrix(preds, gt)
  
  # Extract metrics
  accuracy <- conf_matrix$overall['Accuracy']
  
  # Calculate precision, recall, and F1 score manually
  precision <- sapply(levels(preds), function(level) {
    tp <- sum(preds == level & gt == level)
    fp <- sum(preds == level & gt != level)
    if ((tp + fp) > 0) tp / (tp + fp) else NA
  })
  
  recall <- sapply(levels(preds), function(level) {
    tp <- sum(preds == level & gt == level)
    fn <- sum(preds != level & gt == level)
    if ((tp + fn) > 0) tp / (tp + fn) else NA
  })
  
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  # Calculate macro-averaged precision, recall, and F1 score
  macro_precision <- mean(precision, na.rm = TRUE)
  macro_recall <- mean(recall, na.rm = TRUE)
  macro_f1 <- mean(f1, na.rm = TRUE)
  
  # Store results in the list
  results_list[[i - 2]] <- c(accuracy, macro_precision, macro_recall, macro_f1)
}

# Create a dataframe from the results list
results_df <- do.call(rbind, results_list)

# Set column names using the relevant names from merged_all
colnames(results_df) <- c("Accuracy", "Precision", "Recall", "Macro_F1")

# Add the corresponding column names from merged_all
results_df <- cbind(Column = colnames(merged_all)[3:ncol(merged_all)], results_df)

# Print the results dataframe
print(results_df)

# Save the data as a CSV file
write.csv(results_df, file = "human_ref_summary.csv")

#########plot##################

results_df <- as.data.frame(results_df)
results_df$method <- c("Garfield","Garfield","scArches","scArches","scPoli","scPoli","scGPT","scGPT","QueryMap","QueryMap")
results_df$label <- c("celltype","lineage","celltype","lineage","celltype","lineage","celltype","lineage","celltype","lineage")

results_df <- results_df[,c("Accuracy", "Precision","method","label")]

##reorganzie the dataframe
results_df <- results_df %>%
  pivot_longer(cols = c(Accuracy, Precision), 
               names_to = "Metric", 
               values_to = "Value")

# Ensure 'Value' is numeric
results_df$Value <- as.numeric(results_df$Value) 

# Define the color palette you want to use
color_palette <- c("Garfield" = "#FDB462","scArches" = "#A6CEE3", "scPoli" = "#F898CB", "QueryMap" = "#B2DF8A",  "scGPT" = "#7570B3")

# Plot
# Separate the data into lineage and celltype subsets
lineage_df <- results_df %>% filter(label == "lineage")
celltype_df <- results_df %>% filter(label == "celltype")

# Function to generate bar plots
generate_barplots <- function(data, label_name) {
  metrics <- unique(data$Metric)  # Get unique metrics (Accuracy, Precision)
  
  for (metric in metrics) {
    # Subset the data for the current metric
    subset_df <- data %>% filter(Metric == metric)
    
    # Order 'method' by 'Value' in descending order
    subset_df$method <- reorder(subset_df$method, -subset_df$Value)  # Reorder methods
    
    subset_df$method <- factor(subset_df$method,levels = c("scPoli","scArches", "Garfield", "scGPT", "QueryMap" ))
    
    # Create the plot
    p <- ggplot(subset_df, aes(x = method, y = Value, fill = method)) +
      geom_col(position = position_dodge(width = 0.8)) +  # Adjust dodge width for better spacing
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      scale_fill_manual(values = color_palette) +  # Use the defined color palette
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        legend.position = "none",  # Remove legend since it's redundant
        panel.grid.major = element_blank(),  # Remove grid lines
        panel.grid.minor = element_blank()
      ) +
      labs(title = paste(metric, "Scores for", label_name), x = "", y = "Mean Value")  # Add title and axis labels
    
    # Save the plot as a PDF file
    pdf_file <- paste0(metric, "_", label_name, "_barplot.pdf")
    ggsave(filename = pdf_file, plot = p, width = 8, height = 6, device = "pdf")
    
    # Display the plot
    print(p)
  }
}

# Generate bar plots for lineage
generate_barplots(lineage_df, "Lineage")

# Generate bar plots for celltype
generate_barplots(celltype_df, "Celltype")