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


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250117_label_transfer_Weatherbee_comparasion")
outputDir = getwd()



#load garfield result
garfield <- read.csv("corrected_processed_Weatherbee.h5ad_garfield_query.obs.csv")
unique(garfield$cell_assignment)

weatherbee_sum <- garfield[,c("X","cell_assignment", "transferred_final_anno_unfiltered","transferred_final_lineage_unfiltered"),drop=FALSE]

#load scArches result
scarches <- read.csv("human_Weatherbee_scArches.csv")
scarches_sum <- scarches[,c("X", "scArches_final_anno_pre" ,"scArches_final_lineage_pre"),drop=FALSE]

#load scGPT result
scgpt <- read.csv("human_Weatherbee_scgpt.csv")
scgpt_sum <- scgpt[,c("X", "scgpt_final_anno_pre", "scgpt_final_lineage_pre"),drop=FALSE]

#load querymap result
querymap <- read.csv("human_ref_weatherbee_querymap.csv")

#merge table
merged_all <- merge(weatherbee_sum, scarches_sum, by = 'X', all = TRUE)
merged_all <- merge(merged_all, scgpt_sum, by = 'X', all = TRUE)
merged_all <- merge(merged_all, querymap, by = 'X', all = TRUE)


#only focus on hypo/VE lineage to calculate prediction metrics
merged_all <- merged_all[which(merged_all$cell_assignment%in% c("HYPO/VE","L-EPI","AM-3" )  ),] 


# Define a function to rename values in a column
rename_values <- function(column) {
  column[column %in% c("AVE","VE/YE", "Hypoblast","YS.Endoderm_1","Gut", "ExE_endo", "DE","Endoderm")] <- "HYPO/VE"
  column[column %in% c( "Epi_1","Epi_2", "Epi_3", "Epi_4","epi","Ectoderm")] <- "L-EPI"
  column[column %in% c("Amniotic_epi", "Amnion","non_neuro_ecto")] <- "AM-3"
  return(column)
}


# Apply the renaming function to the relevant columns
columns_to_rename <- colnames(merged_all)[3:10]

for (col in columns_to_rename) {
  merged_all[[col]] <- rename_values(merged_all[[col]])
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
results_df$method <- c("Garfield","Garfield","scArches","scArches","scGPT","scGPT","QueryMap","QueryMap")
results_df$label <- c("celltype","lineage","celltype","lineage","celltype","lineage","celltype","lineage")

results_df <- results_df[,c("Accuracy", "Precision","method","label")]

##reorganzie the dataframe
results_df <- results_df %>%
  pivot_longer(cols = c(Accuracy, Precision), 
               names_to = "Metric", 
               values_to = "Value")


# Define the color palette you want to use
color_palette <- c("Garfield" = "#FDB462", "QueryMap" = "#B2DF8A", "scArches" = "#A6CEE3", "scGPT" = "#7570B3")

# Plot
results_df$Value <- as.numeric(results_df$Value)
results_df$method <- factor(results_df$method, levels = c("Garfield","scArches","QueryMap","scGPT"))

p<- ggplot(results_df, aes(x = Metric, y = Value, fill = method)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~ label, strip.position = "top") +  # Use 'label' for facetting
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = color_palette) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Mean Value")

p

# Save the plot as a PDF
ggsave("label_reansfer_human_ref.pdf", plot = p, device = "pdf", width = 8, height = 6)


