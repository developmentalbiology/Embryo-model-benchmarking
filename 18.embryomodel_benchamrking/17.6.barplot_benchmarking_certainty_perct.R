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


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250806_embryo_model_benchmarking/mean_certainty_metrics_consistent_cells")
outputDir = getwd()

# lineage
df <- read.csv("lineage_overall_mean_certainty_consistent_cells_comparison.csv")

# Assign specific colors to each model
model_colors <- c(
  "corrected_processed_Hislop" = "#ffddc2",       # Red
  "corrected_processed_Liu" = "#f8a18d",         # Blue
  "corrected_processed_Oldak" = "#f25741",       # Green
  "corrected_processed_zheng_2022" = "#e9d767",  # Purple
  "corrected_processed_Pedroza" = "#9ea517",     # Orange
  "corrected_processed_Ai_model" = "#bbe0ce",    # Yellow
  "corrected_processed_Rowan" = "#587c7c",       # Brown
  "corrected_processed_Weatherbee" = "#01405f",  # Pink
  "corrected_processed_zheng_2019" = "#017d85"   # Gray
)


# Create the bar plot
plot <- ggplot(df, aes(x = reorder(Model, Mean_Certainty), y = Mean_Certainty, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) +  # Bar plot
  scale_fill_manual(values = model_colors) +  # Apply custom color palette
  geom_text(aes(label = paste0(round(Mean_Certainty, 2))), 
            position = position_stack(vjust = 0.5), size = 4) +  # Add percentage labels
  labs(
    title = "Certainty Percentage by Model",
    x = "Model",
    y = "Mean_Certainty"
  ) +
  theme_minimal() +  # Use a clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"  # Remove legend since colors are self-explanatory
  ) +
  coord_flip()  # Flip coordinates for horizontal bars

# Export the plot as a PDF
pdf("certainty_consistent_lineage.pdf", width = 10, height = 6)  # Set dimensions of the PDF
print(plot)  # Print the plot to the PDF file
dev.off()  # Close the PDF device



# cell type
df <- read.csv("celltype_overall_mean_certainty_consistent_cells_comparison.csv")

# Create the bar plot
plot <- ggplot(df, aes(x = reorder(Model, Mean_Certainty), y = Mean_Certainty, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) +  # Bar plot
  scale_fill_manual(values = model_colors) +  # Apply custom color palette
  geom_text(aes(label = paste0(round(Mean_Certainty, 2))), 
            position = position_stack(vjust = 0.5), size = 4) +  # Add percentage labels
  labs(
    title = "Certainty Percentage by Model",
    x = "Model",
    y = "Mean_Certainty"
  ) +
  theme_minimal() +  # Use a clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"  # Remove legend since colors are self-explanatory
  ) +
  coord_flip()  # Flip coordinates for horizontal bars

# Export the plot as a PDF
pdf("certainty_consistent_celltype.pdf", width = 10, height = 6)  # Set dimensions of the PDF
print(plot)  # Print the plot to the PDF file
dev.off()  # Close the PDF device


