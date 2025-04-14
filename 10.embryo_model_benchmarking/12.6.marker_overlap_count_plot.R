# Load required libraries
library(tidyverse)
library(funkyheatmap)
library(RColorBrewer)

# Set working directory
setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250325_embryomodel_benchmarking/marker_overlap_metrics")
outputDir = getwd()

# Read data
df_lineage <- read.csv("combined_lineage_metrics.csv")
df_celltype <- read.csv("combined_celltype_metrics.csv")
df_lineage <- df_lineage[, c(1,8), drop=FALSE]
df_celltype <- df_celltype[, c(1,8), drop=FALSE]
df <- merge(df_lineage, df_celltype, by="model")
colnames(df) <- c("model", "lineage", "celltype")


# Assign specific colors to each model
model_colors <- c(
  "Hislop" = "#ffddc2",       # Red
  "Liu" = "#f8a18d",         # Blue
  "Oldak" = "#f25741",       # Green
  "zheng_2022" = "#e9d767",  # Purple
  "Pedroza" = "#9ea517",     # Orange
  "Ai_model" = "#bbe0ce",    # Yellow
  "Rowan" = "#587c7c",       # Brown
  "Weatherbee" = "#01405f",  # Pink
  "zheng_2019" = "#017d85"   # Gray
)


# Create the bar plot
plot <- ggplot(df, aes(x = reorder(model, lineage), y = lineage, fill = model)) +
  geom_bar(stat = "identity", width = 0.7) +  # Bar plot
  scale_fill_manual(values = model_colors) +  # Apply custom color palette
  geom_text(aes(label = paste0(round(lineage, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4) +  # Add percentage labels
  labs(
    title = "Certainty Percentage by Model",
    x = "Model",
    y = "marker overlap"
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
pdf("markeroverlap_count_lineage.pdf", width = 10, height = 6)  # Set dimensions of the PDF
print(plot)  # Print the plot to the PDF file
dev.off()  # Close the PDF device


# Create the bar plot
plot <- ggplot(df, aes(x = reorder(model, celltype), y = celltype, fill = model)) +
  geom_bar(stat = "identity", width = 0.7) +  # Bar plot
  scale_fill_manual(values = model_colors) +  # Apply custom color palette
  geom_text(aes(label = paste0(round(celltype, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4) +  # Add percentage labels
  labs(
    title = "Certainty Percentage by Model",
    x = "Model",
    y = "marker overlap"
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
pdf("markeroverlap_count_celltype.pdf", width = 10, height = 6)  # Set dimensions of the PDF
print(plot)  # Print the plot to the PDF file
dev.off()  # Close the PDF device

