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


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250806_embryo_model_benchmarking/cell_type_correlation_analysis")
outputDir = getwd()

# lineage

df <- read.csv("ranking_composite.csv")
# ranking by Frequency and then Composite_Rank_Score 
df_sorted <- df[order(-df$Frequency, -df$Composite_Rank_Score), ]

# check data after reorder
print(df_sorted)

# 1. Filter rows where Frequency is 9
top_pairs <- df_sorted[df_sorted$Frequency == 9, ]

#top_pairs <- head(df_sorted, 30)

# 2. Create a bar plot of Mean_Correlation between cell type pairs
library(ggplot2)

# Create pair labels first
top_pairs$Pair <- paste(top_pairs$CellType1, top_pairs$CellType2, sep = " - ")

ggplot(top_pairs, aes(x = reorder(Pair, Mean_Correlation), 
                      y = Mean_Correlation, 
                      fill = Mean_Correlation)) +
  geom_col() +
  geom_errorbar(aes(ymin = Mean_Correlation - Std_Correlation,
                    ymax = Mean_Correlation + Std_Correlation),
                width = 0.2) +
  scale_fill_gradient(low = "#6A0DAD", high = "#FFA500",  # Purple to Orange
                      name = "Correlation") +
  labs(title = "Cell Type Correlations (Frequency = 9)") +
  coord_flip() +
  theme_minimal()


ggsave("cell_type_correlations.pdf", width = 10, height = 6, dpi = 300)

##################### for specific cell types ####################################

# define plotting function
plot_top3_correlations <- function(cell_type, 
                                   data = df_sorted, 
                                   color_low = "#4B0092", 
                                   color_high = "#F28E2B",
                                   corr_limits = c(0.6, 1.0),
                                   save_plot = TRUE,
                                   output_dir = "plots",
                                   file_format = "pdf",
                                   width = 8,
                                   height = 6,
                                   dpi = 300) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Subset top 3 pairs involving the specified cell type
  top3_pairs <- data %>%
    filter(CellType1 == cell_type | CellType2 == cell_type) %>%
    head(3)
  
  # 2. Create clean pair labels
  top3_pairs <- top3_pairs %>%
    mutate(Pair = ifelse(CellType1 == cell_type,
                         paste(cell_type, CellType2, sep = " - "),
                         paste(CellType1, cell_type, sep = " - ")))
  
  # 3. Generate the plot
  p <- ggplot(top3_pairs, aes(x = reorder(Pair, Mean_Correlation), 
                              y = Mean_Correlation,
                              fill = Mean_Correlation)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = Mean_Correlation - Std_Correlation,
                      ymax = Mean_Correlation + Std_Correlation),
                  width = 0.2, color = "gray30", linewidth = 0.5) +
    scale_fill_gradient(low = color_low, high = color_high,
                        limits = corr_limits,
                        name = "Correlation") +
    labs(title = paste("Top 3 Correlated Pairs with", cell_type),
         subtitle = "Sorted by Composite Rank Score",
         x = NULL,
         y = "Mean Correlation") +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      panel.grid.major.y = element_blank()
    ) +
    geom_text(aes(label = round(Mean_Correlation, 3)), 
              hjust = -0.2, size = 3.5)
  
  # 4. Save the plot if requested
  if (save_plot) {
    # Create filename-safe version of cell type
    safe_name <- gsub("[^[:alnum:]]", "_", cell_type)
    filename <- paste0("top3_correlations_", safe_name, ".", file_format)
    filepath <- file.path(output_dir, filename)
    
    # Choose the appropriate saving function
    save_func <- switch(file_format,
                        "pdf" = ggsave,
                        "png" = function(...) ggsave(..., dpi = dpi),
                        "tiff" = function(...) ggsave(..., dpi = dpi),
                        "svg" = ggsave,
                        ggsave)  # default
    
    save_func(filepath, plot = p, width = width, height = height)
    message("Plot saved to: ", filepath)
  }
  
  return(p)
}

# PGC
plot_top3_correlations("PGC")
plot_top3_correlations("Cardiac.mesoderm")
plot_top3_correlations("Endothelium")
plot_top3_correlations("Spinal.cord")
plot_top3_correlations("DE")
plot_top3_correlations("Primitive.streak")
plot_top3_correlations("AVE")
plot_top3_correlations("VE")

#######################################################
library(dplyr)
library(ggplot2)

# 1. First create a unified table of all cell type pairs
all_pairs <- bind_rows(
  # Treat CellType1 as the main cell type
  df_sorted %>% 
    rename(main_cell = CellType1, partner_cell = CellType2) %>%
    select(main_cell, partner_cell, Frequency, Mean_Correlation, Std_Correlation, Composite_Rank_Score),
  
  # Treat CellType2 as the main cell type
  df_sorted %>% 
    rename(main_cell = CellType2, partner_cell = CellType1) %>%
    select(main_cell, partner_cell, Frequency, Mean_Correlation, Std_Correlation, Composite_Rank_Score)
) %>%
  distinct()  # Remove any duplicate pairs

# 2. Get top 1 pair for each cell type, ordered by Frequency
top1_table <- all_pairs %>%
  group_by(main_cell) %>%
  arrange(desc(Composite_Rank_Score)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(Frequency))  # Sort by Frequency for plotting

# 3. Create the plot with split labels
ggplot(top1_table, aes(x = reorder(main_cell, Frequency))) +
  # Left-side labels (main cell types)
  geom_text(aes(y = -0.2, label = main_cell), 
            hjust = 1, size = 4, fontface = "bold") +
  # Right-side labels (partner cell types)
  geom_text(aes(y = 1.1, label = partner_cell), 
            hjust = 0, size = 4) +
  # Correlation bars
  geom_col(aes(y = Mean_Correlation, fill = Mean_Correlation), 
           width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Correlation - Std_Correlation,
                    ymax = Mean_Correlation + Std_Correlation),
                width = 0.2, color = "gray30", linewidth = 0.5) +
  # Color scale
  scale_fill_gradient(low = "#4B0092", high = "#F28E2B",
                      limits = c(0.8, 1.0), name = "Correlation") +
  # Adjust plot area
  scale_y_continuous(limits = c(-0.5, 1.2), expand = c(0, 0)) +
  labs(title = "Top Correlated Pair for Each Cell Type",
       subtitle = "Ordered by Frequency | Left: Cell Type, Right: Top Partner",
       x = NULL, y = "Mean Correlation") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

# 4. Save the plot
ggsave("top_pairs_split_labels.pdf", width = 10, height = 8, dpi = 300)