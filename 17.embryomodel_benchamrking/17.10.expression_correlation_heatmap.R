# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)
library(viridis)
library(stringr)
library(scales)

setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250806_embryo_model_benchmarking/expression_correlation_metrics")
outputDir = getwd()

# Define the lineage to cell type mapping in the desired order
lineage_to_celltypes <- list(
  "TE_TrB" = c("TE", "CTB_1", "CTB_2", "STB_1", "STB_2", "STB_3", "EVT_1", "EVT_2"),
  "epi" = c("Epiblast_1", "Epiblast_2", "Epiblast_3", "Ectoderm"),
  "Primitive.streak" = c("Primitive.streak"),
  "NMP" = c("Neuromesodermal.progenitor"),
  "Notochord" = c("Notochord"),
  "PGC" = c("PGC"),
  "ExE_endo" = c("Hypoblast_1", "Hypoblast_2", "AVE", "VE", "YS.endoderm"),
  "Amniotic_ecto" = c("Amniontic.epi", "Amniontic.ectoderm"),
  "neural_ecto" = c("Neural.crest", "Neural.ectoderm.forebrain", "Neural.ectoderm.hindbrain",
                    "Neural.ectoderm.midbrain", "Spinal.cord"),
  "Endoderm" = c("DE", "DE_1", "DE_2", "Gut"),
  "meso_Exe.meso" = c("Paraxial.mesoderm", "Emergent.mesoderm", "Pre-somatic.mesoderm", "Somite",
                      "Rostral.mesoderm", "Lateral.plate.mesoderm_1", "Lateral.plate.mesoderm_2",
                      "Lateral.plate.mesoderm_3", "Cardiac.mesoderm", "Amniotic.mesoderm",
                      "Exe.meso.progenitor", "YS.mesoderm_1", "YS.mesoderm_2"),
  "hemogenic" = c("Hemogenic.endothelial.progenitor", "Endothelium", "Erythroid",
                  "Primitive.megakaryocyte", "Myeloid.progenitor")
)

# Create a flat mapping from cell type to lineage while maintaining order
all_cell_types <- unlist(lineage_to_celltypes)
all_lineages <- rep(names(lineage_to_celltypes), sapply(lineage_to_celltypes, length))
cell_order_df <- data.frame(
  cell_type = all_cell_types,
  lineage = all_lineages,
  order = 1:length(all_cell_types)
)

# Create a flat mapping from cell type to lineage
cell_to_lineage <- data.frame(
  cell_type = unlist(lineage_to_celltypes),
  lineage = rep(names(lineage_to_celltypes), sapply(lineage_to_celltypes, length))
)

# Function to create the heatmap
create_correlation_heatmap <- function(data_file, output_file, 
                                       is_lineage = FALSE, 
                                       correlation_column = "Pearson_r",
                                       title_text = "Correlation Heatmap") {
  
  # Read the data
  data <- read_csv(data_file, show_col_types = FALSE)
  
  # Determine the label column name - in your CSVs it looks like it's "Label"
  id_column <- "Label"
  
  # Clean model names to remove the "corrected_processed_" prefix if present
  data <- data %>%
    mutate(Model = str_replace(Model, "corrected_processed_", ""))
  
  # Determine if we're using lineage or cell type data
  if (is_lineage) {
    # For lineage data, use all lineages defined in the mapping
    all_labels <- names(lineage_to_celltypes)
    # Create ordered lineage list to match the order in lineage_to_celltypes
    lineage_order <- names(lineage_to_celltypes)
  } else {
    # For cell type data, use all cell types defined in the mapping
    all_labels <- unlist(lineage_to_celltypes)
  }
  
  # Filter the data to only include labels from our reference
  data_filtered <- data %>%
    filter(!!sym(id_column) %in% all_labels)
  
  # Get unique models from the data
  models <- unique(data_filtered$Model)
  
  # Create a complete grid with all possible combinations of labels and models
  complete_grid <- expand.grid(
    Label = all_labels,
    Model = models,
    stringsAsFactors = FALSE
  )
  
  # Join the actual data with the complete grid
  plot_data <- complete_grid %>%
    left_join(data_filtered, by = c("Label", "Model"))
  
  # Print summary statistics to understand the range of correlations
  corr_summary <- summary(data_filtered[[correlation_column]])
  cat("Correlation summary statistics:\n")
  print(corr_summary)
  
  # Add a lineage column for cell types
  if (!is_lineage) {
    plot_data <- plot_data %>%
      left_join(cell_to_lineage, by = c("Label" = "cell_type")) %>%
      left_join(cell_order_df %>% select(cell_type, order), by = c("Label" = "cell_type"))
  } else {
    plot_data <- plot_data %>%
      mutate(lineage = Label) %>%
      mutate(order = match(Label, lineage_order))
  }
  
  # Order models by their mean correlation
  model_order <- data_filtered %>%
    group_by(Model) %>%
    summarize(mean_corr = mean(!!sym(correlation_column), na.rm = TRUE)) %>%
    arrange(desc(mean_corr)) %>%
    pull(Model)
  
  # Set factor levels for proper ordering
  plot_data$Model <- factor(plot_data$Model, levels = model_order)
  
  # Order lineages/cell types by their order in the original list
  if (is_lineage) {
    plot_data$Label <- factor(plot_data$Label, levels = lineage_order)
  } else {
    # Keep the original order from lineage_to_celltypes
    plot_data$Label <- factor(plot_data$Label, levels = all_cell_types)
  }
  
  # Create a custom color scale to emphasize differences when values are close to 1
  na_color <- "grey80"  # Color for NA values
  
  # Create the plot with custom color scale
  p <- ggplot(plot_data, aes(x = Model, y = Label, fill = !!sym(correlation_column))) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      name = "Correlation",
      colors = c("darkred", "red", "orange", "yellow",  "skyblue", "blue", "navy"),
      values = rescale(c(0, 0.2, 0.35, 0.5, 0.75, 0.85, 1)),
      limits = c(0, 1),
      na.value = na_color,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        barwidth = 1,
        barheight = 10
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    labs(title = title_text)
  
  # For cell types, add lineage grouping
  if (!is_lineage) {
    # Add alternating background for different lineages
    lineage_df <- data.frame(
      Label = all_cell_types,
      lineage = all_lineages
    )
    
    # Add rectangles for lineage groups
    lineage_groups <- names(lineage_to_celltypes)
    for (i in seq_along(lineage_groups)) {
      # Find the start and end indices for this lineage
      lineage_cells <- lineage_to_celltypes[[i]]
      if (length(lineage_cells) > 0) {
        min_idx <- min(which(all_cell_types %in% lineage_cells))
        max_idx <- max(which(all_cell_types %in% lineage_cells))
        
        # Add rectangle if this is an even-numbered lineage (for alternating pattern)
        if (i %% 2 == 0) {
          p <- p + annotate("rect", 
                            xmin = 0.5, xmax = length(model_order) + 0.5,
                            ymin = min_idx - 0.5, ymax = max_idx + 0.5,
                            alpha = 0.1, fill = "lightblue")
        }
        
        # Add lineage label
        p <- p + annotate("text", 
                          x = 0, y = (min_idx + max_idx) / 2,
                          label = lineage_groups[i], 
                          hjust = 1.2, size = 3)
      }
    }
    
    p <- p + scale_x_discrete(expand = expansion(mult = c(0.15, 0)))
  }
  
  # Print mean correlation per model
  model_means <- plot_data %>%
    group_by(Model) %>%
    summarize(mean_corr = mean(!!sym(correlation_column), na.rm = TRUE)) %>%
    arrange(desc(mean_corr))
  
  cat("\nModel mean correlations:\n")
  print(model_means)
  
  # Create a version with more extreme color scaling that makes differences more apparent
  p_extreme <- ggplot(plot_data, aes(x = Model, y = Label, fill = !!sym(correlation_column))) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      name = "Correlation",
      colors = c("darkred", "red", "orange", "yellow",  "skyblue", "blue", "navy"),
      values = rescale(c(0, 0.2, 0.35, 0.5, 0.75, 0.85, 1)),
      limits = c(0, 1),
      na.value = na_color,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        barwidth = 1,
        barheight = 10
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    labs(title = paste0(title_text, " (Enhanced Contrast)"))
  
  # Add the same lineage grouping to the extreme version
  if (!is_lineage) {
    # Add rectangles for lineage groups
    lineage_groups <- names(lineage_to_celltypes)
    for (i in seq_along(lineage_groups)) {
      # Find the start and end indices for this lineage
      lineage_cells <- lineage_to_celltypes[[i]]
      if (length(lineage_cells) > 0) {
        min_idx <- min(which(all_cell_types %in% lineage_cells))
        max_idx <- max(which(all_cell_types %in% lineage_cells))
        
        # Add rectangle if this is an even-numbered lineage (for alternating pattern)
        if (i %% 2 == 0) {
          p_extreme <- p_extreme + annotate("rect", 
                                            xmin = 0.5, xmax = length(model_order) + 0.5,
                                            ymin = min_idx - 0.5, ymax = max_idx + 0.5,
                                            alpha = 0.1, fill = "lightblue")
        }
        
        # Add lineage label
        p_extreme <- p_extreme + annotate("text", 
                                          x = 0, y = (min_idx + max_idx) / 2,
                                          label = lineage_groups[i], 
                                          hjust = 1.2, size = 3)
      }
    }
    
    p_extreme <- p_extreme + scale_x_discrete(expand = expansion(mult = c(0.15, 0)))
  }
  
  # Save both plots
  ggsave(output_file, p, width = 10, height = ifelse(is_lineage, 8, 16), dpi = 300)
  
  extreme_file <- gsub("\\.pdf$", "_enhanced.pdf", output_file)
  ggsave(extreme_file, p_extreme, width = 10, height = ifelse(is_lineage, 8, 16), dpi = 300)
  
  cat(paste0("Standard colorscale plot saved to: ", output_file, "\n"))
  cat(paste0("Enhanced contrast plot saved to: ", extreme_file, "\n"))
  
  return(list(standard = p, enhanced = p_extreme))
}

# Path to your CSV files
lineage_file <- "all_models_lineage_correlations.csv"
celltype_file <- "all_models_celltype_correlations.csv"

# Create the heatmaps
lineage_heatmaps <- create_correlation_heatmap(
  lineage_file,
  "lineage_correlation_heatmap.pdf",
  is_lineage = TRUE,
  correlation_column = "Pearson_r",
  title_text = "Lineage Correlation Across Models"
)

celltype_heatmaps <- create_correlation_heatmap(
  celltype_file,
  "celltype_correlation_heatmap.pdf",
  is_lineage = FALSE,
  correlation_column = "Pearson_r",
  title_text = "Cell Type Correlation Across Models"
)

# Display the enhanced contrast heatmaps in R Studio
print(lineage_heatmaps$enhanced)
print(celltype_heatmaps$enhanced)