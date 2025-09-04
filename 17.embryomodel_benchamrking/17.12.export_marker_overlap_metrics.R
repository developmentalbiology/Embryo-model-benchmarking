# Load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(viridis)
library(stringr)
library(patchwork) # For combining plots

setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250806_embryo_model_benchmarking/marker_overlap_metrics")
outputDir = getwd()


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



# Set the data directory path where your metrics files are stored
data_dir <- getwd()


# Find all the model-specific CSV files
celltype_files <- list.files(data_dir, pattern = "_celltype_marker_overlap.csv$", full.names = TRUE)
lineage_files <- list.files(data_dir, pattern = "_lineage_marker_overlap.csv$", full.names = TRUE)

cat("Found", length(celltype_files), "cell type files and", length(lineage_files), "lineage files\n")

# Print sample file names for debugging
cat("Sample cell type files:", head(basename(celltype_files)), "\n")
cat("Sample lineage files:", head(basename(lineage_files)), "\n")

# Extract model names from file names
extract_model_name <- function(file_path) {
  file_name <- basename(file_path)
  model_name <- str_replace(file_name, "_scPoli_annotated_celltype_marker_overlap.csv$|_scPoli_annotated_lineage_marker_overlap.csv$", "")
  #model_name <- str_replace(model_name, "^corrected_processed_", "")
  return(model_name)
}

# Get unique model names
model_names <- unique(c(
  sapply(celltype_files, extract_model_name),
  sapply(lineage_files, extract_model_name)
))

cat("Found", length(model_names), "unique models:", paste(model_names, collapse=", "), "\n")

# Function to process a model's CSV files and extract metrics
process_model_csvs <- function(model_name) {
  # Find CSV files for this model
  # FIXED: Use the correct file naming pattern without "_scPoli_annotated_"
  ct_file <- file.path(data_dir, paste0(model_name, "_scPoli_annotated_celltype_marker_overlap.csv"))
  lin_file <- file.path(data_dir, paste0(model_name, "_scPoli_annotated_lineage_marker_overlap.csv"))
  
  # Clean model name (remove "corrected_processed_" if present)
  clean_model_name <- str_replace(model_name, "^corrected_processed_", "")
  
  # Initialize results
  model_metrics <- list()
  
  # Process cell type metrics if file exists
  if (file.exists(ct_file)) {
    cat("Processing cell type file for model:", model_name, "\n")
    # Read CSV
    ct_data <- read_csv(ct_file, show_col_types = FALSE)
    
    # Calculate overall metrics from the CSV data
    overall_ct_metrics <- ct_data %>%
      summarize(
        model = clean_model_name,
        mean_precision = mean(Precision, na.rm = TRUE),
        weighted_precision = weighted.mean(Precision, Query_Markers, na.rm = TRUE),
        mean_recall = mean(Recall, na.rm = TRUE),
        weighted_recall = weighted.mean(Recall, Ref_Markers, na.rm = TRUE),
        mean_f1 = mean(F1_Score, na.rm = TRUE),
        mean_jaccard = mean(Jaccard, na.rm = TRUE),
        total_overlap = sum(Overlap, na.rm = TRUE),
        total_query_markers = sum(Query_Markers, na.rm = TRUE),
        total_ref_markers = sum(Ref_Markers, na.rm = TRUE),
        overlap_ratio = total_overlap / total_query_markers * 100,
        common_groups = n()
      )
    
    model_metrics$celltype <- overall_ct_metrics
    
    # Add per-group celltype data
    model_metrics$celltype_groups <- ct_data %>%
      mutate(model = clean_model_name) %>%
      select(model, Group, Precision, Recall, F1_Score, Jaccard, 
             Query_Markers, Ref_Markers, Overlap, Overlap_Percentage)
  } else {
    cat("Cell type file not found for model:", model_name, "\n")
    cat("Looked for file:", ct_file, "\n")
  }
  
  # Process lineage metrics if file exists
  if (file.exists(lin_file)) {
    cat("Processing lineage file for model:", model_name, "\n")
    # Read CSV
    lin_data <- read_csv(lin_file, show_col_types = FALSE)
    
    # Calculate overall metrics from the CSV data
    overall_lin_metrics <- lin_data %>%
      summarize(
        model = clean_model_name,
        mean_precision = mean(Precision, na.rm = TRUE),
        weighted_precision = weighted.mean(Precision, Query_Markers, na.rm = TRUE),
        mean_recall = mean(Recall, na.rm = TRUE),
        weighted_recall = weighted.mean(Recall, Ref_Markers, na.rm = TRUE),
        mean_f1 = mean(F1_Score, na.rm = TRUE),
        mean_jaccard = mean(Jaccard, na.rm = TRUE),
        total_overlap = sum(Overlap, na.rm = TRUE),
        total_query_markers = sum(Query_Markers, na.rm = TRUE),
        total_ref_markers = sum(Ref_Markers, na.rm = TRUE),
        overlap_ratio = total_overlap / total_query_markers * 100,
        common_groups = n()
      )
    
    model_metrics$lineage <- overall_lin_metrics
    
    # Add per-group lineage data
    model_metrics$lineage_groups <- lin_data %>%
      mutate(model = clean_model_name) %>%
      select(model, Group, Precision, Recall, F1_Score, Jaccard, 
             Query_Markers, Ref_Markers, Overlap, Overlap_Percentage)
  } else {
    cat("Lineage file not found for model:", model_name, "\n")
    cat("Looked for file:", lin_file, "\n")
  }
  
  return(model_metrics)
}

# Process each model
all_model_metrics <- lapply(model_names, process_model_csvs)
names(all_model_metrics) <- model_names

# Combine metrics across all models
combined_celltype_metrics <- bind_rows(lapply(all_model_metrics, function(x) {
  if (!is.null(x$celltype)) return(x$celltype) else return(NULL)
}))

combined_lineage_metrics <- bind_rows(lapply(all_model_metrics, function(x) {
  if (!is.null(x$lineage)) return(x$lineage) else return(NULL)
}))

combined_celltype_groups <- bind_rows(lapply(all_model_metrics, function(x) {
  if (!is.null(x$celltype_groups)) return(x$celltype_groups) else return(NULL)
}))

combined_lineage_groups <- bind_rows(lapply(all_model_metrics, function(x) {
  if (!is.null(x$lineage_groups)) return(x$lineage_groups) else return(NULL)
}))

# Save the combined metrics
write_csv(combined_celltype_metrics, file.path(data_dir, "combined_celltype_metrics.csv"))
write_csv(combined_lineage_metrics, file.path(data_dir, "combined_lineage_metrics.csv"))
write_csv(combined_celltype_groups, file.path(data_dir, "combined_celltype_per_group.csv"))
write_csv(combined_lineage_groups, file.path(data_dir, "combined_lineage_per_group.csv"))

cat("Saved all combined metrics files\n")

# Create plots from the combined metrics
create_metrics_plots <- function() {
  # Check if we have data
  if (nrow(combined_celltype_metrics) > 0) {
    cat("Creating cell type plots...\n")
    
    # 1. precision
    p1 <- ggplot(combined_celltype_metrics, aes(x = reorder(model, mean_precision), y = mean_precision * 100)) +
      geom_bar(stat = "identity", aes(fill = model)) +
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_precision * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker Precision by Model",
           subtitle = "Percentage of model markers found in reference",
           x = "Model", y = "Precision (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    # 2. recall
    p2 <- ggplot(combined_celltype_metrics, aes(x = reorder(model, mean_recall), y = mean_recall * 100)) +
      geom_bar(stat = "identity", aes(fill = model)) +
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_recall * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker Recall by Model",
           subtitle = "Percentage of reference markers found in model",
           x = "Model", y = "Recall (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    # 3. Cell Type Mean F1 Score
    p3 <- ggplot(combined_celltype_metrics, aes(x = reorder(model, mean_f1), y = mean_f1 * 100)) +
      geom_bar(stat = "identity", aes(fill = model))+
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_f1 * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker F1 Score by Model",
           subtitle = "Higher values indicate better marker overlap",
           x = "Model", y = "Mean F1 Score (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 4. Cell Type Mean Jaccard
    p4 <- ggplot(combined_celltype_metrics, aes(x = reorder(model, mean_jaccard), y = mean_jaccard * 100)) +
      geom_bar(stat = "identity", aes(fill = model))+
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_jaccard * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker Jaccard Index by Model",
           subtitle = "Higher values indicate better marker set similarity",
           x = "Model", y = "Mean Jaccard Index (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    # Save combined plot
    combined_ct <- (p1 + p2) / (p3 + p4)
    ggsave(file.path(data_dir, "celltype_marker_metrics.pdf"), combined_ct, width = 12, height = 10)
    
    # Individual plots
    ggsave(file.path(data_dir, "celltype_f1_score.pdf"), p1, width = 8, height = 6)
    ggsave(file.path(data_dir, "celltype_jaccard.pdf"), p2, width = 8, height = 6)
    ggsave(file.path(data_dir, "celltype_overlap_ratio.pdf"), p3, width = 8, height = 6)
    ggsave(file.path(data_dir, "celltype_marker_counts.pdf"), p4, width = 8, height = 6)
    
    # Create cell type heatmap if we have group data
    if (nrow(combined_celltype_groups) > 0) {
      # Select top models based on overall F1 score
      top_models <- combined_celltype_metrics %>%
        arrange(desc(mean_f1)) %>%
        pull(model) %>%
        head(5) 
      
      # Filter data for top models
      top_group_data <- combined_celltype_groups %>%
        filter(model %in% top_models)
      
      # Select top cell types
      top_groups <- combined_celltype_groups %>%
        group_by(Group) %>%
        summarize(mean_f1 = mean(F1_Score, na.rm = TRUE)) %>%
        arrange(desc(mean_f1)) %>%
        head(15) %>%
        pull(Group)
      
      # Filter for top groups
      plot_data <- top_group_data %>%
        filter(Group %in% top_groups)
      
      # Plot F1 scores by cell type for top models
      p5 <- ggplot(plot_data, aes(x = model, y = reorder(Group, F1_Score), fill = F1_Score)) +
        geom_tile() +
        scale_fill_viridis_c(name = "F1 Score (%)", option = "inferno") +
        labs(title = "F1 Score by Cell Type for Top Models",
             x = "Model", y = "Cell Type") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(data_dir, "celltype_f1_heatmap.pdf"), p5, width = 10, height = 8)
    }
  }
  
  if (nrow(combined_lineage_metrics) > 0) {
    cat("Creating lineage plots...\n")
    
    # 1. precision
    p1 <- ggplot(combined_lineage_metrics, aes(x = reorder(model, mean_precision), y = mean_precision * 100)) +
      geom_bar(stat = "identity", aes(fill = model)) +
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_precision * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker Precision by Model",
           subtitle = "Percentage of model markers found in reference",
           x = "Model", y = "Precision (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    # 2. recall
    p2 <- ggplot(combined_lineage_metrics, aes(x = reorder(model, mean_recall), y = mean_recall * 100)) +
      geom_bar(stat = "identity", aes(fill = model)) +
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_recall * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker Recall by Model",
           subtitle = "Percentage of reference markers found in model",
           x = "Model", y = "Recall (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    # 3. Cell Type Mean F1 Score
    p3 <- ggplot(combined_lineage_metrics, aes(x = reorder(model, mean_f1), y = mean_f1 * 100)) +
      geom_bar(stat = "identity", aes(fill = model))+
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_f1 * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker F1 Score by Model",
           subtitle = "Higher values indicate better marker overlap",
           x = "Model", y = "Mean F1 Score (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 4. Cell Type Mean Jaccard
    p4 <- ggplot(combined_lineage_metrics, aes(x = reorder(model, mean_jaccard), y = mean_jaccard * 100)) +
      geom_bar(stat = "identity", aes(fill = model))+
      scale_fill_manual(values = model_colors) +
      geom_text(aes(label = sprintf("%.1f%%", mean_jaccard * 100)), vjust = -0.5, size = 3) +
      labs(title = "Cell Type Marker Jaccard Index by Model",
           subtitle = "Higher values indicate better marker set similarity",
           x = "Model", y = "Mean Jaccard Index (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save combined plot
    combined_lin <- (p1 + p2) / (p3 + p4)
    ggsave(file.path(data_dir, "lineage_marker_metrics.pdf"), combined_lin, width = 12, height = 10)
    
    # Individual plots
    ggsave(file.path(data_dir, "lineage_Precision.pdf"), p1, width = 8, height = 6)
    ggsave(file.path(data_dir, "lineage_Recall.pdf"), p2, width = 8, height = 6)
    ggsave(file.path(data_dir, "lineage_f1_score.pdf"), p3, width = 8, height = 6)
    ggsave(file.path(data_dir, "lineage_jaccard.pdf"), p4, width = 8, height = 6)
    
    # Create lineage heatmap if we have group data
    if (nrow(combined_lineage_groups) > 0) {
      # Filter for lineage groups
      p6 <- ggplot(combined_lineage_groups, aes(x = model, y = Group, fill = F1_Score)) +
        geom_tile() +
        scale_fill_viridis_c(name = "F1 Score (%)", option = "plasma") +
        labs(title = "F1 Score by Lineage Across All Models",
             x = "Model", y = "Lineage") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(data_dir, "lineage_f1_heatmap.pdf"), p6, width = 10, height = 8)
    }
  }
}

# Create the plots
create_metrics_plots()

# Print paths to the generated files
cat("\nGenerated files:\n")
cat("1. Combined Cell Type Metrics: ", file.path(data_dir, "combined_celltype_metrics.csv"), "\n")
cat("2. Combined Lineage Metrics: ", file.path(data_dir, "combined_lineage_metrics.csv"), "\n")
cat("3. Cell Type Visualizations: ", file.path(data_dir, "celltype_marker_metrics.pdf"), "\n")
cat("4. Lineage Visualizations: ", file.path(data_dir, "lineage_marker_metrics.pdf"), "\n")



# total marker overlap counts
p1 <- ggplot(combined_lineage_metrics, aes(x = reorder(model, total_overlap), y = total_overlap)) +
  geom_bar(stat = "identity", aes(fill = model)) +
  scale_fill_manual(values = model_colors) +
  labs(title = "lineage Marker overlap by Model",
       subtitle = "number of overlapping markers",
       x = "Model", y = "counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(file.path(data_dir, "lineage_overlapping_counts.pdf"), p1, width = 8, height = 6)

# total marker overlap counts
p2 <- ggplot(combined_celltype_metrics, aes(x = reorder(model, total_overlap), y = total_overlap )) +
  geom_bar(stat = "identity", aes(fill = model)) +
  scale_fill_manual(values = model_colors) +
  labs(title = "lineage Marker overlap by Model",
       subtitle = "number of overlapping markers",
       x = "Model", y = "counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(file.path(data_dir, "celltype_overlapping_counts.pdf"), p2, width = 8, height = 6)









