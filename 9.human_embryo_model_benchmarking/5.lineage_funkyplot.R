library(funkyheatmap)
library(dplyr)
library(tibble)
library(ggplot2)

##palettes
palettes <- list(
  features = c(FULL = "black", HVG = "black"),
  black = c("black", "black"),
  orange = rev(colorRampPalette(c("#FFA13F", "white"))(9)),
  purples = rev(colorRampPalette(c("#BA7FB5", "white"))(9)),
  navys = rev(colorRampPalette(c("#748EBB", "white"))(9)),
  olives = rev(colorRampPalette(c("#5F9856", "#F7FCF5"))(9)),
  oranges2 = RColorBrewer::brewer.pal(9, "Oranges")
)

setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250119_embryo_model_benchmarking_metrics")
outputDir = getwd()


# Get all _lineage.xlsx files in the current folder
# Define the folder containing the CSV files
folder_folder <- "./output_folder"

# Get all _lineage.csv files in the output folder
file_list <- list.files(path = folder_folder, pattern = "_lineage\\.csv$", full.names = TRUE)

metric_read <- list()

# Loop through the list of files and read them into the list
for (file in file_list) {
  # Read the Excel file into a data frame
  df <- read.csv(file, row.names = 1)
  
  # Convert the data frame to a matrix and store it in the list
  metric_read[[file]] <- as.matrix(df)
}

# Remove the "./output_folder/" and "_lineage.csv" from the file names
names(metric_read) <- gsub("./output_folder/", "", file_list)
names(metric_read) <- gsub("_lineage.csv$", "", names(metric_read))


############funkyplot###################################
generate_funkyplot <- function(metrics, columns_to_rank = 2:8, legend_titles = c("certainty", "coverage", "Similarity","sum"), 
                               palette_list = c("olives","orange", "navys", "oranges2"), column_width = 1, height = 4, width = 7) {
  
  # Rank function
  label_top_3 <- function(scores) {
    ranks <- rank(-scores, ties.method = "min") 
    labels <- ifelse(ranks <= 3, as.character(ranks), "")
    return(labels)
  }
  
  # Create a plot with the rank labels
  metrics_plot <- metrics %>%
    # Create an ID column showing the final rank
    mutate(id = as.character(seq_len(nrow(metrics)))) %>%
    # Create rank labels for specified columns
    mutate(across(
      colnames(metrics)[columns_to_rank],  # Adjust column range if necessary
      label_top_3,
      .names = "cert_{.col}"  # Prefix for new column names
    )) %>%
    as.data.frame()
  
  # Define the column information with a check
  metrics_column_info <- tribble(
    ~id,                          ~id_color,                     ~name,                          ~geom,              ~group,              ~options,
    "id",                         NA,                            "Rank",                      "text",             "certainty",            list(hjust = 0, width = 1),
    "Row.names",                  NA,                            "dataset",                   "text",             "certainty",            list(hjust = 0, width = 6),
    "perct.certain_lineage",      "perct.certain_lineage",       "perct.certain_lineage",     "funkyrect",        "certainty",            list(palette = palette_list[1], width = column_width),
    "perct.certain_celltype_lineage", "perct.certain_celltype_lineage", "perct.certain_celltype_lineage", "funkyrect",  "certainty",            list(palette = palette_list[1], width = column_width),
    "coverage_lineage",           "coverage_lineage",           "coverage_lineage",          "funkyrect",        "coverage",            list(palette = palette_list[2], width = column_width),
    "coverage_celltype_lineage",  "coverage_celltype_lineage",   "coverage_celltype_lineage", "funkyrect",        "coverage",            list(palette = palette_list[2], width = column_width),
    "similarity_lineage",         "similarity_lineage",          "similarity_lineage",        "funkyrect",        "Similarity",          list(palette = palette_list[3], width = column_width),
    "similarity_celltype",        "similarity_celltype",         "similarity_celltype",       "funkyrect",        "Similarity",          list(palette = palette_list[3], width = column_width),
    "sum",                        "sum",                         "sum",                       "funkyrect",         "sum",              list(palette = palette_list[4], width = column_width)
  )
  
  ## Column Group
  metrics_column_groups <- tribble(
    ~group,         ~palette,      ~level1,
    "certainty",     palette_list[1],  "certainty",
    "coverage",      palette_list[2],  "coverage",
    "Similarity",    palette_list[3],  "Similarity",
    "sum",           palette_list[4],  "sum",
  )
  
  # Create metrics_row_info assuming the first column in metrics_plot contains identifiers
  metrics_row_info <- data.frame(id = metrics_plot$id, group = NA_character_)
  
  ## Legends
  legends <- list(
    list(
      title = legend_titles[1],
      palette = palette_list[1],
      geom = "funkyrect",
      labels = c("0", "0.25", "0.5", "0.75", "1"),
      size = seq(0, 1, length.out = 5)  # Scale from 0 to 1
    ),
    list(
      title = legend_titles[2],
      palette = palette_list[2],
      geom = "funkyrect",
      labels = c("0", "0.25", "0.5", "0.75", "1"),
      size = seq(0, 1, length.out = 5)  # Scale from 0 to 1
    ),
    list(
      title = legend_titles[3],
      palette = palette_list[3],
      geom = "funkyrect",
      labels = c("0", "0.25", "0.5", "0.75", "1"),
      size = seq(0, 1, length.out = 5)  # Scale from 0 to 1
    ),
    list(
      title = legend_titles[4],
      palette = palette_list[4],
      geom = "funkyrect",
      labels = c("0", "0.25", "0.5", "0.75", "1"),
      size = seq(0, 1, length.out = 5)  # Scale from 0 to 1
    )
  )
  
  # Create the heatmap
  p <- funky_heatmap(
    data = metrics_plot,
    column_info = metrics_column_info,
    column_groups = metrics_column_groups,  # Adjust this if needed
    row_info = metrics_row_info,
    palettes = palettes,
    legends = legends,
    position_args = position_arguments(
      col_annot_offset = 4
    ),
    scale_column = FALSE
  )
  
  # Save the plot as a PDF
  ggsave("funkyplot_output.pdf", plot = p, device = "pdf", width = width, height = height)
  
  return(p)
}


############export plot for each lineage##################################
# Loop through each dataset in metric_read
for (i in 1:length(metric_read)) {
  
  # Convert the matrix into a data frame with numeric values
  metrics <- as.data.frame(apply(metric_read[[i]], 2, as.numeric))
  
  # Retain row names
  rownames(metrics) <- rownames(metric_read[[i]])
  
  # Add Row.names as a column and reorder the columns
  metrics$Row.names <- rownames(metrics)
  metrics <- metrics[, c(7, 1:6)]
  
  # Calculate the sum of columns 2 to 7
  metrics$sum <- rowSums(metrics[, 2:7])
  
  # Replace NA values with 0 in the entire dataframe
  metrics[is.na(metrics)] <- 0
  
  
  # Normalize the metrics$sum column to the range 0 to 1
  metrics$sum <- metrics$sum / 6
  
  # Generate the plot
  p <- generate_funkyplot(
    metrics, 
    columns_to_rank = 2:8, 
    legend_titles = c("certainty", "coverage", "Similarity", "sum"), 
    palette_list = c("olives", "orange", "navys", "oranges2")
  )
  
  # Save the plot as a PDF
  ggsave(
    paste0(names(metric_read)[i], "_benchmarking.pdf"), 
    plot = p, 
    device = "pdf", 
    width = 7, 
    height = 4
  )
}
