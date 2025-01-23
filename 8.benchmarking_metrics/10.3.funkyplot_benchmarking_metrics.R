library(funkyheatmap)
library(dplyr)
library(tibble)
library(ggplot2)


##palettes
palettes <- list(
  features = c(FULL = "black", HVG = "black"),
  black = c("black", "black"),
  reds = rev(colorRampPalette(c("#FF3333", "white"))(9)),
  greens = rev(colorRampPalette(c("#32A02C", "white"))(9)),
  blues = rev(colorRampPalette(c("#80B1D2", "white"))(9)),
  yellows = rev(colorRampPalette(c("#FAF5B5", "white"))(9)),
  oranges = rev(colorRampPalette(c("#F16913", "#FFF5EB"))(9)),
  purples = rev(colorRampPalette(c("#BA7FB5", "white"))(9)),
  pinks = rev(colorRampPalette(c("#F7CBDF", "white"))(9)),
  browns = rev(colorRampPalette(c("#CA9C91", "white"))(9)),
  greys = rev(colorRampPalette(c("#808080", "white"))(9)),
  cyans = rev(colorRampPalette(c("#8CD0C3", "white"))(9)),
  magentas = rev(colorRampPalette(c("#BCB9D8", "white"))(9)),
  limes = rev(colorRampPalette(c("#B2E08A", "white"))(9)),
  navys = rev(colorRampPalette(c("#748EBB", "white"))(9)),
  olives = rev(colorRampPalette(c("#5F9856", "#F7FCF5"))(9)),
  maroons = rev(colorRampPalette(c("#800000", "white"))(9)),
  teals = rev(colorRampPalette(c("#008080", "white"))(9)),
  orange = rev(colorRampPalette(c("#FFA13F", "white"))(9)),
  golds = rev(colorRampPalette(c("#FFD700", "white"))(9)),
  oranges2 = RColorBrewer::brewer.pal(9, "Oranges"),
  purples2 = RColorBrewer::brewer.pal(9, "Purples"),
  greens2 = rev(RColorBrewer::brewer.pal(9, "Greens")),
  YlOrBr =rev(RColorBrewer::brewer.pal(9, "YlOrBr")),
  greys = RColorBrewer::brewer.pal(9, "Greys"),
  blus2 = RColorBrewer::brewer.pal(9, "Blues"),
  RdPu = RColorBrewer::brewer.pal(9, "RdPu")
)


metrics <- read.csv("embryomodel_metrics.csv", row.names = 1)
metrics$Row.names <- rownames(metrics)
metrics <- metrics[,c(7,1:6)]


############funkyplot###################################

label_top_3 <- function(scores) {
  ranks <- rank(-scores, ties.method = "min") 
  labels <- ifelse(ranks <= 3, as.character(ranks), "")
  return(labels)
}

metrics_plot <- metrics |>
  # Create an ID column showing the final rank
  mutate(id = as.character(seq_len(nrow(metrics)))) |>
  # Create rank labels for specified columns
  mutate(across(
    c(certain_lineage, certain_celltype,
      cov_lineage, cov_celltype, correlation_human_ref_lineage, correlation_human_ref_celltype),
    label_top_3,
    .names = "cert_{.col}"  # Prefix for new column names
  )) |>
  as.data.frame()

# Define the column information with a check
metrics_column_info <- tribble(
  ~id,                          ~id_color,                     ~name,                          ~geom,              ~group,              ~options,
  "id",                         NA,                            "Rank",                      "text",             "certainty",            list(hjust = 0,width = 1),
  "Row.names",                     NA,                            "dataset",                    "text",             "certainty",            list(hjust = 0, width = 6),
  
  "certain_lineage",            "certain_lineage",            "certain_lineage",        "funkyrect",           "certainty",       list(palette = "olives", width = 1),
  "certain_celltype",          "certain_celltype",           "certain_celltype",         "funkyrect",         "certainty",       list(palette = "olives", width = 1),
  
  "cov_lineage",                 "cov_lineage",            "cov_lineage",           "funkyrect",              "coverage",           list(palette = "oranges2", width = 1),
  "cov_celltype",            "cov_celltype",           "cov_celltype",          "funkyrect",              "coverage",           list(palette = "oranges2", width = 1),
  
  "correlation_human_ref_lineage",   "correlation_human_ref_lineage",    "correlation_human_ref_lineage",    "funkyrect",     "Similarity",         list(palette = "navys", width = 1),
  "correlation_human_ref_celltype",   "correlation_human_ref_celltype",    "correlation_human_ref_celltype",  "funkyrect",   "Similarity",         list(palette = "navys", width = 1),
  
)

##column_group
metrics_column_groups <- tribble(
  ~group,         ~palette,      ~level1,
  "certainty",     "olives",      "certainty",
  "coverage",     "oranges2",      "coverage",
  "Similarity",     "navys",       "Similarity"
  
)

metrics_column_groups

# Create metrics_row_info assuming the first column in metrics_plot contains identifiers
metrics_row_info <- data.frame(id = metrics_plot$id, group = NA_character_)


##legends
legends <- list(
  
  list(
    title = "certainty",
    palette = "olives",
    geom = "funkyrect",
    labels = c("0", " ", "0.5", " ", "1"),
    size = c(1, 1, 1, 1, 1)
  ),
  list(
    title = "coverage",
    palette = "oranges2", 
    geom = "rect",
    labels = c("0", " ", "0.5", " ", "1"),
    size = c(1, 1, 1, 1, 1)
  ),
  list(
    title = "Similarity",
    palette = "navys", 
    geom = "rect",
    labels = c("0", " ", "0.5", " ", "1"),
    size = c(1, 1, 1, 1, 1)
  )
)



# Create the heatmap
p <- funky_heatmap(
  data = metrics_plot,
  column_info = metrics_column_info,
  column_groups = metrics_column_groups,  # You can adjust this if needed
  #row_info = metrics_row_info,            # Adjust row info accordingly
  palettes = palettes,
  legends = legends,
  position_args = position_arguments(
    col_annot_offset =4
  ),
  scale_column = FALSE
)

p
# Save the plot as a PDF
ggsave("embryo_model_benchmarking_metrics.pdf", plot = p, device = "pdf", width = 7, height = 4)



