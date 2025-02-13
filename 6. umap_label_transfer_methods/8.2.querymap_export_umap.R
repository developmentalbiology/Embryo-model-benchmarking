
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readr)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250205_querymap_embryomodel_batches")
outputDir = getwd()

# Load data
obj_new <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/R runing/human_embryo_models_filtered_v2.rds")

#assign colors
godsnot_102 <- c(
  "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
  "#997D87", "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF",
  "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
  "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
  "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
  "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
  "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
  "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
  "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
  "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
  "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
  "#0089A3", "#CB7E98", "#A4E804", "#324E72"
)

# Define the ordered labels

final_anno_labels = c('TE','CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'STBs_3', 'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4',
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


# Define humanref color palette
humanref_color_palette <- godsnot_102[1:length(final_anno_labels)]
humanref_color_mapping <- setNames(humanref_color_palette, final_anno_labels)


# Dynamically set levels for lineage
# Generate color palette
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
colourCount = length(lineage_ordered_labels)
lineage_color_palette <- getPalette(colourCount)
# Create a named vector for lineage colors
lineage_color_mapping <- setNames(lineage_color_palette, lineage_ordered_labels)

################################ querymap #######################################

# Get the list of all CSV files in the directory
csv_files <- list.files(path = outputDir, pattern = "_human_ref_querymap.csv", full.names = TRUE)

# Read all CSV files into a list
csv_list <- lapply(csv_files, read.csv)

# Optionally, assign names to the list based on the CSV file names
names(csv_list) <- gsub("_human_ref_querymap.csv", "", basename(csv_files))

# add querymap results to the rds

attri <- c("querymap_human_ref_final_anno", "querymap_human_ref_final_lineage")


# Loop through each CSV in csv_list and process the corresponding Seurat object
for (i in 1:length(csv_list)) {
  
  # Extract data from the i-th CSV in csv_list
  data <- csv_list[[i]]
  
  # Extract the corresponding Seurat object from obj_new using the name of the i-th CSV file
  obj_name <- names(csv_list)[i]  # Get the name of the i-th CSV file
  obj <- obj_new[[obj_name]]       # Extract the corresponding Seurat object from obj_new
  
  # Set rownames for data and select required attributes
  rownames(data) <- data$X
  data <- data[, attri]  # Ensure 'attri' is defined in your code
  
  # Match row names between data and Seurat object meta data
  data <- data[match(rownames(obj@meta.data), rownames(data)), ]
  
  # Add the querymap annotations to the Seurat object
  obj$querymap_human_ref_final_anno <- data$querymap_human_ref_final_anno
  obj$querymap_human_ref_final_lineage <- data$querymap_human_ref_final_lineage
  
  # Create UMAP plots for the annotations
  p1 <- DimPlot(obj, reduction = "umap", group = "querymap_human_ref_final_anno", cols = humanref_color_mapping) +
    ggtitle(paste0("querymap_human_ref_final_anno_", obj_name))
  
  p2 <- DimPlot(obj, reduction = "umap", group = "querymap_human_ref_final_lineage", cols = lineage_color_mapping) +
    ggtitle(paste0("querymap_human_ref_final_lineage_", obj_name))
  
  # Save p1 (final_anno) plot as a PDF
  ggsave(
    filename = paste0("querymap_human_ref_final_anno_", obj_name, ".pdf"),
    plot = p1,
    device = "pdf",
    width = 8,   # Adjust the width as needed
    height = 4   # Adjust the height as needed
  )
  
  # Save p2 (final_lineage) plot as a PDF
  ggsave(
    filename = paste0("querymap_human_ref_final_lineage_", obj_name, ".pdf"),
    plot = p2,
    device = "pdf",
    width = 5,   # Adjust the width as needed
    height = 4   # Adjust the height as needed
  )
}
