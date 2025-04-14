library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readr)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250207_embryomodel_batch_scgpt_scarches")

# Load data
obj_new <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/R runing/human_embryo_models_filtered_v2.rds")

#########################check expression###############################

#load markers
genes_to_plot = c("CDX1","CDX2","GATA3", "TBXT", "SP5", "MIXL1","GABRP","TFAP2C",
                  "EPCAM","NANOS3", "SOX17", "PDGFRA", "PODXL", "AFP", "GJB1", "DCN"
                  
                     )


#####################plot all umaps for each gene on one plot############################

#Create a folder named "plots" if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}
# Loop through each gene
for (gene in genes_to_plot) {
  # Create an empty list to store UMAP plots for this gene
  umap_plots <- list()
  
  # Loop through each dataset
  for (table_name in names(obj_new)) {
    dataset <- obj_new[[table_name]]
    
    # Create a UMAP plot for the current gene
    plot <- FeaturePlot(dataset, features = gene, reduction = "umap", combine = FALSE)[[1]] +
      ggtitle(paste(gene, "in", table_name))  # Add title to each plot
    
    umap_plots[[table_name]] <- plot  # Store the plot in the list
  }
  
  # Combine all plots for this gene into one
  combined_plot <- gridExtra::grid.arrange(grobs = umap_plots, ncol = 2)  # Adjust ncol as needed
  
  # Save the combined plot for the gene in the "plots" folder
  ggsave(paste0("plots/", gene, "_umap_plots.pdf"), plot = combined_plot, width = 6, height = 12)  # Adjust size as needed
  
}
