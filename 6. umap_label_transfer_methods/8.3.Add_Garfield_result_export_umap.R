
.libPaths(c("/home/liuxiaodongLab/fanxueying/R_LIBS",
            "/home/liuxiaodongLab/liuyifang/miniconda3/envs/R-4.4.0/lib/R/library"
            ))


library(devtools)
library(RColorBrewer) 
library(reshape2)
library(plotly)
library(Seurat)
library(dplyr)
library(gplots)


setwd("/storage/liuxiaodongLab/fanxueying/mayanalysis/Garfield_run/20250118_Garfield_embryo_model_batch")
outputDir = getwd()


# load data
obj <- readRDS("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/human_embryo_models_filtered.rds")

# remove Simunovic dataset
obj <- obj[!names(obj) %in% "Simunovic"]


# add garfield human_ref results
garfield_files <- list.files(outputDir, pattern = "\\.csv$", full.names = TRUE)

garfield_tables <- list()
for (file in garfield_files)  {
  
  garfield <- read.csv(file)
  
  # Extract the desired substring from the file name
  file_name <- basename(file)  # Get the file name without the path
  new_name <- sub("corrected_processed_", "", file_name)  # Remove "corrected_processed_"
  new_name <- sub(".h5ad", "", new_name)  # # Remove ".h5ad"
  new_name <- sub("(_garfield).*", "", new_name)
  
  
  
  # Store it in the list with the new name
  garfield_tables[[new_name]] <- garfield
  
}

# define threshold
threshold=0.2

#check whether names all match
all(names(garfield_tables) %in% names(obj) )


#add garfield results into obj 

attri <- c("transferred_final_anno_unfiltered","transferred_final_anno_uncert", "transferred_final_lineage_unfiltered","transferred_final_lineage_uncert")

for (i in 1:length(obj)) {
  
  dataset <- obj[[i]]
  name <- names(obj[i])
  
  garf <- garfield_tables[[name]]
  rownames(garf) <- garf$X
  
  garf <- garf[,attri]
  
  #remove filtered cells from obj.rds
  dataset <- subset(dataset, cells = rownames(garf))
  #match rownames
  garf <- garf[match(rownames(dataset@meta.data),rownames(garf) ),]
  colnames(garf) <-paste0("human_ref_",colnames(garf))
  
  meta <- cbind(dataset@meta.data,garf )
  dataset@meta.data <- meta
  
  obj[[i]] <- dataset
  
}

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

for (i in 1:length(obj)) {
  
  dataset <- obj[[i]]
  name <- names(obj[i])
  dataset <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000)
  dataset <- ScaleData(dataset)
  dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))
  dataset <- FindNeighbors(dataset, dims = 1:10) 
  dataset <- FindClusters(dataset, resolution = 0.5)
  dataset <- RunUMAP(dataset, dims = 1:10)

  obj[[i]] <- dataset
  
  # Create plots
  # Create a PDF to save all plots
  pdf(file = paste0("UMAP_all_", name, ".pdf"), width = 15, height = 7)
  
  # Create plots
  p1 <- DimPlot(dataset, reduction = "umap") +
    ggtitle(paste0("UMAP ", name))
  p2 <- DimPlot(dataset, reduction = "umap", group = "human_ref_transferred_final_anno_unfiltered",cols = humanref_color_mapping) +
    ggtitle(paste0("human_ref_transferred_reanno_unfiltered_", name))
  p3 <- DimPlot(dataset, reduction = "umap", group = "human_ref_transferred_final_lineage_unfiltered",cols = lineage_color_mapping) +
    ggtitle(paste0("human_ref_transferred_lineage_unfiltered_", name))
 
  # Print plots to PDF
  ggsave(filename = paste0("UMAP_", name, ".pdf"), plot = p1, device = "pdf")
  ggsave(filename = paste0("human_ref_transferred_reanno_unfiltered_", name, ".pdf"), plot = p2, device = "pdf")
  ggsave(filename = paste0("human_ref_transferred_lineage_unfiltered_", name, ".pdf"), plot = p3, device = "pdf")
  
}


saveRDS(obj, file = "/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/human_embryo_models_filtered_v2.rds")



