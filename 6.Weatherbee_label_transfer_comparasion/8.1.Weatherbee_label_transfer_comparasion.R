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




# load data
obj <- readRDS("human_embryo_models_filtered_v2.rds")
obj <- obj$Weatherbee


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
lineage_ordered_labels = c("Endoderm","epi","ExE_endo","Exe_meso", "hemogenic","mesoderm", "neural_ecto", "non_neuro_ecto", "PGC", "TE_TrB",  "Gastru")


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

################################ Garfield#######################################

attri <- c("transferred_final_anno_unfiltered", "transferred_final_lineage_unfiltered")
garf <- read.csv("corrected_processed_Weatherbee.h5ad_garfield_query.obs.csv")

rownames(garf) <- garf$X

garf <- garf[,attri]


#match rownames
garf <- garf[match(rownames(obj@meta.data),rownames(garf) ),]
colnames(garf) <-paste0("human_ref_",colnames(garf))

obj$human_ref_transferred_final_anno_unfiltered <- garf$human_ref_transferred_final_anno_unfiltered
obj$human_ref_transferred_final_lineage_unfiltered <- garf$human_ref_transferred_final_lineage_unfiltered

p2 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_anno_unfiltered",cols = humanref_color_mapping) +
  ggtitle(paste0("human_ref_final_anno_unfiltered_", "Weatherbee"))

p3 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_lineage_unfiltered",cols = lineage_color_mapping) +
  ggtitle(paste0("human_ref_final_lineage_unfiltered_", "Weatherbee"))

# Save p2 as a PDF
ggsave(
  filename = "garfield_final_anno_Weatherbee.pdf",
  plot = p2,
  device = "pdf",
  width = 12,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

# Save p3 as a PDF
ggsave(
  filename = "garfield_final_lineage_Weatherbee.pdf",
  plot = p3,
  device = "pdf",
  width = 6,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)


sort(table(obj$human_ref_transferred_final_anno_unfiltered), decreasing = TRUE)

################################ scArches#######################################

attri <- c("scArches_final_anno_pre", "scArches_final_lineage_pre")
garf <- read.csv("human_Weatherbee_scArches.csv")

rownames(garf) <- garf$X

garf <- garf[,attri]


#match rownames
garf <- garf[match(rownames(obj@meta.data),rownames(garf) ),]
colnames(garf) <-paste0("human_ref_",colnames(garf))

obj$human_ref_transferred_final_anno_unfiltered <- garf$human_ref_scArches_final_anno_pre
obj$human_ref_transferred_final_lineage_unfiltered <- garf$human_ref_scArches_final_lineage_pre

p2 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_anno_unfiltered",cols = humanref_color_mapping) +
  ggtitle(paste0("human_ref_final_anno_unfiltered_", "Weatherbee"))

p3 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_lineage_unfiltered",cols = lineage_color_mapping) +
  ggtitle(paste0("human_ref_final_lineage_unfiltered_", "Weatherbee"))

# Save p2 as a PDF
ggsave(
  filename = "scArches_final_anno_Weatherbee.pdf",
  plot = p2,
  device = "pdf",
  width = 8,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

# Save p3 as a PDF
ggsave(
  filename = "scArches_final_lineage_Weatherbee.pdf",
  plot = p3,
  device = "pdf",
  width = 5,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

sort(table(obj$human_ref_transferred_final_anno_unfiltered), decreasing = TRUE)
################################ querymap #######################################

attri <- c("querymap_human_ref_final_anno", "querymap_human_ref_final_lineage")
garf <- read.csv("human_ref_weatherbee_querymap.csv")

rownames(garf) <- garf$X

garf <- garf[,attri]


#match rownames
garf <- garf[match(rownames(obj@meta.data),rownames(garf) ),]
colnames(garf) <-paste0("human_ref_",colnames(garf))

obj$human_ref_transferred_final_anno_unfiltered <- garf$human_ref_querymap_human_ref_final_anno
obj$human_ref_transferred_final_lineage_unfiltered <- garf$human_ref_querymap_human_ref_final_lineage

p2 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_anno_unfiltered",cols = humanref_color_mapping) +
  ggtitle(paste0("human_ref_final_anno_unfiltered_", "Weatherbee"))

p3 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_lineage_unfiltered",cols = lineage_color_mapping) +
  ggtitle(paste0("human_ref_final_lineage_unfiltered_", "Weatherbee"))

# Save p2 as a PDF
ggsave(
  filename = "querymap_final_anno_Weatherbee.pdf",
  plot = p2,
  device = "pdf",
  width = 8,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

# Save p3 as a PDF
ggsave(
  filename = "querymap_final_lineage_Weatherbee.pdf",
  plot = p3,
  device = "pdf",
  width = 5,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

################################ scGPT #######################################

attri <- c("scgpt_final_anno_pre", "scgpt_final_lineage_pre")
garf <- read.csv("human_Weatherbee_scgpt.csv")

rownames(garf) <- garf$X

garf <- garf[,attri]


#match rownames
garf <- garf[match(rownames(obj@meta.data),rownames(garf) ),]
colnames(garf) <-paste0("human_ref_",colnames(garf))

obj$human_ref_transferred_final_anno_unfiltered <- garf$human_ref_scgpt_final_anno_pre
obj$human_ref_transferred_final_lineage_unfiltered <- garf$human_ref_scgpt_final_lineage_pre

p2 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_anno_unfiltered",cols = humanref_color_mapping) +
  ggtitle(paste0("human_ref_final_anno_unfiltered_", "Weatherbee"))

p3 <- DimPlot(obj, reduction = "umap", group = "human_ref_transferred_final_lineage_unfiltered",cols = lineage_color_mapping) +
  ggtitle(paste0("human_ref_final_lineage_unfiltered_", "Weatherbee"))

# Save p2 as a PDF
ggsave(
  filename = "scgpt_final_anno_Weatherbee.pdf",
  plot = p2,
  device = "pdf",
  width = 8,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

# Save p3 as a PDF
ggsave(
  filename = "scgpt_final_lineage_Weatherbee.pdf",
  plot = p3,
  device = "pdf",
  width = 5,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

######################export expression#################################


p <- FeaturePlot(obj,features = c("GATA6", "SOX17", "PDGFRA", "GATA4", "PODXL", "AFP", "GJB1","DCN"
                                  ) )

p <- FeaturePlot(obj,features = c(, "COL6A1","COL6A2","ALDH1A2", "MEOX2", "FST","LUM", "POSTN", "HAND2", "ANXA1", "CREB3L1", "IGF2", "PLAGL1", "NID2", "FRMD4B", "OAF", "VCAN", "TEK", "IGFBP3",
                                  "PDPN", "COL3A1", "IGFBP7", "NPY", "DCN","HES7", "TBX6", "PRRX1","FOXF1") )

# Save p as a PDF
ggsave(
  filename = "marker expression.pdf",
  plot = p,
  device = "pdf",
  width = 12,   # Adjust the width as needed
  height = 9   # Adjust the height as needed
)







