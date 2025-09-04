setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250709_monkey_annotation_comparasion")
outputDir = getwd()

library(Seurat)
library(tidyverse)
library(Matrix)
library(dplyr)


monkey <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/data/monkey_clustering.rds")

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
final_anno_labels = c(  'TE', 'CTB_1','CTB_2','CTB_3','CTB_4', 'STB_1','STB_2', 'EVT_1','EVT_2', 
                        'Epiblast_1','Epiblast_2', 'Ectoderm_1','Ectoderm_2',
                        'Amniotic.epi', 'Amniontic.ectoderm_1', 'Amniontic.ectoderm_2',
                        'PGC',
                        'Primitive.streak',
                        'Neuromesodermal.progenitor',
                        'Neural.crest','Neural.ectoderm.brain', 'Spinal.cord', 'Ectoderm_3','Ectoderm_4',
                        'Paraxial.mesoderm', 'Emergent.mesoderm','Pre-somatic.mesoderm', 'Somite','Rostral.mesoderm',
                        'Lateral.plate.mesoderm_1',  'Lateral.plate.mesoderm_2', 'Cardiac.mesoderm_1','Cardiac.mesoderm_2','Allantois', 
                        'Connecting.stalk','Amniotic.mesoderm', 'Exe.meso.progenitor_1', 'Exe.meso.progenitor_2',
                        'Pre-YS.mesoderm','YS.mesoderm',
                        'Hypoblast','AVE','VE_YE','YS.endoderm_1','YS.endoderm_2', 
                        'DE','Gut_1','Gut_2',
                        'Notochord_1', 'Notochord_2',
                        'Hemogenic.endothelial.progenitor', 'Endothelium', 'Erythroid','Primitive.megakaryocyte', 'Myeloid.progenitor'
                        
                        
)


# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)


# summarize sc-data cell type composition

meta <-monkey@meta.data

# Summarize the frequency and calculate percentage
freq_summary <- meta %>%
  group_by(stage, orig.ident, reanno) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, orig.ident, reanno)

# load spatial data

cs910 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250708_RCTD_monkey-optimization/monkey_Jia_CS910.Rds")

meta_cs910 <- cs910@meta.data
colnames(meta_cs910)[which(colnames(meta_cs910)=="Carnegie_Stage")] <- "stage"
meta_cs910$orig.ident <- "in house"
colnames(meta_cs910)[which(colnames(meta_cs910)=="final_anno")] <- "reanno"

# error check
setdiff(meta_cs910$reanno,final_anno_labels)
setdiff(final_anno_labels,meta_cs910$reanno)

# calculate spatial summary

freq_summary_cs910 <- meta_cs910 %>%
  group_by(stage, orig.ident, reanno) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, orig.ident, reanno)


# merge summary
freq_summary <- rbind(freq_summary, freq_summary_cs910)


# Reshape data to wide format
wide_df <- freq_summary %>%
  select(stage, orig.ident, reanno, percentage) %>%
  pivot_wider(names_from = c(stage, orig.ident), values_from = percentage, values_fill = list(percentage = NA))

# Convert all columns to numeric (except 'reanno')
wide_df_numeric <- wide_df %>%
  mutate(across(where(is.character), as.factor)) %>%  # Convert character to factor
  mutate(across(where(is.factor), as.numeric)) %>%  # Convert factor to numeric
  replace(is.na(.), 0)  # Replace NA values with 0


# Calculate the correlation matrix of the percentage values
correlation_matrix <- cor(wide_df_numeric[,-1])  # Exclude the 'reanno' column which is not numeric


library(ComplexHeatmap)
library(circlize)

# Make sure cor_matrix is a matrix (if it's not already)
cor_matrix <- as.matrix(correlation_matrix)

# Create the heatmap with dendrograms
p <- Heatmap(cor_matrix, 
             name = "Correlation",                   # Name for the color bar
             cluster_rows = TRUE,                    # Cluster rows with dendrogram
             cluster_columns = TRUE,                 # Cluster columns with dendrogram
             show_row_dend = TRUE,                   # Show row dendrogram
             show_column_dend = TRUE,                # Show column dendrogram
             col = colorRamp2(c( -0.27, 0.5, 1), c("blue", "white", "red")),  # Color scale
             heatmap_legend_param = list(title = "Correlation"),
             column_names_rot = 45,                  # Rotate column names for better visibility
             row_names_gp = gpar(fontsize = 12),     # Adjust row font size
             column_names_gp = gpar(fontsize = 12)   # Adjust column font size
)

pdf(paste0(outputDir, "/monkey_datasets_correlation_spatial.pdf"), width = 10, height = 10)
p

dev.off()





