
setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_annotation_comparasion")
outputDir = getwd()

library(Seurat)
library(tidyverse)
library(Matrix)
library(dplyr)


human <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/data/human_clustering_20250724_v3.rds")

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
final_anno_labels = c('TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
                      'Epiblast_1','Epiblast_2','Epiblast_3','Ectoderm',
                      'Amniontic.epi','Amniontic.ectoderm',
                      'PGC',
                      'Primitive.streak',
                      'Neuromesodermal.progenitor',
                      'Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 'Neural.ectoderm.midbrain','Spinal.cord',
                      'Paraxial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 'Rostral.mesoderm', 'Lateral.plate.mesoderm_1',
                      'Lateral.plate.mesoderm_2','Lateral.plate.mesoderm_3','Cardiac.mesoderm','Amniotic.mesoderm','Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2',
                      'Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm',
                      'DE','Gut',
                      'Notochord',
                      'Hemogenic.endothelial.progenitor','Endothelium','Erythroid','Primitive.megakaryocyte','Myeloid.progenitor'
)


# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)

## correct "Primitive Streak" in orig_anno and orig_sub_anno
human$orig_anno[which(human$orig_anno=="Primitive Streak")] <- "Primitive.streak"  
human$orig_sub_anno[which(human$orig_sub_anno=="Primitive Streak")] <- "Primitive.streak"  

# summarize sc-data cell type composition

meta <-human@meta.data

# Summarize the frequency and calculate percentage
freq_summary <- meta %>%
  group_by(stage, orig.ident, reanno) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, orig.ident, reanno)

# load spatial data

cs8 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3/human_cs8_reanno.Rds")
cs7 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3/human_cs7_reanno.Rds")

meta_cs7 <- cs7@meta.data
meta_cs7$stage <- "CS7"
meta_cs7$orig.ident <- "Cui_2025"
colnames(meta_cs7)[which(colnames(meta_cs7)=="final_anno")] <- "reanno"

meta_cs8 <- cs8@meta.data
meta_cs8 <- cs8@meta.data
meta_cs8$stage <- "CS8"
meta_cs8$orig.ident <- "Xiao_2024"
colnames(meta_cs8)[which(colnames(meta_cs8)=="final_anno")] <- "reanno"

# calculate spatial summary

freq_summary_cs7 <- meta_cs7 %>%
  group_by(stage, orig.ident, reanno) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, orig.ident, reanno)

freq_summary_cs8 <- meta_cs8 %>%
  group_by(stage, orig.ident, reanno) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, orig.ident, reanno)

# merge summary
freq_summary <- rbind(freq_summary, freq_summary_cs7,freq_summary_cs8)


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
             col = colorRamp2(c( -0.25, 0.5, 1), c("blue", "white", "red")),  # Color scale
             heatmap_legend_param = list(title = "Correlation"),
             column_names_rot = 45,                  # Rotate column names for better visibility
             row_names_gp = gpar(fontsize = 12),     # Adjust row font size
             column_names_gp = gpar(fontsize = 12)   # Adjust column font size
)

pdf(paste0(outputDir, "/human_datasets_correlation_spatial.pdf"), width = 10, height = 10)
p

dev.off()






