
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

# ============================================================================
# 0. distribution of reanno by stage
# ============================================================================
df <- as.data.frame(table(human$stage, human$reanno))
colnames(df)[1:2] <- c("stage", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)
df$stage <- factor(df$stage, levels = c("E6","E6_IVC","E7_IVC","E8_IVC","E9_IVC","E10_IVC","E11_IVC","E12_IVC","E13_IVC","E14_IVC","CS7","CS9","CS10") )


df <- do.call(rbind, lapply(unique(df$stage), function(x) {
  
  df <- df[which(df$stage==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=stage, y=perct, fill=reanno)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_reanno_stage_proportion.pdf"), width = 12, height = 4)
p
dev.off()


# ============================================================================
# 1. distribution of reanno by stage
# ============================================================================

# Define the color palette
stage_color = c(
  CS7="#1f77b4",    # Blue
  E6_IVC="#aa40fc", # Purple
  CS9 = "#ff7f0e",  # Orange
  E7_IVC="#8c564b",  # Brown
  CS10="#279e68",    # Green
  E8_IVC="#e377c2",  # Pink
  E6="#d62728",      # Red
  E9_IVC= "#b5bd61",  # Yellow-Green
  E10_IVC= "#17becf", # Cyan
  E11_IVC="#aec7e8", # Light Blue
  E12_IVC="#ffbb78", # Orange-Yellow
  E13_IVC="#98df8a", # Light Green
  E14_IVC= "#ff9896"  # Light Pink
)

df <- as.data.frame(table(human$stage, human$reanno))
colnames(df)[1:2] <- c("stage", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)
df$stage <- factor(df$stage, levels = c("E6","E6_IVC","E7_IVC","E8_IVC","E9_IVC","E10_IVC","E11_IVC","E12_IVC","E13_IVC","E14_IVC","CS7","CS9","CS10") )


df <- do.call(rbind, lapply(unique(df$reanno), function(x) {
  
  df <- df[which(df$reanno==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=reanno, y=perct, fill=stage)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = stage_color) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_stage_reanno_proportion.pdf"), width = 12, height = 4)
p
dev.off()


# ============================================================================
# 2. compare orig_anno and reanno
# ============================================================================

df <- as.data.frame(table(human$orig_anno, human$reanno))
colnames(df)[1:2] <- c("orig_anno", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)


df <- do.call(rbind, lapply(unique(df$orig_anno), function(x) {
  
  df <- df[which(df$orig_anno==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=orig_anno, y=perct, fill=reanno)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/orig_anno_reanno_proportion.pdf"), width = 12, height = 4)
p
dev.off()

# ============================================================================
# 3. compare orig_sub_anno and reanno
# ============================================================================

df <- as.data.frame(table(human$orig_sub_anno, human$reanno))
colnames(df)[1:2] <- c("orig_sub_anno", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)


df <- do.call(rbind, lapply(unique(df$orig_sub_anno), function(x) {
  
  df <- df[which(df$orig_sub_anno==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=orig_sub_anno, y=perct, fill=reanno)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/orig_sub_anno_reanno_proportion.pdf"), width = 12, height = 4)
p
dev.off()

# ============================================================================
# 4. compare harmo.anno and lineage
# ============================================================================  

# Define the color palette
lineage_color = c(
  Amniotic_ecto="#1f77b4",    # Blue
  Notochord="#aa40fc", # Purple
  Endoderm = "#ff7f0e",  # Orange
  PGC="#8c564b",  # Brown
  ExE_endo="#279e68",    # Green
  Primitive.streak="#e377c2",  # Pink
  NMP="#d62728",      # Red
  TE_TrB= "#b5bd61",  # Yellow-Green
  epi= "#17becf", # Cyan
  hemogenic="#aec7e8", # Light Blue
  meso_Exe.meso="#ffbb78", # Orange-Yellow
  neural_ecto="#98df8a" # Orange-Yellow
)


df <- as.data.frame(table(human$harmo.anno, human$lineage))
colnames(df)[1:2] <- c("harmo.anno", "lineage")


df <- do.call(rbind, lapply(unique(df$harmo.anno), function(x) {
  
  df <- df[which(df$harmo.anno==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=harmo.anno, y=perct, fill=lineage)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = lineage_color) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/harmo.anno_lineage_proportion.pdf"), width = 8, height = 4)
p
dev.off()


# ============================================================================
# 5. cell type distribution by datasets
# ============================================================================  


# Define the color palette
dataset_color = c(
  Ai="#1f77b4",    # Blue
  xiang="#aa40fc", # Purple
  Yuan_2025 = "#ff7f0e",  # Orange
  zeng="#8c564b",  # Brown
  mole="#279e68",    # Green
  zhou="#e377c2",  # Pink
  tyser="#d62728"    # Red
)

df <- as.data.frame(table(human$orig.ident, human$reanno))
colnames(df)[1:2] <- c("orig.ident", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)


df <- do.call(rbind, lapply(unique(df$reanno), function(x) {
  
  df <- df[which(df$reanno==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=reanno, y=perct, fill=orig.ident)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = dataset_color) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/dataset_reanno_proportion.pdf"), width = 12, height = 4)
p
dev.off()

# ============================================================================
# 6. cell type distribution by stage and datasets
# ============================================================================  

# Get the unique 'orig.ident' values
orig_idents <- unique(human@meta.data$orig.ident)

# Define the desired overall order of stages
stage_order <- c("E6","E6_IVC","E7_IVC","E8_IVC","E9_IVC","E10_IVC","E11_IVC","E12_IVC","E13_IVC","E14_IVC","CS7","CS9","CS10")

# Directory to save the plots
output_dir <- "percentage_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# Initialize lists to store the plots

percentage_plots <- list()

# Loop through each unique 'orig.ident' and generate the plots
for (i in seq_along(orig_idents)) {
  ident <- orig_idents[i]
  
  # Subset the Seurat object's metadata for the current 'orig.ident'
  ident_data <- human@meta.data[human@meta.data$orig.ident == ident,]
  
  # Create a data frame for plotting
  ident_p1 <- as.data.frame(table(ident_data$stage, ident_data$reanno))
  colnames(ident_p1)[1:2] <- c("stage", "reanno")
  
  # Determine the available stages and reorder based on the overall order
  available_stages <- intersect(stage_order, ident_p1$stage)
  ident_p1$stage <- factor(ident_p1$stage, levels = available_stages)
  
  # Calculate percentage
  ident_p1 <- do.call(rbind, lapply(unique(ident_p1$stage), function(x) {
    df <- ident_p1[ident_p1$stage == x,]
    df$perct <- df$Freq / sum(df$Freq)
    return(df)
  }))
  
  # Generate percentage plot and store in the list
  p <- ggplot(data = ident_p1, aes(x = stage, y = perct, fill = reanno)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Percentage plot for", ident)) +
    scale_fill_manual(values = color_palette) + 
    theme(text = element_text(size = 30),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save plot to file
  file_name <- file.path(output_dir, paste0("percentage_plot_", ident, ".pdf"))
  ggsave(filename = file_name, plot = p, width = 12, height = 8, dpi = 300)
  
  # Store plot in list if needed later
  percentage_plots[[i]] <- p
}

# Optional: print all plots to console
print(percentage_plots)

# ============================================================================
# 7. correlation of cell type composition by stage and dataset
# ============================================================================ 

meta <-human@meta.data

# Summarize the frequency and calculate percentage
freq_summary <- meta %>%
  group_by(stage, orig.ident, reanno) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, orig.ident, reanno)


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
             col = colorRamp2(c( 0,0.5, 1), c("blue", "white", "red")),  # Color scale
             heatmap_legend_param = list(title = "Correlation"),
             column_names_rot = 45,                  # Rotate column names for better visibility
             row_names_gp = gpar(fontsize = 12),     # Adjust row font size
             column_names_gp = gpar(fontsize = 12)   # Adjust column font size
)

pdf(paste0(outputDir, "/human_datasets_correlation.pdf"), width = 10, height = 10)
p

dev.off()