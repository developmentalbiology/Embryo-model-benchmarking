########import packages and seurat objects of dataset###########
suppressMessages(library(plotly))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(genefilter))
library(readxl)
library(devtools)
library(RColorBrewer)
library(tidyr)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250107_reference_comparasion")
outputDir = getwd()

# load crossspecies reference data
crossspecies <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/R runing/human_monkey_ref.rds")

# load human reference data
human <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/R runing/human_ref.rds")

# correct labels
human$reanno[which(human$reanno=="Primitive streak")] <- "Primitive.streak"
crossspecies$crossspecies.anno[which(crossspecies$crossspecies.anno=="Primitive streak")]<- "Primitive.streak"

#correct crossspecies$reanno
#crossspecies$reanno[which(crossspecies$species=="Homo sapiens" & crossspecies$reanno=="YS.mesoderm")] <- "YS.mesoderm_1"
#crossspecies$reanno[which(crossspecies$species=="Macaca fascicularis" & crossspecies$reanno=="Primitive.streak")] <- "Primitive streak"

# Define the ordered labels
humanref_ordered_labels <- c(
  'TE', 'CTBs', 'STBs_1','STBs_2','STBs_3', 'EVTs_1', 'EVTs_2', 
  'Epi_1','Epi_2','Epi_3','Epi_4','Primitive.streak','Notochord', 'PGC',
  'Nascent mesoderm',
  'Hypoblast', 'AVE','VE/YE','YS.Endoderm', 'Exe.meso progenitor', 'YS.mesoderm', 'Exe.endothelium',
  'Amnion', 'Amniotic_epi', 'Ectoderm','Neural crest','Neural tube_1', 'Neural tube_2','DE',
  'Paraxial mesoderm','Intermediate mesoderm','Emergent mesoderm',
  'Lateral plate mesoderm_1','Lateral plate mesoderm_2','Lateral plate mesoderm_3',
  'pre-somatic mesoderm','Somite',
  'Rostral mesoderm_1','Rostral mesoderm_2',
  'Rostral mesoderm_3', 'Cardiac myocyte',
  'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid','Myeloid progenitor'
)


# Define the ordered labels
crossspecies_ordered_labels = c('CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4',
                                'Epi_1', 'Epi_2', 'Epi_3',
                                'Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm',  'Exe.endothelium', 
                                'Amnion', 'Amniotic_epi',  'Ectoderm_1', 'Ectoderm_2',
                                'Neural tube', 'Neural crest',
                                'Primitive.streak', 'Nascent mesoderm','PGC',
                                'Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
                                'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4',
                                'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm',
                                'Cardiac myocyte', 
                                'Notochord', 'DE', 'Gut',
                                'Hypoblast', 'AVE', 'YS.Endoderm_1', 'YS.Endoderm_2', 
                                'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor'
                                
)

# scanpy color palatte
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

# correct labels
human$reanno[which(human$reanno=="Primitive streak")] <- "Primitive.streak"

df <- crossspecies@meta.data
df <- df[which(df$species=="Homo sapiens"),]
df <- df[,c("crossspecies.anno", "crossspecie.lineage")]

# error check
setdiff(rownames(human@meta.data), rownames(df))
setdiff(rownames(df),rownames(human@meta.data))

# correct cell names
diff_df_to_human <- setdiff(rownames(df), rownames(human@meta.data))

# Remove "_1" from those specific row names in df that match diff_df_to_human
rownames(df)[rownames(df) %in% diff_df_to_human] <- gsub("_1$", "", rownames(df)[rownames(df) %in% diff_df_to_human])

df <- df[match(rownames(human@meta.data), rownames(df)),]

human$crossspecies.anno <- df$crossspecies.anno
human$crossspecie.lineage <- df$crossspecie.lineage

# Retrieve unique reanno levels

human$crossspecies.anno <- factor(human$crossspecies.anno, levels = crossspecies_ordered_labels, ordered = TRUE)
reanno_levels <- levels(human$crossspecies.anno)

# Determine the number of levels in reanno
num_levels <- length(reanno_levels)

# Check if the palette has enough colors, if not, extend it
if (num_levels > length(godsnot_102)) {
  extended_palette <- colorRampPalette(godsnot_102)(num_levels)
} else {
  extended_palette <- godsnot_102[1:num_levels]
}

# Assign colors to reanno levels
assigned_colors <- setNames(extended_palette, reanno_levels)

# Run DimPlot with the assigned colors
DimPlot(human, reduction = "umap", group.by = "crossspecies.anno") +
  scale_color_manual(values = assigned_colors) +
  ggtitle("UMAP Plot with Assigned Colors") +
  theme(text = element_text(size = 16))

##################step1. check annotations between human_ref and cross-species_ref#####################

# Define function for percentage calculation
calculate_percentage <- function(data, group_col1, group_col2) {
  # Summarize data by the specified grouping columns
  data_sum <- data %>%
    group_by(!!sym(group_col1), !!sym(group_col2)) %>%
    summarize(frequency = n(), .groups = 'drop') %>%
    arrange(!!sym(group_col1), !!sym(group_col2))
  
  # Calculate percentage within each group_col1
  data_sum <- data_sum %>%
    group_by(!!sym(group_col1)) %>%
    mutate(perct = frequency / sum(frequency)) %>%
    ungroup()
  
  return(data_sum)
}

# Assuming 'meta' is already defined as human@meta.data
meta <- human@meta.data

# Calculate percentage
percentage <- calculate_percentage(meta, "reanno", "crossspecies.anno")

# Ensure the factors are ordered as desired
percentage$reanno <- factor(percentage$reanno, levels = humanref_ordered_labels, ordered = TRUE)
percentage$crossspecies.anno <- factor(percentage$crossspecies.anno, levels = crossspecies_ordered_labels, ordered = TRUE)

# Generate the plot
p <- ggplot(data=percentage, aes(x=reanno, y=perct, fill=crossspecies.anno)) +
  geom_bar(stat="identity", position="stack") +  # Using position="stack" to stack the bars
  ggtitle("reanno") +
  scale_fill_manual(values = assigned_colors) +  # Use the correctly named color palette
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)
  )  

# Print the plot
print(p)


# Calculate percentage
percentage <- calculate_percentage(meta, "crossspecies.anno","reanno")

# Ensure the factors are ordered as desired
percentage$crossspecies.anno <- factor(percentage$crossspecies.anno, levels = crossspecies_ordered_labels, ordered = TRUE)
percentage$reanno <- factor(percentage$reanno, levels = humanref_ordered_labels, ordered = TRUE)

# Retrieve unique reanno levels
human$reanno <- factor(human$reanno, levels = humanref_ordered_labels, ordered = TRUE)
reanno_levels <- levels(human$reanno)

# Determine the number of levels in reanno
num_levels <- length(reanno_levels)

# Check if the palette has enough colors, if not, extend it
if (num_levels > length(godsnot_102)) {
  extended_palette <- colorRampPalette(godsnot_102)(num_levels)
} else {
  extended_palette <- godsnot_102[1:num_levels]
}

# Assign colors to reanno levels
assigned_colors <- setNames(extended_palette, reanno_levels)

# Generate the plot
p <- ggplot(data=percentage, aes(x=crossspecies.anno, y=perct, fill=reanno)) +
  geom_bar(stat="identity", position="stack") +  # Using position="stack" to stack the bars
  ggtitle("reanno") +
  scale_fill_manual(values = assigned_colors) +  # Use the correctly named color palette
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)
  )  

# Print the plot
print(p)

##############step2. compare cell number#######################

df_count <- data.frame( table(human$lineage),table(human$crossspecie.lineage)) 
df_count <- df_count[,-3]
colnames(df_count) <- c("lineage", "human_ref", "crossspecies_ref")

# reorganize df
df_count <- df_count %>%
  pivot_longer(cols = c(human_ref, crossspecies_ref), 
               names_to = "reference", 
               values_to = "count")

# assign levels
df_count$reference <- factor(df_count$reference, levels = c("human_ref", "crossspecies_ref")) 
df_count$lineage <- factor(df_count$lineage, levels = c("TE_TrB", "epi", "Gastru", "Notochord", "PGC", "ExE_endo", "Exe_meso", "non_neuro_ecto",
                                                        "neural_ecto","Endoderm", "mesoderm","hemogenic"))  

# Create the plot
p <- ggplot(data = df_count, aes(x = lineage, y = count, fill = reference)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("Reanno") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("human_ref" = "#FFA780", "crossspecies_ref" = "#53579C")) +   # Adjust with actual group names
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.x = element_blank(),  # Optional: remove x-axis title
    axis.title.y = element_text(size = 12),  # Optional: customize y-axis title font size
    legend.position = "top"  # Optional: adjust legend position
  ) +
  scale_x_discrete(drop = FALSE) +  # Ensures missing clusters are not dropped
  scale_y_continuous(
    trans = 'log10',  # Apply log scale transformation to the y-axis
    breaks = c(1, 10, 100, 1000, 10000, 15000),  # Custom breaks (log scale)
    labels = c("1", "10", "100", "1k", "10k", "15k")  # Custom labels
  )
# Print the plot
print(p)

# Save the ggplot barplot to a PDF
ggsave("cell_count_comparasion.pdf", plot = p, width = 7, height = 6)


##############step3. compare annotations per lineage#######################

df_count <- human@meta.data[, c("lineage", "reanno")]

# Group by lineage and count unique reanno
df_count <- df_count %>%
  group_by(lineage) %>%
  summarise(unique_anno_count = n_distinct(reanno)) %>%
  arrange(lineage) 

colnames(df_count) <- c("lineage", "count")
df_count$reference <- "human_ref"

df_count2 <- human@meta.data[, c("crossspecie.lineage", "crossspecies.anno")]

# Group by final_lineage and count unique final_anno
df_count2 <- df_count2 %>%
  group_by(crossspecie.lineage) %>%
  summarise(unique_anno_count = n_distinct(crossspecies.anno)) %>%
  arrange(crossspecie.lineage) 

colnames(df_count2) <- c("lineage", "count")
df_count2$reference <- "crossspecies_ref"

df_count <- rbind(df_count, df_count2)


# assign levels
df_count$reference <- factor(df_count$reference, levels = c("human_ref", "crossspecies_ref")) 
df_count$lineage <- factor(df_count$lineage, levels = c("TE_TrB", "epi", "Gastru", "Notochord", "PGC", "ExE_endo", "Exe_meso", "non_neuro_ecto",
                                                        "neural_ecto","Endoderm", "mesoderm","hemogenic"))  

# Create the plot
p <- ggplot(data = df_count, aes(x = lineage, y = count, fill = reference)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("annotations per lineage") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("human_ref" = "#FFA780", "crossspecies_ref" = "#53579C")) +   # Adjust with actual group names
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.x = element_blank(),  # Optional: remove x-axis title
    axis.title.y = element_text(size = 12),  # Optional: customize y-axis title font size
    legend.position = "top"  # Optional: adjust legend position
  ) +
  scale_x_discrete(drop = FALSE) 
# Print the plot
print(p)

# Save the ggplot barplot to a PDF
ggsave("annotation_per_lineage_comparasion.pdf", plot = p, width = 7, height = 6)


##################step4. comparing cluster_purity#####################

# Clustering labels (annotations)  
library(bluster)


# Extract PCA embeddings, as the feature matrix
human <- NormalizeData(human, normalization.method = "LogNormalize", scale.factor = 10000)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)
human <- ScaleData(human)
human <- RunPCA(human, features = VariableFeatures(object = human))
human <- FindNeighbors(human, dims = 1:10) 
human <- FindClusters(human, resolution = 0.5)

# Function to generate mean_purity dataframe based on a dynamic clustering label
generate_mean_purity <- function(seurat_object, cluster_column = "final_lineage") {
  # Ensure the specified clustering column exists in the Seurat object
  if (!cluster_column %in% colnames(seurat_object@meta.data)) {
    stop(paste("Column", cluster_column, "not found in Seurat object's meta.data"))
  }
  
  # Extract clustering labels dynamically based on the provided column name
  clustering_labels <- as.factor(seurat_object@meta.data[[cluster_column]])
  
  # Get the PCA embeddings from the Seurat object
  human_features <- Embeddings(seurat_object, reduction = "pca")
  
  # Calculate neighbor purity (assuming neighborPurity function is defined or loaded)
  np <- neighborPurity(human_features, clustering_labels)
  
  # Convert neighbor purity result to a data frame
  pure.data <- as.data.frame(np)
  
  # Convert 'maximum' column to a factor
  pure.data$maximum <- factor(pure.data$maximum)
  
  # Add the clustering labels to the purity data frame
  pure.data$cluster <- clustering_labels
  
  # Calculate the mean purity by cluster
  mean_purity_by_cluster <- pure.data %>%
    group_by(cluster) %>%
    summarize(mean_purity = mean(purity, na.rm = TRUE))
  
  # Return the mean purity dataframe
  return(mean_purity_by_cluster)
}

mean_purity_human_anno <- generate_mean_purity(human, cluster_column = "reanno")
mean_purity_human_anno$group <- "human_reanno"

mean_purity_crossspecie_anno <- generate_mean_purity(human, cluster_column = "crossspecies.anno")
mean_purity_crossspecie_anno$group <- "crossspecie_anno"

purity <- rbind(mean_purity_human_anno,mean_purity_crossspecie_anno)

########plot#######################
# Ensure the 'group' factor has correct levels
purity$group <- factor(purity$group, levels = c("human_reanno", "crossspecie_anno"))  # Replace with your actual group names

# Ensure the 'cluster' factor has correct levels (include all unique clusters)
purity$cluster <- factor(purity$cluster, levels = unique(purity$cluster))  # Include all clusters in the factor levels

# Create the plot
p <- ggplot(data = purity, aes(x = cluster, y = mean_purity, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("Reanno") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("human_reanno" = "#FFA780", "crossspecie_anno" = "#53579C")) +   # Adjust with actual group names
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.x = element_blank(),  # Optional: remove x-axis title
    axis.title.y = element_text(size = 12),  # Optional: customize y-axis title font size
    legend.position = "top"  # Optional: adjust legend position
  ) +
  scale_x_discrete(drop = FALSE)  # Ensures missing clusters are not dropped

# Print the plot
print(p)

# Save the ggplot barplot to a PDF
ggsave("purity.pdf", plot = p, width = 12, height = 6)


##################step5. adjust for the final annotation#####################

human$final_anno <- human$crossspecies.anno
human$final_anno <- as.character(human$final_anno)

human$final_anno[which(human$reanno %in% c("TE") )] <- "TE"
human$final_anno[which(human$reanno %in% c("STBs_1") )] <- "STBs_1"
human$final_anno[which(human$reanno %in% c("STBs_2") )] <- "STBs_2"
human$final_anno[which(human$reanno %in% c("STBs_3") )] <- "STBs_3"

human$final_anno[which(human$reanno %in% c("Epi_3") )] <- "Epi_3"
human$final_anno[which(human$reanno %in% c("Epi_4") )] <- "Epi_4"

human$final_anno[which(human$reanno %in% c("Primitive.streak") )] <- "Primitive.streak"
human$final_anno[which(human$reanno %in% c("Nascent mesoderm") )] <- "Nascent mesoderm"
human$final_anno[which(human$reanno %in% c("Emergent mesoderm") )] <- "Emergent mesoderm"
human$final_anno[which(human$reanno %in% c("Intermediate mesoderm") )] <- "Intermediate mesoderm"
human$final_anno[which(human$reanno %in% c("Paraxial mesoderm") )] <- "Paraxial mesoderm"
human$final_anno[which(human$reanno %in% c("VE/YE") )] <- "VE/YE"

human$final_anno[which(human$final_anno %in% c("pre-somatic mesoderm") )] <- "NA"
human$final_anno[which(human$reanno %in% c("pre-somatic mesoderm") )] <- "pre-somatic mesoderm"
human$final_anno[which(human$final_anno %in% c("NA") )] <- "Lateral plate mesoderm_5"
human$final_anno[which(human$crossspecies.anno%in% c("Somite") )] <- "Somite"

human$final_anno[which(human$final_anno %in% c("Cardiac myocyte") )] <- "NA"
human$final_anno[which(human$reanno %in% c("Cardiac myocyte") )] <- "Cardiac myocyte"
human$final_anno[which(human$final_anno %in% c("NA") )] <- "Lateral plate mesoderm_4"


# Retrieve unique reanno levels

human$final_anno <- factor(human$final_anno)
reanno_levels <- levels(human$final_anno)

# Determine the number of levels in reanno
num_levels <- length(reanno_levels)

# Check if the palette has enough colors, if not, extend it
if (num_levels > length(godsnot_102)) {
  extended_palette <- colorRampPalette(godsnot_102)(num_levels)
} else {
  extended_palette <- godsnot_102[1:num_levels]
}

# Assign colors to reanno levels
assigned_colors <- setNames(extended_palette, reanno_levels)


DimPlot(human, group.by = "final_anno")+
  scale_color_manual(values = assigned_colors) +
  ggtitle("UMAP Plot with Assigned Colors") +
  theme(text = element_text(size = 16))


# final_anno
mean_purity_final_anno <- generate_mean_purity(human, cluster_column = "final_anno")
mean_purity_final_anno$group <- "final_anno"

purity <- rbind(purity,mean_purity_final_anno)


##################step6. plot all three#####################
# Ensure the 'group' factor has correct levels
purity$group <- factor(purity$group, levels = c("human_reanno", "crossspecie_anno", "final_anno"))  # Replace with your actual group names

# Ensure the 'cluster' factor has correct levels (include all unique clusters)
cluster_levels <- c("TE", "CTBs","CTBs_1", "CTBs_2",  "CTBs_3", "STBs_1", "STBs_2", "STBs_3", "EVTs_1", "EVTs_2", "EVTs_3", "EVTs_4",
                    "Epi_1",  "Epi_2", "Epi_3","Epi_4", 
                    "Primitive.streak","Nascent mesoderm", 
                    "Notochord", "PGC", 
                    "Hypoblast", "AVE", "VE/YE", "YS.Endoderm", "YS.Endoderm_1", "YS.Endoderm_2",
                    "Exe.meso progenitor", "Allantois_1", "Allantois_2", "pre-YS.mesoderm", "YS.mesoderm", "Exe.endothelium",
                    "Amniotic_epi", "Amnion", "Ectoderm", "Ectoderm_1", "Ectoderm_2",
                   "Neural crest",  "Neural tube", "Neural tube_1","Neural tube_2",
                   "DE", "Gut",
                   "Emergent mesoderm", "Intermediate mesoderm",  "Paraxial mesoderm", 
                   "Lateral plate mesoderm_1", "Lateral plate mesoderm_2","Lateral plate mesoderm_3", "Lateral plate mesoderm_4",
                   "Lateral plate mesoderm_5", "pre-somatic mesoderm", "Somite",  "Rostral mesoderm", "Rostral mesoderm_1", "Rostral mesoderm_2", "Rostral mesoderm_3",
                   "Cardiac myocyte", 
                   "Hemogenic endothelial progenitors", "Endothelium", "Erythroid", "Myeloid progenitor"
                   )


purity$cluster <- factor(purity$cluster, levels = cluster_levels )  # Include all clusters in the factor levels


# Create the plot
p <- ggplot(data = purity, aes(x = cluster, y = mean_purity, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("Reanno") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("human_reanno" = "#FFA780", "crossspecie_anno" = "#53579C", "final_anno"="#6DC8CF")) +   # Adjust with actual group names
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.x = element_blank(),  # Optional: remove x-axis title
    axis.title.y = element_text(size = 12),  # Optional: customize y-axis title font size
    legend.position = "top"  # Optional: adjust legend position
  ) +
  scale_x_discrete(drop = FALSE)  # Ensures missing clusters are not dropped

# Print the plot
print(p)

# Save the ggplot barplot to a PDF
ggsave("purity_human_ordered.pdf", plot = p, width = 12, height = 6)

########plot final lineage #######################

unique(human$final_anno)

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

# error check 
setdiff(final_anno_labels,unique(human$final_anno))

human$final_anno <- factor(human$final_anno, levels = final_anno_labels, ordered = TRUE)

# add final_lineage label
# Define lineages as a list of vectors
lineages <- list(
  TE_TrB = c('TE', 'CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'STBs_3','EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4'),
  epi = c('Epi_1', 'Epi_2', 'Epi_3', 'Epi_4'),
  Gastru = c('Primitive.streak', 'Nascent mesoderm'),
  Notochord = c('Notochord'),
  PGC = c('PGC'),
  ExE_endo = c('Hypoblast', 'AVE', 'VE/YE', 'YS.Endoderm_1', 'YS.Endoderm_2'),
  Exe_meso = c('Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm', 'Exe.endothelium'),
  non_neuro_ecto = c('Amnion', 'Amniotic_epi', 'Ectoderm_1', 'Ectoderm_2'),
  neural_ecto = c('Neural tube', 'Neural crest'),
  Endoderm = c('DE', 'Gut'),
  mesoderm = c('Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
               'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4', 'Lateral plate mesoderm_5',
               'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm', 'Cardiac myocyte'),
  hemogenic = c('Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor')
)

# Ensure 'lineage' column exists and initialize it to 'Unknown'
human$final_lineage <- 'Unknown'

# Loop through each lineage and assign the corresponding cells based on 'reanno'
for (lineage in names(lineages)) {
  # Get the annotations for the current lineage
  annotations <- lineages[[lineage]]
  
  # Update the 'lineage' column for cells whose 'reanno' matches any of the annotations
  human$final_lineage[human$final_anno %in% annotations] <- lineage
}

#error check
human$final_anno[which(human$final_lineage=="Unknown")]

# plot annotation

reanno_levels <- levels(human$final_anno)

# Determine the number of levels in reanno
num_levels <- length(reanno_levels)

# Check if the palette has enough colors, if not, extend it
if (num_levels > length(godsnot_102)) {
  extended_palette <- colorRampPalette(godsnot_102)(num_levels)
} else {
  extended_palette <- godsnot_102[1:num_levels]
}

# Assign colors to reanno levels
assigned_colors <- setNames(extended_palette, reanno_levels)


DimPlot(human, group.by = "final_anno")+
  scale_color_manual(values = assigned_colors) +
  ggtitle("UMAP Plot with Assigned Colors") +
  theme(text = element_text(size = 16))


mean_purity_final_lineage <- generate_mean_purity(human, cluster_column = "final_lineage")
mean_purity_final_lineage$group <- "final_lineage"

mean_purity_human_lineage <- generate_mean_purity(human, cluster_column = "lineage")
mean_purity_human_lineage$group <- "human_lineage"

mean_purity_crossspecie_lineage <- generate_mean_purity(human, cluster_column = "crossspecie.lineage")
mean_purity_crossspecie_lineage$group <- "crossspecie_lineage"

purity <- rbind(mean_purity_final_lineage,mean_purity_human_lineage,mean_purity_crossspecie_lineage)

########plot#######################
# Ensure the 'group' factor has correct levels
purity$group <- factor(purity$group, levels = c("human_lineage", "crossspecie_lineage","final_lineage"))  # Replace with your actual group names

# Ensure the 'cluster' factor has correct levels (include all unique clusters)
purity$cluster <- factor(purity$cluster, levels = unique(purity$cluster))  # Include all clusters in the factor levels

# Create the plot
p <- ggplot(data = purity, aes(x = cluster, y = mean_purity, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("Reanno") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("human_lineage" = "#FFA780", "crossspecie_lineage" = "#53579C","final_lineage"="#6DC8CF")) +   # Adjust with actual group names
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.x = element_blank(),  # Optional: remove x-axis title
    axis.title.y = element_text(size = 12),  # Optional: customize y-axis title font size
    legend.position = "top"  # Optional: adjust legend position
  ) +
  scale_x_discrete(drop = FALSE)  # Ensures missing clusters are not dropped

# Print the plot
print(p)


saveRDS(
  object = human,
  file = "human_20250108.Rds"
)

library(openxlsx)
# Export the Seurat meta.data to an Excel file
write.xlsx(human@meta.data, file = "human_20250108_metadata.xlsx")


