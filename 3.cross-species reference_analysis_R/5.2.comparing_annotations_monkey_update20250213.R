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

# subset monkey data
monkey <- subset(crossspecies, species == "Macaca fascicularis")
rm(crossspecies)

# check anno
setdiff(unique(monkey$reanno),unique(monkey$crossspecies.anno))
setdiff(unique(monkey$crossspecies.anno),unique(monkey$reanno))

# correct anno
monkey$crossspecies.anno[which(monkey$crossspecies.anno=="Primitive streak")]<- "Primitive.streak"

# get monkey umap embedding

# load monkey reference data
monkey_umap <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/R runing/monkey_ref.rds")

meta <-  monkey@meta.data

meta <- meta[match(meta$X, monkey_umap$X),]

monkey_umap$final_anno <- meta$final_anno


# Extract UMAP embeddings from monkey_umap
umap_embeddings <- Embeddings(monkey_umap, reduction = "umap")

ids <- c("TCATTACCAGAGTGTG-1_1_1_1_1_1_2", "ACTCTCGCATTACGGT-1_1_1_1_2", "AGCCAATTCGGTCACG-1_1_1_1_2",  "CTACTATAGCTGACAG-1_1_1_1_2",    
         "GGAATGGTCCACGGAC-1_1_1_1_2", "GGAGGATTCTATCCAT-1_1_1_1_2", "AGTGTTGAGCGTGAGT-1_1_1_2",   "ATAGACCAGCCGTCGT-1_1_1_2",      
         "CTGGCAGCAGACACAG-1_1_1_2",  "CTTACCGTCTCTAAGG-1_1_1_2", "TGTCCCACAATTGCTG-1_1_1_2",  "AGGCCACCAAGTTCGT-1_1_2_2",      
         "GTTCCGTAGAGAGCAA-1_1_2_2" ) 
# Remove the last occurrence of "_2" from each string
ids_cleaned <- sub("(_2)$", "", ids)

rownames(umap_embeddings)[rownames(umap_embeddings) %in% ids_cleaned] <- 
  paste0(rownames(umap_embeddings)[rownames(umap_embeddings) %in% ids_cleaned], "_2")


setdiff(rownames(umap_embeddings),monkey$X )
setdiff(monkey$X,rownames(umap_embeddings) )

umap_embeddings <- umap_embeddings[match(rownames(umap_embeddings), monkey$X),]

# Add UMAP embeddings to monkey Seurat object
monkey@reductions[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings, key = "UMAP_", assay = DefaultAssay(monkey))

saveRDS(monkey, file = "D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250107_reference_comparasion/monkey_20250108.Rds")

#########################
# plot
#scanpy color palatte
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
monkey_ordered_labels <- c(
  'TE', 'CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'STBs_3', 'STBs_4', 'STBs_5',
  'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4', 'EVTs_5',
  'Epi_1', 'Epi_2', 'Epi_3', 'Epi_4', 'Epi_5', 'Epi_6', 'Amniotic_epi', 'Amnion_1', 'Amnion_2', 'Amnion_3', 'Amniotic mesoderm',
  'Primitive.streak', 'Epithelium', 'Ectoderm_1', 'Ectoderm_2', 'PGC_1', 'PGC_2',
  'Hypoblast_1', 'Hypoblast_2', 'VE/YE', 'AVE', 'Gut_1', 'Gut_2',
  'DE', 'Notochord', 'YS.Endoderm_1', 'YS.Endoderm_2', 'YS.Endoderm_3', 'Exe.endothelium',
  'pre-Allantois', 'Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm',
  'Neural tube_1', 'Neural tube_2', 'Neural crest',
  'Nascent mesoderm', 'Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
  'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4',
  'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm',
  'Cardiac myocyte_1', 'Cardiac myocyte_2',
  'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor'
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

# error check
setdiff(unique(monkey$reanno), monkey_ordered_labels)
setdiff(unique(monkey$crossspecies.anno), crossspecies_ordered_labels)

# Retrieve unique reanno levels

monkey$reanno <- factor(monkey$reanno, levels = monkey_ordered_labels, ordered = TRUE)
reanno_levels <- levels(monkey$reanno)

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
DimPlot(monkey, reduction = "umap", group.by = "reanno") +
  scale_color_manual(values = assigned_colors) +
  ggtitle("UMAP Plot with Assigned Colors") +
  theme(text = element_text(size = 16))


# Retrieve unique reanno levels

monkey$crossspecies.anno <- factor(monkey$crossspecies.anno, levels = crossspecies_ordered_labels, ordered = TRUE)
reanno_levels <- levels(monkey$crossspecies.anno)

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
DimPlot(monkey, reduction = "umap", group.by = "crossspecies.anno") +
  scale_color_manual(values = assigned_colors) +
  ggtitle("UMAP Plot with Assigned Colors") +
  theme(text = element_text(size = 16))

##############step1. compare cell number#######################
df_count <- data.frame( table(monkey$lineage),table(monkey$crossspecie.lineage)) 
df_count <- df_count[,-3]
colnames(df_count) <- c("lineage", "monkey_ref", "crossspecies_ref")

# reorganize df
df_count <- df_count %>%
  pivot_longer(cols = c(monkey_ref, crossspecies_ref), 
               names_to = "reference", 
               values_to = "count")

# assign levels
df_count$reference <- factor(df_count$reference, levels = c("monkey_ref", "crossspecies_ref")) 
df_count$lineage <- factor(df_count$lineage, levels = c("TE_TrB", "epi", "Gastru", "Notochord", "PGC", "ExE_endo", "Exe_meso", "non_neuro_ecto",
                                                        "neural_ecto","Endoderm", "mesoderm","hemogenic"))  

# Create the plot
p <- ggplot(data = df_count, aes(x = lineage, y = count, fill = reference)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("Reanno") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("monkey_ref" = "#FFA780", "crossspecies_ref" = "#53579C")) +   # Adjust with actual group names
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
    breaks = c(1, 10, 100, 1000, 10000, 30000),  # Custom breaks (log scale)
    labels = c("1", "10", "100", "1k", "10k", "30k")  # Custom labels
  )
# Print the plot
print(p)

# Save the ggplot barplot to a PDF
ggsave("cell_count_comparasion_monkey.pdf", plot = p, width = 7, height = 6)

##############step2. compare annotations per lineage#######################
df_count <- monkey@meta.data[, c("lineage", "reanno")]

# Group by lineage and count unique reanno
df_count <- df_count %>%
  group_by(lineage) %>%
  summarise(unique_anno_count = n_distinct(reanno)) %>%
  arrange(lineage) 

colnames(df_count) <- c("lineage", "count")
df_count$reference <- "monkey_ref"

df_count2 <- monkey@meta.data[, c("crossspecie.lineage", "crossspecies.anno")]

# Group by final_lineage and count unique final_anno
df_count2 <- df_count2 %>%
  group_by(crossspecie.lineage) %>%
  summarise(unique_anno_count = n_distinct(crossspecies.anno)) %>%
  arrange(crossspecie.lineage) 

colnames(df_count2) <- c("lineage", "count")
df_count2$reference <- "crossspecies_ref"

df_count <- rbind(df_count, df_count2)


# assign levels
df_count$reference <- factor(df_count$reference, levels = c("monkey_ref", "crossspecies_ref")) 
df_count$lineage <- factor(df_count$lineage, levels = c("TE_TrB", "epi", "Gastru", "Notochord", "PGC", "ExE_endo", "Exe_meso", "non_neuro_ecto",
                                                        "neural_ecto","Endoderm", "mesoderm","hemogenic"))  

# Create the plot
p <- ggplot(data = df_count, aes(x = lineage, y = count, fill = reference)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("annotations per lineage") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("monkey_ref" = "#FFA780", "crossspecies_ref" = "#53579C")) +   # Adjust with actual group names
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
ggsave("annotation_per_lineage_comparasion_monkey.pdf", plot = p, width = 7, height = 6)



#########step3. compare cluster purity################

# Install and load required package
#if (!require("cluster")) install.packages("cluster", dependencies = TRUE)
library(cluster)

# Extract PCA embeddings (or other embeddings like UMAP)
# Use the PCA embeddings as the feature matrix
monkey <- NormalizeData(monkey, normalization.method = "LogNormalize", scale.factor = 10000)
monkey <- FindVariableFeatures(monkey, selection.method = "vst", nfeatures = 2000)
monkey <- ScaleData(monkey)
monkey <- RunPCA(monkey, features = VariableFeatures(object = monkey))
monkey <- FindNeighbors(monkey, dims = 1:10) 
monkey <- FindClusters(monkey, resolution = 0.5)


# Clustering labels (annotations)  
library(bluster)

# Function to generate mean_purity dataframe based on a dynamic clustering label
generate_mean_purity <- function(seurat_object, cluster_column = "final_lineage") {
  # Ensure the specified clustering column exists in the Seurat object
  if (!cluster_column %in% colnames(seurat_object@meta.data)) {
    stop(paste("Column", cluster_column, "not found in Seurat object's meta.data"))
  }
  
  # Extract clustering labels dynamically based on the provided column name
  clustering_labels <- as.factor(seurat_object@meta.data[[cluster_column]])
  
  # Get the PCA embeddings from the Seurat object
  features <- Embeddings(seurat_object, reduction = "pca")
  
  # Calculate neighbor purity (assuming neighborPurity function is defined or loaded)
  np <- neighborPurity(features, clustering_labels)
  
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

mean_purity_monkey <- generate_mean_purity(monkey, cluster_column = "reanno")
mean_purity_monkey$group <- "monkey"
mean_purity_crossspecies <- generate_mean_purity(monkey, cluster_column = "crossspecies.anno")
mean_purity_crossspecies$group <- "crossspecies"

purity <- rbind(mean_purity_monkey,mean_purity_crossspecies)

########plot#######################
# Ensure the 'group' factor has correct levels
purity$group <- factor(purity$group, levels = c("monkey", "crossspecies"))  # Replace with your actual group names

# Ensure the 'cluster' factor has correct levels (include all unique clusters)
purity$cluster <- factor(purity$cluster, levels = unique(purity$cluster))  # Include all clusters in the factor levels

# Create the plot
p <- ggplot(data = purity, aes(x = cluster, y = mean_purity, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("Reanno") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("monkey" = "#FFA780", "crossspecies" = "#53579C")) +   # Adjust with actual group names
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

##########step4. add final annotation######################################

monkey$final_anno <- monkey$crossspecies.anno
monkey$final_anno <- as.character(monkey$final_anno)

monkey$final_anno[which(monkey$reanno %in% c("TE") )] <- "TE"

monkey$final_anno[which(monkey$reanno %in% c("Cardiac myocyte_2") )] <- "Cardiac myocyte_2"
monkey$final_anno[which(monkey$final_anno %in% c("Cardiac myocyte") )] <- "Cardiac myocyte_1"

monkey$final_anno[which(monkey$reanno %in% c("Hypoblast_2") )] <- "Hypoblast_2"
monkey$final_anno[which(monkey$final_anno %in% c("Hypoblast") )] <- "Hypoblast_1"

monkey$final_anno[which(monkey$reanno %in% c("Gut_1") )] <- "Gut_1"
monkey$final_anno[which(monkey$final_anno %in% c("Gut") )] <- "Gut_2"

monkey$final_anno[which(monkey$reanno %in% c("Amnion_1") )] <- "Amnion_1"
monkey$final_anno[which(monkey$reanno %in% c("Amnion_2") )] <- "Amnion_2"
monkey$final_anno[which(monkey$reanno %in% c("Amnion_3") )] <- "Amnion_3"
monkey$final_anno[which(monkey$final_anno == "Amnion" )] <- "Amnion_1"

monkey$final_anno[which(monkey$reanno == "VE/YE" )] <- "VE/YE"

monkey$final_anno[which(monkey$reanno == "STBs_3" )] <- "STBs_3"
monkey$final_anno[which(monkey$reanno == "STBs_4" )] <- "STBs_4"

########plot final lineage#######################

unique(monkey$final_anno)

# Define the ordered labels
final_anno_labels = c('TE','CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'STBs_3', 'STBs_4', 'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4',
                      'Epi_1', 'Epi_2', 'Epi_3',
                      'Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm',  'Exe.endothelium', 
                      'Amniotic_epi','Amnion_1','Amnion_2','Amnion_3', 'Ectoderm_1', 'Ectoderm_2',
                      'Neural tube', 'Neural crest',
                      'Primitive.streak', 'Nascent mesoderm','PGC',
                      'Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
                      'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4',
                      'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm',
                      'Cardiac myocyte_1', 'Cardiac myocyte_2',
                      'Notochord', 'DE', 'Gut_1','Gut_2',
                      'Hypoblast_1','Hypoblast_2', 'AVE', 'VE/YE', 'YS.Endoderm_1', 'YS.Endoderm_2', 
                      'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor'
                      
)

monkey$final_anno <- factor(monkey$final_anno, levels = final_anno_labels, ordered = TRUE)

# Retrieve unique reanno levels
reanno_levels <- levels(monkey$final_anno)

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
DimPlot(monkey, reduction = "umap", group.by = "final_anno") +
  scale_color_manual(values = assigned_colors) +
  ggtitle("UMAP Plot with Assigned Colors") +
  theme(text = element_text(size = 16))

# add final_lineage label
# Define lineages as a list of vectors
lineages <- list(
  TE_TrB = c('TE', 'CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2',  'STBs_3', 'STBs_4', 'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4'),
  epi = c('Epi_1', 'Epi_2', 'Epi_3'),
  Gastru = c('Primitive.streak', 'Nascent mesoderm'),
  Notochord = c('Notochord'),
  PGC = c('PGC'),
  ExE_endo = c( 'Hypoblast_1','Hypoblast_2', 'AVE', 'VE/YE', 'YS.Endoderm_1', 'YS.Endoderm_2'),
  Exe_meso = c('Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm', 'Exe.endothelium'),
  non_neuro_ecto = c('Amniotic_epi','Amnion_1','Amnion_2','Amnion_3', 'Ectoderm_1', 'Ectoderm_2'),
  neural_ecto = c('Neural tube', 'Neural crest'),
  Endoderm = c('DE', 'Gut_1','Gut_2'),
  mesoderm = c('Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
               'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4', 'Lateral plate mesoderm_5',
               'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm', 'Cardiac myocyte_1','Cardiac myocyte_2'),
  hemogenic = c('Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor')
)

# Ensure 'lineage' column exists and initialize it to 'Unknown'
monkey$final_lineage <- 'Unknown'

# Loop through each lineage and assign the corresponding cells based on 'reanno'
for (lineage in names(lineages)) {
  # Get the annotations for the current lineage
  annotations <- lineages[[lineage]]
  
  # Update the 'lineage' column for cells whose 'reanno' matches any of the annotations
  monkey$final_lineage[monkey$final_anno %in% annotations] <- lineage
}

#error check
monkey$final_anno[which(monkey$final_lineage=="Unknown")]


saveRDS(
  object = monkey,
  file = "monkey_20250108.Rds"
)

##########step5. plot mean purity ######################################
mean_purity_final <- generate_mean_purity(monkey, cluster_column = "final_anno")
mean_purity_final$group <- "final_anno"

purity <- rbind(purity,mean_purity_final)

# Ensure the 'group' factor has correct levels
purity$group <- factor(purity$group, levels = c("monkey", "crossspecies","final_anno"))  # Replace with your actual group names

# Ensure the 'cluster' factor has correct levels (include all unique clusters)
cluster_levels <- c("TE","CTBs_1", "CTBs_2",  "CTBs_3", "STBs_1", "STBs_2", "STBs_3", "STBs_4", "EVTs_1", "EVTs_2", "EVTs_3", "EVTs_4","EVTs_5",
                    "Epi_1",  "Epi_2", "Epi_3","Epi_4", "Epi_5",
                    "Primitive.streak","Nascent mesoderm", 
                    "Notochord", "PGC", "PGC_1", "PGC_2",
                    "Hypoblast", "Hypoblast_1","Hypoblast_2", "AVE", "VE/YE", "YS.Endoderm_1", "YS.Endoderm_2","YS.Endoderm_3",
                    "pre-Allantois", "Allantois_1", "Allantois_2",   "pre-YS.mesoderm", "YS.mesoderm", "Exe.endothelium",
                    "Amniotic_epi", "Amnion", "Amnion_1", "Amnion_2","Amnion_3", "Ectoderm_1", "Ectoderm_2",
                    "Neural crest", "Neural tube", "Neural tube_1","Neural tube_2",
                    "DE","Gut","Gut_1","Gut_2",
                    "Emergent mesoderm", "Intermediate mesoderm",  "Paraxial mesoderm", 
                    "Lateral plate mesoderm_1", "Lateral plate mesoderm_2","Lateral plate mesoderm_3", "Lateral plate mesoderm_4",
                    "pre-somatic mesoderm", "Somite",  "Rostral mesoderm", "Cardiac myocyte", "Cardiac myocyte_1", "Cardiac myocyte_2", 
                    "Hemogenic endothelial progenitors", "Endothelium", "Erythroid", "Myeloid progenitor"
                    
                    
)



purity$cluster <- factor(purity$cluster, levels = cluster_levels )  # Include all clusters in the factor levels

# Create the plot
p <- ggplot(data = purity, aes(x = cluster, y = mean_purity, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +   # 'dodge' separates bars for each group
  ggtitle("Reanno") +
  # Customize color palette for the groups
  scale_fill_manual(values = c("monkey" = "#FFA780", "crossspecies" = "#53579C","final_anno"="#6DC8CF")) +   # Adjust with actual group names
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
ggsave("monkey annotation.pdf", plot = p, width = 12, height = 6)

