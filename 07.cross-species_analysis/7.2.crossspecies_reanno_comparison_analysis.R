setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250725_cross-species_integration_v3")

library(Seurat)
library(tidyverse)
library(Matrix)
library(dplyr)
library(ggplot2)

outputDir = getwd()

# load data
crossspecies <- readRDS("human_monkey_crossspecies_20250725.rds")


# ============================================================================
# 1.check and correct base annotation
# ============================================================================  

human_reanno <- unique(crossspecies$reanno[which(crossspecies$species=="Homo sapiens")])
monkey_reanno <- unique(crossspecies$reanno[which(crossspecies$species=="Macaca fascicularis")])


human_reanno <- gsub("_[0-9]+$", "",human_reanno)
monkey_reanno <- gsub("_[0-9]+$", "",monkey_reanno)

setdiff(human_reanno, monkey_reanno)
setdiff(monkey_reanno,human_reanno)

# correct monkey naming
crossspecies$reanno <- gsub("Neural.ectoderm.fore_midbrain", "Neural.ectoderm.brain", crossspecies$reanno)
crossspecies$reanno[grep("Neural.ectoderm.forebrain",crossspecies$reanno) ] <- "Neural.ectoderm.brain"
crossspecies$reanno[grep("Neural.ectoderm.hindbrain",crossspecies$reanno) ] <- "Neural.ectoderm.brain"
crossspecies$reanno[grep("Neural.ectoderm.midbrain",crossspecies$reanno) ] <- "Neural.ectoderm.brain"

crossspecies$reanno[which(crossspecies$reanno=="Amniontic.epi")] <- "Amniotic.epi"
crossspecies$reanno[which(crossspecies$reanno=="VE")] <- "VE_YE"


# ============================================================================
# 2.add base annotation
# ============================================================================  

crossspecies$reanno_brief <- gsub("_[0-9]+$", "",crossspecies$reanno)

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
crossspecies_labels = c("TE", "CTB", "STB","EVT",
                        "Epiblast", "Ectoderm" ,
                        "Amniotic.epi" , "Amniontic.ectoderm",
                        "PGC" ,
                        "Primitive.streak",
                        "Neuromesodermal.progenitor",
                        "Neural.crest" ,"Neural.ectoderm.brain","Spinal.cord", 
                        "Paraxial.mesoderm", "Emergent.mesoderm", "Pre-somatic.mesoderm" , "Somite", "Rostral.mesoderm",
                        "Lateral.plate.mesoderm","Cardiac.mesoderm" ,"Connecting.stalk","Amniotic.mesoderm" ,"Exe.meso.progenitor","Allantois","Pre-YS.mesoderm", "YS.mesoderm",
                        "Hypoblast", "AVE","VE_YE","YS.endoderm" ,
                        "DE", "Gut", 
                        "Notochord" , 
                        "Hemogenic.endothelial.progenitor","Endothelium","Erythroid","Myeloid.progenitor" , "Primitive.megakaryocyte"
                                 
                        
)


crossspecies$reanno_brief <- factor(crossspecies$reanno_brief, levels = crossspecies_labels)

# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(crossspecies_labels)], crossspecies_labels)

p <- DimPlot(crossspecies, group = 'reanno_brief', reduction = "umap",, cols = color_palette)

# Save p1 (final_anno) plot as a PDF
ggsave(
  filename = paste0("crossspecies_reanno_brief.pdf"),
  plot = p,
  device = "pdf",
  width = 9,   # Adjust the width as needed
  height = 4   # Adjust the height as needed
)

# ============================================================================
# 3.calculate specie proportion by reanno_brief
# ============================================================================  

#summarize crossspecies annotations by species
meta <- crossspecies@meta.data

#define function for calculation of percentage
calculate_percentage <- function(data, group_col1, group_col2) {
  # Summarize data by the specified grouping columns
  data_sum <- data %>%
    group_by(!!sym(group_col1), !!sym(group_col2)) %>%
    summarize(frequency = n(), .groups = 'drop') %>%
    arrange(!!sym(group_col1), !!sym(group_col2))
  
  # Calculate percentage for each unique value in group_col1
  percentage <- do.call(rbind, lapply(unique(data_sum[[group_col1]]), function(x) {
    df <- data_sum[which(data_sum[[group_col1]] == x), ]
    df$perct <- df$frequency / sum(df$frequency)
    return(df)
  }))
  
  return(percentage)
}

##calculate cell_type percentage by species
percentage <- calculate_percentage(meta, "reanno_brief", "species")

percentage$reanno_brief <- factor(percentage$reanno_brief, levels = crossspecies_labels)

# Define your custom color palette
color_species <- c("Homo sapiens"="#1F77B4", "Macaca fascicularis" = "#FF7F0E")

# stacked bar plot
p <- ggplot(data=percentage, aes(x=reanno_brief, y=perct, fill=species)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_species) + 
  theme(text=element_text(size=10),axis.text.x = element_text(angle = 90, hjust = 1))  
p
# Save the plot as a PDF
ggsave("crossspeces_ref_species_summary.pdf", plot = p, device = "pdf", width = 7, height = 4)

# separate bar plot
p <- ggplot(data = percentage, aes(x = reanno_brief, y = perct)) +
  geom_bar(stat = "identity", aes(fill = species), show.legend = FALSE) +
  facet_grid(species ~ ., scales = "free_y", space = "free_y") +  # Species in rows
  ggtitle("Cell Type Composition by Species") +
  scale_fill_manual(values = color_species) + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cell Type", y = "Percentage") +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(size = 10, face = "bold", angle = 90),  # Horizontal facet labels
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

p

ggsave("crossspeces_ref_species_summary_2.pdf", plot = p, device = "pdf", width = 9, height = 6)


# ============================================================================
# 4.summarize  reanno_brief proportion by stage and dataset
# ============================================================================  

final_anno_labels <- unique(crossspecies$reanno_brief)

meta <- crossspecies@meta.data
meta <-meta[-which(meta$orig.ident %in% c("Gong","Zhai_2023", "zhou" )),]

freq_summary <- meta %>%
  group_by(stage, species, reanno_brief) %>%
  summarise(frequency = n())%>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, reanno_brief)

# ============================================================================
# 5.load monkey spatial data, summarize reanno_brief proportion by stage and dataset
# ============================================================================  

cs910 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250708_RCTD_monkey-optimization/monkey_Jia_CS910.Rds")
cs910$species <- "Macaca fascicularis"

# correct monkey naming
#cs910$final_anno <- gsub("Ectodrm", "Ectoderm", cs910$final_anno)
#cs910$final_anno <- gsub("Amnion", "Amniontic.ectoderm", cs910$final_anno)
#cs910$final_anno <- gsub("Cardiac.myocyte", "Cardiac.mesoderm", cs910$final_anno)
#cs910$final_anno <- gsub("Neural.ectoderm.fore_midbrain", "Neural.ectoderm.brain", cs910$final_anno)
#cs910$final_anno[grep("NMP",cs910$final_anno) ] <- "Neuromesodermal.progenitor"
#cs910$final_anno[grep("Axial.mesoderm",cs910$final_anno) ] <- "Paraxial.mesoderm"
#saveRDS(cs910, file ="D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250708_RCTD_monkey-optimization/monkey_Jia_CS910.Rds" )

# correct monkey naming
crossspecies$reanno <- gsub("Neural.ectoderm.fore_midbrain", "Neural.ectoderm.brain", crossspecies$reanno)
crossspecies$reanno[grep("Neural.ectoderm.forebrain",crossspecies$reanno) ] <- "Neural.ectoderm.brain"
crossspecies$reanno[grep("Neural.ectoderm.hindbrain",crossspecies$reanno) ] <- "Neural.ectoderm.brain"
crossspecies$reanno[grep("Neural.ectoderm.midbrain",crossspecies$reanno) ] <- "Neural.ectoderm.brain"

crossspecies$reanno[which(crossspecies$reanno=="Amniontic.epi")] <- "Amniotic.epi"
crossspecies$reanno[which(crossspecies$reanno=="VE")] <- "VE_YE"


cs910$reanno_brief <- gsub("_[0-9]+$", "",cs910$final_anno)

# error check 
setdiff(unique(cs910$reanno_brief), unique(crossspecies$reanno_brief))


meta_cs910 <- cs910@meta.data
colnames(meta_cs910)[which(colnames(meta_cs910)=="Carnegie_Stage")] <- "stage"

# calculate spatial summary
freq_summary_cs910 <- meta_cs910 %>%
  group_by(stage, species, reanno_brief) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, reanno_brief)


# ============================================================================
# 6.load human spatial data, summarize reanno_brief proportion by stage and dataset
# ============================================================================  

cs8 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3/human_cs8_reanno.Rds")
cs7 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3/human_cs7_reanno.Rds")

process_data <- function(data, species = "Homo sapiens") {
  # Step 1: add specie column 
  data$species <- species
  
  # Step 2: correct cell type naming
  data$final_anno <- gsub("Neural.ectoderm.fore_midbrain", "Neural.ectoderm.brain", data$final_anno)
  data$final_anno[grep("Neural.ectoderm.forebrain", data$final_anno)] <- "Neural.ectoderm.brain"
  data$final_anno[grep("Neural.ectoderm.hindbrain", data$final_anno)] <- "Neural.ectoderm.brain"
  data$final_anno[grep("Neural.ectoderm.midbrain", data$final_anno)] <- "Neural.ectoderm.brain"
  
  data$final_anno[which(data$final_anno == "Amniontic.epi")] <- "Amniotic.epi"
  data$final_anno[which(data$final_anno == "VE")] <- "VE_YE"
  
  # Step 3: add "brief" column with removing suffix
  data$reanno_brief <- gsub("_[0-9]+$", "", data$final_anno)
  
  # return data
  return(data)
}

# process data
cs7 <- process_data(cs7)
cs8 <- process_data(cs8)

# error check 
setdiff(unique(cs7$reanno_brief), unique(crossspecies$reanno_brief))
setdiff(unique(cs8$reanno_brief), unique(crossspecies$reanno_brief))

cs7$stage <- "CS7"
cs8$stage <- "CS8"

meta_cs7 <- cs7@meta.data
meta_cs8 <- cs8@meta.data

# calculate spatial summary

freq_summary_cs7 <- meta_cs7 %>%
  group_by(stage, species, reanno_brief) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, reanno_brief)

freq_summary_cs8 <- meta_cs8 %>%
  group_by(stage, species, reanno_brief) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, reanno_brief)


# ============================================================================
# 7.summarize reanno_brief distribution
# ============================================================================  

# merge summary
freq_summary <- rbind(freq_summary, freq_summary_cs7,freq_summary_cs8, freq_summary_cs910)

freq_summary$group <- paste0(freq_summary$stage, "_", freq_summary$species )

# Step 1: summarize frequency
freq_summary <- freq_summary %>%
  group_by(group, reanno_brief) %>%
  summarise(
    frequency = sum(frequency),
    .groups = 'drop'
  )

# Step 2: calculate total frequency by group
group_totals <- freq_summary %>%
  group_by(group) %>%
  summarise(total = sum(frequency), .groups = 'drop')

# Step 3: calculate percentage by group
freq_summary_pct <- freq_summary %>%
  left_join(group_totals, by = "group") %>%
  mutate(percentage = (frequency / total) * 100)

freq_summary_pct$reanno_brief <- factor(freq_summary_pct$reanno_brief, levels = crossspecies_labels)

p <- ggplot(data=freq_summary_pct, aes(x=group, y=percentage, fill=reanno_brief)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_monkey_reanno_brief_by_stage_proportion.pdf"), width = 7, height = 4)
p
dev.off()

# ============================================================================
# 8.distribution correlation matrix
# ============================================================================  

# Reshape data to wide format
wide_df <- freq_summary_pct %>%
  select(group, reanno_brief, percentage) %>%
  pivot_wider(names_from = group, values_from = percentage, values_fill = list(percentage = NA))

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

pdf(paste0(outputDir, "/crossspecies_correlation_spatial.pdf"), width = 10, height = 10)
p
dev.off()

# ============================================================================
# 9.summarize lineage distribution
# ============================================================================  

final_anno_labels <- unique(crossspecies$reanno_brief)

# define lineage colors
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

meta <- crossspecies@meta.data
meta <-meta[-which(meta$orig.ident %in% c("Gong","Zhai_2023", "zhou" )),]

freq_summary <- meta %>%
  group_by(stage, species, lineage) %>%
  summarise(frequency = n())%>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, lineage)

# assign lineage label to spatial data 
assign_lineage <- function(obj, anno_col = "reanno_brief") {
  
  # Step 1: define lineages
  lineages <- list(
    TE_TrB = c("TE", "CTB", "STB","EVT" ),
    epi = c("Epiblast", "Ectoderm"),
    Primitive.streak = c("Primitive.streak"),
    NMP = c("Neuromesodermal.progenitor"),
    Notochord = c("Notochord"),
    PGC = c("PGC"),
    ExE_endo = c("Hypoblast", "AVE","VE_YE","YS.endoderm"),
    Amniotic_ecto = c("Amniotic.epi" , "Amniontic.ectoderm"),
    neural_ecto = c("Neural.crest" ,"Neural.ectoderm.brain", "Spinal.cord"),
    Endoderm = c("DE", "Gut"),
    meso_Exe.meso = c("Paraxial.mesoderm", "Emergent.mesoderm", "Pre-somatic.mesoderm" , "Somite", "Rostral.mesoderm",
                      "Lateral.plate.mesoderm","Cardiac.mesoderm" ,"Connecting.stalk","Amniotic.mesoderm" ,"Exe.meso.progenitor","Allantois","Pre-YS.mesoderm", "YS.mesoderm"),
    hemogenic = c("Hemogenic.endothelial.progenitor","Endothelium","Erythroid","Myeloid.progenitor" , "Primitive.megakaryocyte")
  
    
  )
  
  # Step 2: add 'lineage' column with default value 'Unknown'
  obj[["lineage"]] <- "Unknown"
  
  # Step 3: add lineage label according to the anno_col
  for (lineage in names(lineages)) {
    annotations <- lineages[[lineage]]
    
    # get anno_col 
    metadata_values <- obj@meta.data[[anno_col]]
    
    # extract cell names
    cell_names <- rownames(obj@meta.data)
    
    # matching cell names
    cells_to_assign <- cell_names[metadata_values %in% annotations]
    
    # assign lineage label
    obj@meta.data[["lineage"]][which(cell_names %in%  cells_to_assign)] <- lineage
  }
  
  # return
  return(obj)
}

# human cs7, cs8

cs7 <- assign_lineage(cs7, anno_col = "reanno_brief")
cs8 <- assign_lineage(cs8, anno_col = "reanno_brief")

meta_cs7 <- cs7@meta.data
meta_cs8 <- cs8@meta.data

meta_cs7$stage <- "CS7"
meta_cs7$species <- "Homo sapiens"

meta_cs8$stage <- "CS8"
meta_cs8$species <- "Homo sapiens"


freq_summary_cs7 <- meta_cs7 %>%
  group_by(stage, species, lineage) %>%
  summarise(frequency = n())%>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, lineage)

freq_summary_cs8 <- meta_cs8 %>%
  group_by(stage, species, lineage) %>%
  summarise(frequency = n())%>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, lineage)

# monkey cs910
cs910 <- assign_lineage(cs910, anno_col = "reanno_brief")

meta_cs910 <- cs910@meta.data
colnames(meta_cs910)[which(colnames(meta_cs910)=="Carnegie_Stage")] <- "stage"

freq_summary_cs910 <- meta_cs910 %>%
  group_by(stage, species, lineage) %>%
  summarise(frequency = n()) %>%
  mutate(total = sum(frequency),  # Calculate total frequency per 'stage' and 'orig.ident'
         percentage = (frequency / total) * 100) %>%
  arrange(stage, species, lineage)

# merge summary
freq_summary <- rbind(freq_summary, freq_summary_cs7,freq_summary_cs8, freq_summary_cs910)

freq_summary$group <- paste0(freq_summary$stage, "_", freq_summary$species )

# Step 1: calculate frequency
freq_summary <- freq_summary %>%
  group_by(group, lineage) %>%
  summarise(
    frequency = sum(frequency),
    .groups = 'drop'
  )

# Step 2: calculate total frequency by group
group_totals <- freq_summary %>%
  group_by(group) %>%
  summarise(total = sum(frequency), .groups = 'drop')

# Step 3: calculate percentage by group
freq_summary_pct <- freq_summary %>%
  left_join(group_totals, by = "group") %>%
  mutate(percentage = (frequency / total) * 100)


p <- ggplot(data=freq_summary_pct, aes(x=group, y=percentage, fill=lineage)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = lineage_color) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_monkey_lineage_by_stage_proportion.pdf"), width = 7, height = 4)
p
dev.off()




