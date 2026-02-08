setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250725_cross-species_integration_v3")

library(Seurat)
library(tidyverse)
library(Matrix)
library(dplyr)
library(ggplot2)

outputDir = getwd()

# load data
data <- read.csv("human_monkey_correlation_jaccard.csv")

process_data <- function(data) {
  
  # define a function to correct naming
  clean_celltype <- function(x) {
    # 1. correct celltype naming
    x <- gsub("Ectodrm", "Ectoderm", x)
    x <- gsub("Amniontic.ectoderm", "Amnion", x)
    x <- gsub("Cardiac.myocyte", "Cardiac.mesoderm", x)
    x <- gsub("Neural.ectoderm.fore_midbrain", "Neural.ectoderm.brain", x)
    x[grep("Neuromesodermal.progenitor", x)] <- "NMP"
    
    x[grep("Neural.ectoderm.forebrain", x)] <- "Neural.ectoderm.brain"
    x[grep("Neural.ectoderm.hindbrain", x)] <- "Neural.ectoderm.brain"
    x[grep("Neural.ectoderm.midbrain", x)] <- "Neural.ectoderm.brain"
    
    x[which(x == "Amniontic.epi")] <- "Amniotic.epi"
    x[which(x == "VE")] <- "VE/YE"
    x[which(x == "Myocyte.progenitor")] <- "Rostral.mesoderm"
    
    #  correct naming specifically for human human_celltype
    if ("human_celltype" %in% names(data)) {
      x[which(x == "Connecting.stalk" & data$species == "Homo sapiens")] <- "Lateral.plate.mesoderm_3"
    }
    
    # 2. remove suffix（such as STB_1 -> STB）
    x <- gsub("_(\\d+)$", "", x)
    
    return(x)
  }
  
  # create new human_brief and monkey_brief
  if ("human_celltype" %in% names(data)) {
    data$human_brief <- clean_celltype(data$human_celltype)
  }
  
  if ("monkey_celltype" %in% names(data)) {
    data$monkey_brief <- clean_celltype(data$monkey_celltype)
  }
  
  return(data)
}


# process data
data <- process_data(data)



# 1. draw heatmap

# use tidyr::pivot_wider to convert data format
library(tidyr)
wide_data <- data %>%
  pivot_wider(
    names_from = monkey_celltype,
    values_from = c(correlation, jaccard_index),
    id_cols = human_celltype
  )

# convert wide_data to long format
long_data <- wide_data %>%
  pivot_longer(
    cols = -human_celltype,
    names_to = c(".value", "monkey_celltype"),
    names_pattern = "(correlation|jaccard_index)_(.+)"
  )

# set ggplot theme
theme_set(theme_minimal())

# 4. generate plot
library(viridis)  
ggplot(long_data, aes(x = monkey_celltype, y = human_celltype)) +
  geom_point(
    aes(size = jaccard_index, color = correlation),
    shape = 16,  
    stroke = 0   
  ) +
  scale_color_viridis( name = "Correlation") +  #
  scale_size_continuous(range = c(2, 15), name = "Jaccard Index") +  
  labs(
    title = "Heatmap of Correlation and Jaccard Index",
    x = "Monkey Cell Type",
    y = "Human Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


# 2. draw dotplot

library(ggplot2)
library(dplyr)

# Step 1: define colors
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
  "#772600", "#D790FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC",
  "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72"
)

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

# assign colors
color_palette <- setNames(godsnot_102[seq_along(crossspecies_labels)], crossspecies_labels)

# Step 2: only keep rows which human_brief == monkey_brief 
matched_df <- data[data$human_brief == data$monkey_brief, , drop = FALSE]
rownames(matched_df) <- NULL  

# Step 3: plot
p <- ggplot(matched_df, aes(x = correlation, y = jaccard_index, color = human_brief)) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  labs(
    title = "Correlation vs Jaccard Index (Human Brief = Monkey Brief)",
    x = "Pearson Correlation",
    y = "Jaccard Index",
    color = "Cell Type"
  ) +
  xlim(0.4, 0.76) +
  ylim(0, 0.4) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(
    legend.position = c(1, 1),  
    legend.justification = c(1, 1),
    legend.title = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),  
    legend.key.size = unit(6, "points"),  
    legend.text = element_text(size = 8)  
  )



print(p)

# Step 4: save the plot
ggsave("correlation_vs_jaccard.pdf", plot = p, width = 13, height = 4)












