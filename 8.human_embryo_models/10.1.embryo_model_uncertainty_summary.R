library(readxl)
library(devtools)
library(RColorBrewer) 
library(reshape2)
library(plotly)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250118_Garfield_embryo_model_batch")
outputDir = getwd()


# List all .csv files in the directory
folderpath <- "./data"
garfield_files <- list.files(folderpath, pattern = "\\.csv$", full.names = TRUE)

garfield_tables <- list()
for (file in garfield_files)  {
  
  garfield <- read.csv(file)
  
  # Extract the desired substring from the file name
  file_name <- basename(file)  # Get the file name without the path
  new_name <- sub("corrected_processed_", "", file_name)  # Remove "corrected_processed_"
  new_name <- sub(".h5ad", "", new_name)  # Remove "corrected_processed_"
  new_name <- sub("(_garfield).*", "", new_name)
  
  
  
  # Store it in the list with the new name
  garfield_tables[[new_name]] <- garfield
  
}


# check lineage uncertainty
# Loop through each table in the list
for (table_name in names(garfield_tables)) {
  garfield <- garfield_tables[[table_name]]
  
  # Check if 'transferred_lineage_uncert' exists in the table
  if ("transferred_final_lineage_uncert" %in% colnames(garfield)) {
    min_value <- min(garfield$transferred_final_lineage_uncert, na.rm = TRUE)
    max_value <- max(garfield$transferred_final_lineage_uncert, na.rm = TRUE)
    
    # Print the results
    cat("Table:", table_name, "Min:", min_value, "Max:", max_value, "\n")
  } else {
    cat("Table:", table_name, "does not have 'transferred_lineage_uncert' column.\n")
  }
}


# check reanno uncertainty
# Loop through each table in the list
for (table_name in names(garfield_tables)) {
  garfield <- garfield_tables[[table_name]]
  
  # Check if 'transferred_reanno_uncert' exists in the table
  if ("transferred_final_anno_uncert" %in% colnames(garfield)) {
    min_value <- min(garfield$transferred_final_anno_uncert, na.rm = TRUE)
    max_value <- max(garfield$transferred_final_anno_uncert, na.rm = TRUE)
    
    # Print the results
    cat("Table:", table_name, "Min:", min_value, "Max:", max_value, "\n")
  } else {
    cat("Table:", table_name, "does not have 'transferred_reanno_uncert' column.\n")
  }
}

######plot percentage##############

# define threshold
threshold=0.2

perct_results <- list()

# Loop through each table in the garfield_tables list
for (table_name in names(garfield_tables)) {
  garfield <- garfield_tables[[table_name]]
  
  # Assign groups based on the threshold
  garfield$group <- "uncertain"
  garfield$group[garfield$transferred_final_lineage_uncert < threshold] <- "certain lineage"
  garfield$group[garfield$transferred_final_anno_uncert < threshold] <- "certain cell type"
  
  # Create a data frame of counts
  df <- as.data.frame(table(garfield$group))
  colnames(df) <- c("group", "count")
  df$perct <- df$count / sum(df$count)
  df$dataset <- table_name  # Use the table name as the dataset identifier
  
  perct_results[[table_name]] <- df  # Store the result in the list
}

# Combine all results into a single data frame
perct_results <- do.call(rbind, perct_results)

perct_results$group <- factor(perct_results$group, levels=c("uncertain","certain lineage","certain cell type"))


# Generate the bar plot
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
colors <- getPalette(length(unique(perct_results$group) ))
# Assign colors to each group
colors <- c("uncertain" = colors[3], 
            "certain lineage" = colors[2], 
            "certain cell type" = colors[1])

p <-  ggplot(data=perct_results, aes(x=dataset, y=perct, fill=group)) +
  geom_bar(stat="identity", position="fill") +  # Use position="fill" for percentage stacking
  ggtitle("Prediction Summary") +
  scale_fill_manual(values=colors) + 
  theme(text=element_text(size=10))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

p

# Save the plot as a PDF
ggsave("Prediction Summary_human_ref.pdf", plot = p, device = "pdf", width = 6, height = 4)



######plot cell proportion per dataset with certain lineage##############

proportion_results <- list()
for (table_name in names(garfield_tables)) {
  
  garfield <- garfield_tables[[table_name]]
  
  # Assign groups based on the threshold
  garfield$group <- "uncertain"
  garfield$group[garfield$transferred_final_lineage_uncert < threshold] <- "certain lineage"
  
  df <- as.data.frame(table(garfield$transferred_final_lineage_unfiltered))
  colnames(df) <- c("lineage", "count")
  df$perct <- df$count / sum(df$count)
  df$dataset <- table_name
  
  proportion_results[[table_name]] <- df  # Store the result in the list
  
}

# Combine all results into a single data frame
proportion_results <- do.call(rbind, proportion_results)

proportion_results$lineage <- factor(proportion_results$lineage, levels = c("Endoderm","epi","ExE_endo","Exe_meso", "hemogenic","mesoderm", "neural_ecto", "non_neuro_ecto","Notochord",  "PGC",  "TE_TrB", "Gastru"), ordered = TRUE)

getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
colourCount = length(unique(proportion_results$lineage))

p <- ggplot(data=proportion_results, aes(x=dataset, y=perct, fill=lineage)) +
  geom_bar(stat="identity", position="fill") +  # Use position="fill" for percentage stacking
  ggtitle("certain lineage") +
  scale_fill_manual(values=getPalette(colourCount)) + 
  theme(text=element_text(size=10)) 


p

# Save the plot as a PDF
ggsave("embryo_model_certain_lineage.pdf", plot = p, device = "pdf", width = 6, height = 4)


######plot cell proportion per dataset with certain cell type##############
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

# Define humanref color palette
humanref_color_palette <- godsnot_102[1:length(final_anno_labels)]
humanref_color_mapping <- setNames(humanref_color_palette, final_anno_labels)

proportion_results <- list()
for (table_name in names(garfield_tables)) {
  
  garfield <- garfield_tables[[table_name]]
  
  # Assign groups based on the threshold
  garfield$group <- "uncertain"
  garfield$group[garfield$transferred_final_anno_uncert < threshold] <- "certain cell type"
  
  garfield <- garfield[which(garfield$group=="certain cell type"),]
  
  df <- as.data.frame(table(garfield$transferred_final_anno_unfiltered))
  colnames(df) <- c("reanno", "count")
  df$perct <- df$count / sum(df$count)
  df$dataset <- table_name
  
  proportion_results[[table_name]] <- df  # Store the result in the list
  
}

# Combine all results into a single data frame
proportion_results <- do.call(rbind, proportion_results)

#error check
all(final_anno_labels %in% unique(proportion_results$reanno))
all(unique(proportion_results$reanno) %in% final_anno_labels)

# Check for labels in ordered_labels that are not in proportion_results$reanno
mismatched_labels <- setdiff(final_anno_labels,unique(proportion_results$reanno))


# Get the unique clusters from your dataset
proportion_results$reanno <- factor(proportion_results$reanno, levels = final_anno_labels, ordered = TRUE)

p <- ggplot(data=proportion_results, aes(x=dataset, y=perct, fill=reanno)) +
  geom_bar(stat="identity", position="fill") +  # Use position="fill" for percentage stacking
  ggtitle("certain cell_type") +
  scale_fill_manual(values=humanref_color_mapping) + 
  theme(text=element_text(size=10)) 

p

# Save the plot as a PDF
ggsave("embryo_model_certain_celltype.pdf", plot = p, device = "pdf", width = 10, height = 4)





