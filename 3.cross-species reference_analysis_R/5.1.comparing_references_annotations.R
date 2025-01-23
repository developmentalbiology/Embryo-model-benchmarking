
########import packages and seurat objects of dataset###########
library(plotly)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(genefilter)
library(readxl)
library(devtools)
library(RColorBrewer)
library(tidyr)

# load human reference data
human <- readRDS("human_ref.rds")

# correct human$reanno
human$reanno[which(human$reanno=="Primitive.streak")] <- "Primitive streak"


# load monkey reference data
monkey <- readRDS("monkey_ref.rds")

# correct monkey$reanno
monkey$reanno[which(monkey$reanno=="Primitive.streak")] <- "Primitive streak"
monkey$reanno[which(monkey$reanno=="Node")] <- "Notochord"
monkey$reanno[which(monkey$reanno=="Paraxial.meso")] <- "Paraxial mesoderm"


# load crossspecies reference data
crossspecies <- readRDS("human_monkey_ref.rds")

# Define the ordered labels
crossspecies_labels = c('CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4',
                        'Epi_1', 'Epi_2', 'Epi_3',
                        'Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm',  'Exe.endothelium', 
                        'Amnion', 'Amniotic_epi',  'Ectoderm_1', 'Ectoderm_2',
                        'Neural tube', 'Neural crest',
                        'Primitive streak', 'Nascent mesoderm','PGC',
                        'Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
                        'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4',
                        'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm',
                        'Cardiac myocyte', 
                        'Notochord', 'DE', 'Gut',
                        'Hypoblast', 'AVE', 'YS.Endoderm_1', 'YS.Endoderm_2', 
                        'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor'
                   
)

# error check
all(crossspecies_labels %in% unique(crossspecies$crossspecies.anno))
all(unique(crossspecies$crossspecies.anno) %in% crossspecies_labels)

# summarize crossspecies annotations by species
meta <- crossspecies@meta.data

# define function for calculation of percentage
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

# calculate cell_type percentage by species
percentage <- calculate_percentage(meta, "crossspecies.anno", "species")

percentage$crossspecies.anno <- factor(percentage$crossspecies.anno, levels = crossspecies_labels)
#percentage$species<- factor(percentage$species, levels = c("Macaca fascicularis","Homo sapiens" ))

# Define your custom color palette
color_palette <- c("Homo sapiens"="#1F77B4", "Macaca fascicularis" = "#FF7F0E")


p <- ggplot(data=percentage, aes(x=crossspecies.anno, y=perct, fill=species)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),axis.text.x = element_text(angle = 90, hjust = 1))  
p
# Save the plot as a PDF
ggsave("crossspeces_ref_species_summary.pdf", plot = p, device = "pdf", width = 7, height = 4)


# calculate lineage percentage by species
percentage <- calculate_percentage(meta, "crossspecie.lineage", "species")

percentage$crossspecie.lineage <- factor(percentage$crossspecie.lineage)

p <- ggplot(data=percentage, aes(x=crossspecie.lineage, y=perct, fill=species)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),axis.text.x = element_text(angle = 90, hjust = 1))  
p
# Save the plot as a PDF
ggsave("crossspeces_ref_species_lineage_summary.pdf", plot = p, device = "pdf", width = 7, height = 4)


# plot
p <- DimPlot(crossspecies, reduction = "umap", group.by = "species",cols = color_palette)
p

# Save the plot as a PDF
ggsave("crossspeces_ref_species_umap.pdf", plot = p, device = "pdf", width = 10, height = 7)

##################check cell types that only have monkey cells#####################
celltypes <- percentage$crossspecies.anno[which(percentage$perct==1)]  
celltypes

##############################compare human annotations#############################################
meta <- crossspecies@meta.data
meta <- meta[which(meta$species=="Homo sapiens"),c("crossspecies.anno","crossspecie.lineage")]

meta_human <- human@meta.data
meta_human <- meta_human[,c("reanno","lineage")]

# error check
rownames(meta)[which(rownames(meta)!=rownames(meta_human))]
rownames(meta_human)[which(rownames(meta)!=rownames(meta_human))]

# fix cell names
rownames(meta_human)[which(rownames(meta)!=rownames(meta_human))] <- paste0(rownames(meta_human)[which(rownames(meta)!=rownames(meta_human))],"_1")
any(rownames(meta)!=rownames(meta_human))

meta <- merge(meta,meta_human, by='row.names',all=TRUE)

# compare annotation levels
human_reanno <- unique(meta$reanno)
crossspecies_human_reanno <- unique(meta$crossspecies.anno)

human_reanno_brief <- gsub("_[0-9]+$", "",human_reanno)
crossspecies_human_reanno_brief <- gsub("_[0-9]+$", "",crossspecies_human_reanno)

# error check
setdiff(human_reanno_brief, crossspecies_human_reanno_brief)
setdiff(crossspecies_human_reanno_brief,human_reanno_brief)

all(meta$lineage %in% meta$crossspecie.lineage)
all(meta$crossspecie.lineage %in% meta$lineage)

# add brief columns
meta$human_reanno_brief <- gsub("_[0-9]+$", "", meta$reanno)
meta$crossspecies_human_reanno_brief <- gsub("_[0-9]+$", "", meta$crossspecies.anno)


meta$consis_reanno <- ifelse(meta$human_reanno_brief == meta$crossspecies_human_reanno_brief, TRUE, FALSE)
meta$consis_lineage <- ifelse(meta$lineage == meta$crossspecie.lineage, TRUE, FALSE)


# Calculate the percentages for 'consis_reanno'
reanno_table <- table(meta$consis_reanno)
reanno_percentages <- prop.table(reanno_table) * 100

# Calculate the percentages for 'consis_lineage'
lineage_table <- table(meta$consis_lineage)
lineage_percentages <- prop.table(lineage_table) * 100

# Combine both percentage results into a dataframe
summary_df <- data.frame(
  Category = c("consis_reanno", "consis_lineage"),
  True_Percentage = c(reanno_percentages["TRUE"], lineage_percentages["TRUE"]),
  False_Percentage = c(reanno_percentages["FALSE"], lineage_percentages["FALSE"])
)

# plot
# Reorganize the dataframe into long format using pivot_longer
summary_df <- pivot_longer(summary_df, 
                        cols = c("True_Percentage" , "False_Percentage"), 
                        names_to = "variable", 
                        values_to = "percentage")

# Define your custom color palette
color_palette <- c("True_Percentage"="#FFB500", "False_Percentage" = "#c7c7c7")

p <- ggplot(data=summary_df, aes(x=Category, y=percentage, fill=variable)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette)+
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
# Save the plot as a PDF
ggsave("human_anno_compare.pdf", plot = p, device = "pdf", width = 3, height = 4)


#############################calculate human consistent and inconsistent annotation############################################

meta_check_lineage <- meta[which(meta$consis_lineage=="FALSE"),]
meta_ok_lineage <- meta[which(meta$consis_lineage=="TRUE"),]

meta_check_lineage_table <- as.data.frame(table(meta_check_lineage$lineage))
meta_ok_lineage_table <- as.data.frame(table(meta_ok_lineage$lineage))

# Rename the columns for clarity
colnames(meta_check_lineage_table) <- c("lineage", "count")
colnames(meta_ok_lineage_table) <- c("lineage", "count")

meta_check_lineage_table$group <- "inconsistent"
meta_ok_lineage_table$group <- "consistent"

# Merge the two dataframes by lineage
meta_check_all <- rbind(meta_check_lineage_table,meta_ok_lineage_table)

# Calculate the total counts for each lineage
meta_check_all$total_count <- ave(meta_check_all$count, meta_check_all$lineage, FUN = sum)

# Calculate the percentage by dividing the count by the total count
meta_check_all$percentage <- (meta_check_all$count / meta_check_all$total_count) * 100


# Define your custom color palette
color_palette <- c("TE_TrB"="#E377C2", 
                   "ExE_endo"="#FF7F0E",
                   "epi"="#B5BD61",
                   "mesoderm"="#AEC7E8", 
                   "neural_ecto"="#FFBB78",
                   "Gastru"="#D62728", 
                   "hemogenic"="#17BECF", 
                   "non_neuro_ecto"="#98DF8A",
                   "Endoderm"="#1F77B4",
                   "Notochord"="#AA40FC", 
                   "Exe_meso"= "#279E68",
                   "PGC"= "#8C564B"   )


# Initialize fill_color with default value
meta_check_all$fill_color <- "grey"  # Set default to grey for inconsistent groups

# Assign colors based on lineage for consistent groups
for (lineage in names(color_palette)) {
  meta_check_all$fill_color[meta_check_all$lineage == lineage & 
                              meta_check_all$group == "consistent"] <- color_palette[lineage]
}


meta_check_all$lineage <- factor(meta_check_all$lineage, levels = c("TE_TrB", "epi","Gastru", "Notochord", "PGC", "ExE_endo", "Exe_meso", "non_neuro_ecto", "neural_ecto", "Endoderm", "mesoderm","hemogenic" ))
  
# Create the plot
p <- ggplot(meta_check_all, aes(x=lineage, y=percentage, fill=fill_color)) +
  geom_bar(stat="identity", position="stack") +
  labs(title="human_reanno check", x="Lineage", y="Percentage") +
  scale_fill_identity() +  # Use the fill_color directly
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p

# Save the plot as a PDF
ggsave("human_reanno_check.pdf", plot = p, device = "pdf", width = 7, height = 4)


################further analysis on human exe_mesoderm cells##################################

exe_meso <- meta[which(meta$lineage=="Exe_meso" & meta$consis_lineage=="FALSE"),]

# Calculate percentage
df <- as.data.frame(table(exe_meso$crossspecies_human_reanno_brief))
df$percentage <- df$Freq / sum(df$Freq)
df$group <- "Exe_meso"


exe_meso_sub <- exe_meso[,c("reanno","crossspecies.anno" )]

# Plot the Sankey diagram

library(ggforce)
library(dplyr)
library(tidyr)

# Prepare the data
dataset <- exe_meso_sub
label_all <- dataset
ann_names <- colnames(label_all)
ann_size <- length(ann_names)

# Create count matrix
count_mat <- label_all %>% dplyr::count(across(all_of(ann_names)))

# Prepare data for plotting
test_gr <- gather_set_data(count_mat, 1:ann_size)
test_gr$x <- factor(test_gr$x, labels = ann_names)

# Check if test_gr has rows
if (nrow(test_gr) == 0) {
  stop("The test_gr data is empty. Check your data.")
}

my_color <- c("Exe.meso progenitor"="#1ABC9C", 
              "YS.mesoderm"="#1E8449",
              "Exe.endothelium"="#A3E4D7" )


# Plot
P <- ggplot(test_gr, aes(x = x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = reanno), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(color = "red", angle = 0, nudge_x = 0.07, hjust = 0, size = 4) +
  theme_classic(base_size = 20) + 
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 15), # hjust = 1, move the x-legend left
    plot.margin = margin(10, 10, 30, 10)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Increase space on the top
  scale_x_discrete(expand = expansion(mult = c(0, 0.2))) + # Increase space on the right
  labs(x = "", y = "") +
  scale_fill_manual(values = my_color) +
  annotate(geom = "segment", x = 4.25, xend = 4.25, y = 10, yend = 61, lwd = 2)

P

ggsave(filename = "exe_meso.pdf", plot = P, width = 10, height = 7, units = "in", dpi = 300)


##############################compare monkey annotations#############################################
meta <- crossspecies@meta.data
meta <- meta[which(meta$species=="Macaca fascicularis"),c("crossspecies.anno","crossspecie.lineage")]

meta_monkey <- monkey@meta.data
meta_monkey <- meta_monkey[,c("reanno","lineage")]

# error check
rownames(meta)[which(rownames(meta)!=rownames(meta_monkey))]
rownames(meta_monkey)[which(rownames(meta)!=rownames(meta_monkey))]

# fix cell names
rownames(meta_monkey)[which(rownames(meta)!=rownames(meta_monkey))] <- paste0(rownames(meta_monkey)[which(rownames(meta)!=rownames(meta_monkey))],"_1")
any(rownames(meta)!=rownames(meta_monkey))

meta <- merge(meta,meta_monkey, by='row.names',all=TRUE)

# compare annotation levels
monkey_reanno <- unique(meta$reanno)
crossspecies_monkey_reanno <- unique(meta$crossspecies.anno)

monkey_reanno_brief <- gsub("_[0-9]+$", "",monkey_reanno)
crossspecies_monkey_reanno_brief <- gsub("_[0-9]+$", "",crossspecies_monkey_reanno)

# error check
setdiff(monkey_reanno_brief, crossspecies_monkey_reanno_brief)
setdiff(crossspecies_monkey_reanno_brief,monkey_reanno_brief)

all(meta$lineage %in% meta$crossspecie.lineage)
all(meta$crossspecie.lineage %in% meta$lineage)

# add brief columns
meta$monkey_reanno_brief <- gsub("_[0-9]+$", "", meta$reanno)
meta$crossspecies_monkey_reanno_brief <- gsub("_[0-9]+$", "", meta$crossspecies.anno)


meta$consis_reanno <- ifelse(meta$monkey_reanno_brief == meta$crossspecies_monkey_reanno_brief, TRUE, FALSE)
meta$consis_lineage <- ifelse(meta$lineage == meta$crossspecie.lineage, TRUE, FALSE)

# Calculate the percentages for 'consis_reanno'
reanno_table <- table(meta$consis_reanno)
reanno_percentages <- prop.table(reanno_table) * 100

# Calculate the percentages for 'consis_lineage'
lineage_table <- table(meta$consis_lineage)
lineage_percentages <- prop.table(lineage_table) * 100

# Combine both percentage results into a dataframe
summary_df <- data.frame(
  Category = c("consis_reanno", "consis_lineage"),
  True_Percentage = c(reanno_percentages["TRUE"], lineage_percentages["TRUE"]),
  False_Percentage = c(reanno_percentages["FALSE"], lineage_percentages["FALSE"])
)

# plot
# Reorganize the dataframe into long format using pivot_longer
summary_df <- pivot_longer(summary_df, 
                           cols = c("True_Percentage" , "False_Percentage"), 
                           names_to = "variable", 
                           values_to = "percentage")

# Define your custom color palette
color_palette <- c("True_Percentage"="#FFB500", "False_Percentage" = "#c7c7c7")

p <- ggplot(data=summary_df, aes(x=Category, y=percentage, fill=variable)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette)+
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p
# Save the plot as a PDF
ggsave("monkey_anno_compare.pdf", plot = p, device = "pdf", width = 3, height = 4)

##############################calculate monkey consistent and inconsistent annotation############################################

meta_check_lineage <- meta[which(meta$consis_lineage=="FALSE"),]
meta_ok_lineage <- meta[which(meta$consis_lineage=="TRUE"),]

meta_check_lineage_table <- as.data.frame(table(meta_check_lineage$lineage))
meta_ok_lineage_table <- as.data.frame(table(meta_ok_lineage$lineage))

# Rename the columns for clarity
colnames(meta_check_lineage_table) <- c("lineage", "count")
colnames(meta_ok_lineage_table) <- c("lineage", "count")

meta_check_lineage_table$group <- "inconsistent"
meta_ok_lineage_table$group <- "consistent"

# Merge the two dataframes by lineage
meta_check_all <- rbind(meta_check_lineage_table,meta_ok_lineage_table)
  
# Calculate the total counts for each lineage
meta_check_all$total_count <- ave(meta_check_all$count, meta_check_all$lineage, FUN = sum)

# Calculate the percentage by dividing the count by the total count
meta_check_all$percentage <- (meta_check_all$count / meta_check_all$total_count) * 100


# Define your custom color palette
color_palette <- c("TE_TrB"="#E377C2", 
                   "ExE_endo"="#FF7F0E",
                   "epi"="#B5BD61",
                   "mesoderm"="#AEC7E8", 
                   "neural_ecto"="#FFBB78",
                   "Gastru"="#D62728", 
                   "hemogenic"="#17BECF", 
                   "non_neuro_ecto"="#98DF8A",
                   "Endoderm"="#1F77B4",
                   "Notochord"="#AA40FC", 
                   "Exe_meso"= "#279E68",
                   "PGC"= "#8C564B"   )


# Initialize fill_color with default value
meta_check_all$fill_color <- "grey"  # Set default to grey for inconsistent groups

# Assign colors based on lineage for consistent groups
for (lineage in names(color_palette)) {
  meta_check_all$fill_color[meta_check_all$lineage == lineage & 
                              meta_check_all$group == "consistent"] <- color_palette[lineage]
}

meta_check_all$lineage <- factor(meta_check_all$lineage, levels = c("TE_TrB", "epi","Gastru", "Notochord", "PGC", "ExE_endo", "Exe_meso", "non_neuro_ecto", "neural_ecto", "Endoderm", "mesoderm","hemogenic" ))

# Create the plot
p <- ggplot(meta_check_all, aes(x=lineage, y=percentage, fill=fill_color)) +
  geom_bar(stat="identity", position="stack") +
  labs(title="human_reanno check", x="Lineage", y="Percentage") +
  scale_fill_identity() +  # Use the fill_color directly
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p

#Save the plot as a PDF
ggsave("monkey_reanno_check.pdf", plot = p, device = "pdf", width = 7, height = 4)


###########################further compare human and monkey ExE.meso#########################################################

human_meta <- human@meta.data[which(human@meta.data$lineage=="Exe_meso"),]
monkey_meta <- monkey@meta.data[which(monkey@meta.data$lineage=="Exe_meso"| monkey@meta.data$reanno %in% c("Lateral plate mesoderm_1",
                                                                                                           "Lateral plate mesoderm_2",
                                                                                                           "Lateral plate mesoderm_3",
                                                                                                           "Lateral plate mesoderm_4")    ),]

#define function for percentage calculation
calculate_percentage <- function(data, stage, reanno) {
  data_sum <- data %>%
    group_by(stage, reanno) %>%
    summarize(frequency = n(), .groups = 'drop') %>%
    arrange(stage, reanno)
  percentage <- do.call(rbind, lapply(unique(data_sum$reanno), function(x) {
    
    df <- data_sum[which(data_sum$reanno==x),]
    df$perct <- df$frequency/sum(df$frequency)
    
    return(df)
  }))
  
  return(percentage)
}


human_meta_sum <- calculate_percentage(human_meta, "stage", "reanno")
monkey_meta_sum <- calculate_percentage(monkey_meta, "stage", "reanno")


ggplot(data=human_meta_sum, aes(x=reanno, y=perct, fill=stage)) +
  geom_bar(stat="identity") +
  ggtitle("reanno")


ggplot(data=monkey_meta_sum, aes(x=reanno, y=perct, fill=stage)) +
  geom_bar(stat="identity") +
  ggtitle("reanno")

######################
library(ggpubr)

zeileis_25 <- c(
  "#023fa5",
  "#7d87b9",
  "#bec1d4",
  "#d6bcc0",
  "#bb7784",
  "#8e063b",
  "#4a6fe3",
  "#8595e1",
  "#b5bbe3",
  "#e6afb9",
  "#e07b91",
  "#d33f6a",
  "#11c638",
  "#8dd593",
  "#c6dec7",
  "#ead3c6",
  "#f0b98d",
  "#ef9708",
  "#0fcfc0",
  "#9cded6",
  "#d5eae7",
  "#f3e1eb",
  "#f6c4e1",
  "#f79cd4",
  "#7f7f7f"
)

# Define the ordered levels for the 'stage' column
ordered_stages <- c("E11_IVC", "E12_IVC", "E13_IVC", "E14_IVC", 
                    "E15_IVC", "E16_IVC", "E17_IVC", "E18_IVC", "E19_IVC", "E20_IVC", "E21_IVC", 
                    "E22_IVC", "E23_IVC", "E24_IVC", "E25_IVC", 
                    "E13", "E14","CS7", "CS8", "CS9", "CS11")

# Convert 'stage' column to a factor with the specified order
human_meta_sum$stage <- factor(human_meta_sum$stage, levels = ordered_stages)
monkey_meta_sum$stage <- factor(monkey_meta_sum$stage, levels = ordered_stages)

# Define a color palette that has at least as many colors as unique stages
if (length(ordered_stages) > length(zeileis_25)) {
  stop("Not enough colors for each unique stage.")
}

# Create a named color vector based on the ordered stages
stage_colors <- setNames(zeileis_25[1:length(ordered_stages)], ordered_stages)


data=human_meta_sum

unique_reannos <- unique(data$reanno)
for (anno in unique_reannos) {
  # Filter data for the current reanno
  data_subset <- subset(data, reanno == anno)
  
  p <- ggdonutchart(data_subset, "frequency", label = "reanno",
                    fill = "stage", color = "white",
                    palette = stage_colors )+
    ggtitle(paste("Donut Chart for", reanno)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Print or save the plot
  print(p)  # Display in the R console
  # save each plot as a file
  ggsave(filename = paste0("donut_chart_human_", anno, ".pdf"), plot = p)
  
  
}

