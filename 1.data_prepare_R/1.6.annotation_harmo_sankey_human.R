setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/1. dataset_prepare")
outputDir = getwd()

human <- readRDS("human_all_2025.rds")

# fix colnames
colnames(human$yuan@meta.data)[11] <- "sub_anno"
colnames(human$yuan@meta.data)[12] <- "harmo.anno"

# load library and data
library(Seurat)
library(ggforce)
library(dplyr)
library(tidyr)

# generate label metadata
extract_metadata_labels <- function(raw_list) {
  
  label_list <- lapply(names(raw_list), function(dataset_name) {
    data <- raw_list[[dataset_name]]
    
    df <- tibble::tibble(
      dataset     = dataset_name,
      orig.anno   = data$anno,
      orig.sub_anno = data$sub_anno,
      harmo.anno  = data$harmo.anno
    )
    
    unique(df)
  })
  
  # Combine all into one dataframe
  label_all <- dplyr::bind_rows(label_list)
  
  return(label_all)
}

label_all <- extract_metadata_labels(raw_list = human)

label_all <- label_all[,-1]


label_all$orig.anno[which(label_all$orig.anno=="NA")]<- "Unknown"
label_all$orig.sub_anno[which(label_all$orig.sub_anno=="NA")]<- "Unknown"
label_all$harmo.anno[which(label_all$harmo.anno=="NA")]<- "Unknown"


# Prepare the data
label_all$orig.anno <- factor(label_all$orig.anno, levels = c(
  "Epiblast", "Hypoblast", "TE_CTB_STB", "NMP", "NSC",
  "MesoPro", "Myoblast", "LPM", "Neural Crest", "Endothelium",
  "Epithelium", "Erythroid", "Immune Cell", "Endoderm",
  "Notochord", "IPC", "DA N", "Glu N", "Chondroblast",
  "TE_AME", "AVE_YE", "Primitive.streak", "ExE.Mesoderm",
  "Hemogenic Endothelial Progenitors", "ExE Mesoderm",
  "Erythroblasts", "Primitive Streak", "Advanced Mesoderm",
  "Non-Neural Ectoderm", "Emergent Mesoderm", "Nascent Mesoderm",
  "Axial Mesoderm", "ICM", "EPI", "CTB", "STB", "EVT",
  "PSA-EPI", "TE", "PE", "MIX", "Unknown"
))

label_all$orig.sub_anno <- factor(label_all$orig.sub_anno, levels = c(
  "Epiblast", "Hypoblast", "TE_CTB_STB", "Lateral plate Mesoderm 1",
  "Neural tbue 1", "Intermediate Mesoderm", "Somite",
  "Lateral plate Mesoderm 2", "NMP", "Paraxial mesoderm",
  "Neural Crest", "Endothelium", "Epithelium", "Presomitic Mesoderm",
  "Myocyte Progenitor", "Neural tbue 3", "Mesenchyme",
  "Neural tbue 2", "Endoderm", "Erythroid", "Cardiac Myocyte",
  "Notochord", "TE_AME", "AVE_YE", "Primitive.streak",
  "ExE.Mesoderm", "Blood Progenitors", "YS Mesoderm",
  "Erythro-Myeloid Progenitors", "Hemogenic Endothelium",
  "Erythroblasts", "Primitive Streak", "Advanced Mesoderm",
  "Myeloid Progenitors", "YS Endoderm", "DE(P)", "NNE",
  "Emergent Mesoderm", "DE(NP)", "Amnion", "PGC",
  "Nascent Mesoderm", "Axial Mesoderm", "ICM", "EPI",
  "CTB", "STB", "EVT", "PSA-EPI", "TE", "PE", "MIX", "Unknown"
))

label_all$harmo.anno <- factor(label_all$harmo.anno, levels = c(
  "Epiblast", "Hypoblast", "TE/TrB", "Mesoderm", "Neural.Ectoderm",
  "Primitive.streak", "Hemogenic.Endothelial.Progenitors",
  "Non-Neural.Ectoderm", "Endoderm", "Notochord", "ExE.Endoderm",
  "ExE.Mesoderm", "PGC", "ICM", "Unknown"
))

ann_names <- colnames(label_all)
ann_size <- length(ann_names)

# Create count matrix
count_mat <- label_all %>% dplyr::count(across(all_of(ann_names)))

# Check if count matrix has rows
if (nrow(count_mat) == 0) {
  stop("The count matrix is empty. Check your data.")
}


colors_15 <- c(
  "#A349A4", "#FFFF33", "#E7298A", "#1B9E77", "#FDB462",
  "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#8DD3C7",
  "#9F000F", "#BEBADA", "#FB8072", "#80B1D3",  "#091833"
 
)

# Assign names
names(colors_15) <- unique(label_all$harmo.anno)

# Map colors to unique values in each column
color_list <- list()
count <- 1
for (ann in ann_names) {
  unique_values <- unique(label_all[[ann]])
  count_max <- count + length(unique_values) - 1
  color_list[[ann]] <- setNames(color[count:count_max], unique_values)
  count <- count_max + 1
}

# Prepare data for plotting
test_gr <- gather_set_data(count_mat, 1:ann_size)
test_gr$x <- factor(test_gr$x, labels = ann_names)

# Check if test_gr has rows
if (nrow(test_gr) == 0) {
  stop("The test_gr data is empty. Check your data.")
}

# Plot the Sankey diagram
P <- ggplot(test_gr, aes(x = x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = factor(harmo.anno)), alpha = 0.8, axis.width = 0.15) +
  geom_parallel_sets_axes(axis.width = 0.15) +
  geom_parallel_sets_labels(color = "black", angle = 0, nudge_x = 0.07, hjust = 0, size = 3) +
  scale_fill_manual(values = colors_15, name = "harmo.anno") +
  labs(x = "", y = "") +
  theme_sankey

ggsave(filename = "human_label_sankey.pdf", plot = P, width = 10, height = 7, units = "in", dpi = 300)
