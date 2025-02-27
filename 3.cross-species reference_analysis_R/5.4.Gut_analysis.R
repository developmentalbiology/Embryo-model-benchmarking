library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readr)
library(tidyr)
library(dplyr)
library(reshape2)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250206_human_reference_annotation_comparasion")
outputDir = getwd()

#load human reference data
human <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final code/20250107_reference_comparasion/human_20250108.Rds")

# Gut
meta <- human@meta.data

Gut <- meta[which(meta$final_anno=="Gut" ),]

Gut <- Gut[,c("reanno","final_anno" )]

# Plot the Sankey diagram

library(ggforce)
library(dplyr)
library(tidyr)

# Prepare the data
dataset <- Gut
label_all <- dataset
ann_names <- colnames(label_all)
ann_size <- length(ann_names)

# Create count matrix
count_mat <- label_all %>% dplyr::count(across(all_of(ann_names)))

count_mat <- count_mat %>%
  mutate(
    reanno = as.character(reanno),  # Convert reanno to a character
    final_anno = as.character(final_anno)  # Convert final_anno to a character
  )


# Prepare data for plotting
test_gr <- gather_set_data(count_mat, 1:ann_size)
test_gr$x <- factor(test_gr$x, labels = ann_names)

# Check if test_gr has rows
if (nrow(test_gr) == 0) {
  stop("The test_gr data is empty. Check your data.")
}

my_color <- c("Amnion"="#6a3a4c", 
              "DE"="#00846f",
              "Lateral plate mesoderm_1" ="#b903aa",
              "Gut"="#372101"
              
)


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
  scale_fill_manual(values = my_color) 

P

ggsave(filename = "Sanky_Gut.pdf", 
       plot = P, device = "pdf", width = 5, height = 8)


# marker expression

human_sub <- subset(human, subset = final_anno %in% c("DE", "Gut") )

##calculate humanref_celltype mean
human_sub$final_anno <- as.factor(human_sub$final_anno)
Idents(human_sub) <- human_sub$final_anno
# Calculate average expression using AggregateExpression
avgexp <- AggregateExpression(human_sub, assays = "RNA", slot = "data")
avgexp <- avgexp %>% as.data.frame()


features = c("SOX17","GATA4","PDGFRA","FOXA2", 
             "FOXA1", "HOXA10", "CXCL12", "CDX2")  


# calculate de genes
#DE_Gut.markers <- FindMarkers(human_sub, ident.1 = "Gut", ident.2 = "DE")

#DE_Gut.markers %>%
#  dplyr::filter(p_val_adj < 0.0001) %>%
#  dplyr::filter(abs(avg_log2FC) > 2) -> DE_Gut.markers


#features =rownames(DE_Gut.markers)


plot_gene <- avgexp[features,]

library(pheatmap)
library(RColorBrewer)

RdBu_palette <- colorRampPalette(c(rev(brewer.pal(5, "Blues")), 'white', brewer.pal(5, "Reds")))
#RdBu_palette <- c(rev(brewer.pal(5, "Blues")), 'white', brewer.pal(5, "Reds"))

num_colors <- 100
colors <- RdBu_palette(num_colors)


p <- pheatmap(
  plot_gene,
  color = colors,
  cluster_rows = TRUE,  # Disable clustering for rows (genes)
  cluster_cols = TRUE,  # Disable clustering for columns (timepoints)
  main = "Log Normalized Gene Expression Heatmap",  # Title of the heatmap
  scale = "column",  
  fontsize_row = 8,  # Font size for gene names (rows)
  fontsize_col = 8   # Font size for timepoints (columns)
)

p

ggsave(filename = "heatmap_DE_Gut.pdf", 
       plot = p, device = "pdf", width = 5, height = 8)


