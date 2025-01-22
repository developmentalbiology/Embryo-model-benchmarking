library(Seurat)
library(tidyverse)
library(spacexr)
library(Matrix)


# RCTD results
#reload dataset
cs8 <- readRDS("cs8_human_embryo.rds")
RCTD <- readRDS("RCTD_cs8_human_final_anno_fullmode_celltype_20250108.rds")


#get xy coordinates
xy <- cs8@meta.data[,c("x","y")]
xy <- as.matrix(xy)

#RCTD results
results <- RCTD@results
norm_weights <- normalize_weights(results$weights)
cell_type_names <- RCTD@cell_type_info$info[[2]]
spatialRNA <- RCTD@spatialRNA
xy <- xy[row.names(spatialRNA@coords), ]
spatialRNA@coords <- as.data.frame(xy)


##furter analysis
df <- as.data.frame(norm_weights)
df <- df[match(rownames(df), rownames(cs8@meta.data)),]

# Use apply to find the maximum value per row
df2 <- data.frame(
  max_value = apply(df, 1, max),                       # Maximum value per row
  max_column = colnames(df)[apply(df, 1, which.max)]   # Column name with the max value
)


#plot percentage
perct <- as.data.frame(table(df2$max_column))

getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
colourCount = length(unique(df2$max_column))

ggplot(data=perct, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") +
  ggtitle("human_CS8_RCTD-reanno") +
  scale_fill_manual(values = getPalette(colourCount))+
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

cs8$celltype <- df2$max_column
cs8$max_score <- df2$max_value

# Plot UMAP with 'max_score' 
# Convert xy matrix to a dimension reduction object
colnames(xy) <- c("spatial_1", "spatial_2")
# Create the spatial reduction object
spatial_reduction <- CreateDimReducObject(embeddings = xy, key = "spatial_", assay = DefaultAssay(cs8))

# Add the 'spatial' reduction to the Seurat object
cs8@reductions[["spatial"]] <- spatial_reduction


#subset 3x3

#Extract UMAP coordinates from the Seurat object
umap_coords <- Embeddings(cs8, reduction = "spatial")
range(umap_coords[, 1])
range(umap_coords[, 2])

# Subset cells based on x and y coordinates
cells_to_keep <- rownames(umap_coords[umap_coords[, 1] >= 9000 & umap_coords[, 1] <= 17000 &
                                        umap_coords[, 2] >= 8000 & umap_coords[, 2] <= 19800, ])

# Subset cells based on x and y coordinates
#cells_to_keep <- rownames(umap_coords[umap_coords[, 1] >= 9000 & umap_coords[, 1] <= 17000 &
#                                        umap_coords[, 2] >= 28000 & umap_coords[, 2] <= 40000, ])


# Subset the Seurat object
cs8_subset <- subset(cs8, cells = cells_to_keep)


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
                      'Hypoblast', 'AVE', 'VE_YE', 'YS.Endoderm_1', 'YS.Endoderm_2', 
                      'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor'
                      
)


# Reorder the 'first_type' column according to the custom list and remove those that don't exist
cs8_subset$celltype <- factor(cs8_subset$celltype, levels = final_anno_labels, ordered = TRUE)

# Optionally, you can drop unused levels
cs8_subset$celltype <- droplevels(cs8_subset$celltype)


# Get the unique clusters from your dataset
unique_clusters <- unique(cs8_subset$celltype)

# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)


# Generate the plot with the custom palette
p <- DimPlot(cs8_subset, reduction = "spatial", group.by = "celltype") +
  scale_color_manual(values = color_palette) +  # Apply the custom color palette
  ggtitle("UMAP Plot with Custom Color Palette") 

pdf(paste0(outputDir, "/human_CS8_RCTD_full_top.pdf"), width = 12, height = 7)
p

dev.off()

#####################export original annotation################################
# Define cell types
clusters <- c("AM","AM.EXE.Meso", "Endo","Epi/Ecto", 
                "Ery", "Gast/PS", "HEP",  "Meso", "Noto",
                "Visceral.Endo", "YS.Endo", "YS.EXE.Meso-A", 
                "YS.EXE.Meso-B" 
                )

# Define corresponding colors
colors <- c("#ffff00", "#1ce6ff", "#ff34ff", "#ff4a46", "#008941", 
            "#006fa6", "#a30059", "#ffdbe5", "#7a4900", "#0000a6", 
            "#63ffac", "#b79762", "#004d43")

# Assign colors to cell types as a named vector
cell_type_colors <- setNames(colors, clusters)

p <- DimPlot(cs8_subset, reduction = "spatial", group.by = "clusters") +
  scale_color_manual(values = cell_type_colors) +  # Apply the custom color palette
  ggtitle("UMAP Plot with Custom Color Palette") 

pdf(paste0(outputDir, "/human_CS8_RCTD_full_top_original.pdf"), width = 12, height = 7)
p

dev.off()


#######################plot norm.weight#########################################
# Make plots
plot_puck_wrapper <- function(puck, plot_val, cell_type = NULL, minUMI = 0, maxUMI = 200000, min_val = NULL, max_val = NULL, title = NULL, my_cond = NULL) {
  UMI_filter = (puck@nUMI >= minUMI) & (puck@nUMI < maxUMI)
  ylimit = NULL
  if(!is.null(my_cond))
    my_cond = UMI_filter & my_cond
  else
    my_cond = UMI_filter
  if(!is.null(cell_type))
    my_cond = my_cond & (puck@cell_labels == cell_type)
  if(!is.null(min_val))
    my_cond = my_cond & (plot_val > min_val)
  if(!is.null(max_val)) {
    epsilon = 0.00001
    plot_val[plot_val >= max_val - epsilon] = max_val - epsilon
    if(!is.null(min_val))
      ylimit = c(min_val, max_val)
  }
  
  p <- plot_puck_continuous(puck, names(which(my_cond)), plot_val, title = title, ylimit = ylimit)
  
  # Add the color scale
  p <- p +
    ggplot2::scale_color_gradientn(colors = c("blue", "yellow", "red"), limits = c(0, 1))
  
  # Add theme to change the background color to black
  p <- p +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "black"),
      plot.background = ggplot2::element_rect(fill = "black"),
      plot.title = ggplot2::element_text(color = "white", size = 36)
      # text = ggplot2::element_text(color = "white"),
      # axis.title.x=element_blank(),
      # axis.title.y=element_blank(),
      # axis.ticks=element_blank(),
      # axis.text.x=element_blank(),
      # axis.text.y=element_blank(),
      # legend.position="none"
    )
  
  return(p)
}

plots <- vector(mode = "list", length = length(cell_type_names))
weights <- norm_weights
puck <- spatialRNA

for (i in 1:length(cell_type_names)) {
  cell_type = cell_type_names[i]
  plot_var <- weights[,cell_type]; names(plot_var) = rownames(weights)
  if(sum(weights[,cell_type]) > 0)
    plots[[i]] <- plot_puck_wrapper(puck, plot_var, NULL, minUMI = 0, maxUMI = 200000, min_val = 0, max_val = 1, title = cell_type)
}

pdf("human_CS8_RCTD_score.pdf", width = 25, height = 10)
invisible(lapply(plots, print))
dev.off()


# Loop through each cell type and save as a separate PDF
for (i in 1:length(cell_type_names)) {
  cell_type <- cell_type_names[i]
  plot_var <- weights[, cell_type]
  names(plot_var) <- rownames(weights)
  
  if (sum(weights[, cell_type]) > 0) {
    plot <- plot_puck_wrapper(
      puck,
      plot_var,
      NULL,
      minUMI = 0,
      maxUMI = 200000,
      min_val = 0,
      max_val = 1,
      title = cell_type
    )
    
    # Save the plot as a PDF
    pdf(file = paste0("human_CS8_RCTD_score_", cell_type, ".pdf"), width = 25, height = 10)
    print(plot)
    dev.off()
  }
}


#############################################################################
# Subset the Seurat object
cs8 <- NormalizeData(cs8)
cs8_subset <- subset(cs8, cells = cells_to_keep)

cs8$celltype_expression <- "NA"

######################################calculate distribution percenatge############################################

# Label cells as "HEP" 
threshold <- 1

cs8_expression <- FetchData(cs8, vars = c("PECAM1", "GATA1", "GP1BB","ITGA2B","NFE2"), layer = "data")
# Calculate row means for the selected genes and store in a new column
cs8_expression[["Expression_Mean"]] <- rowMeans(cs8_expression)

cell.labels <- rownames(cs8_expression)[which(cs8_expression$Expression_Mean>threshold)]

cs8$celltype_expression[which(rownames(cs8@meta.data) %in% cell.labels)] <- "H.E.P"



# Label cells as "Erythroid" 
threshold <- 3.5

cs8_expression <- FetchData(cs8, vars = c("KLF1","GYPB","HBE1","HBZ"), layer = "data")
# Calculate row means for the selected genes and store in a new column
cs8_expression[["Expression_Mean"]] <- rowMeans(cs8_expression)

cell.labels <- rownames(cs8_expression)[which(cs8_expression$Expression_Mean>threshold)]

cs8$celltype_expression[which(rownames(cs8@meta.data) %in% cell.labels)] <- "Erythroid"



# plot expression_celltype
cs8_subset <- subset(cs8, cells = cells_to_keep)

colors <- c("H.E.P" = "blue", "Erythroid" = "red", "NA" = "grey")

p <- DimPlot(cs8_subset, reduction = "spatial", group.by = "celltype_expression",cols = colors )

pdf(paste0(outputDir, "/expression_celltype.pdf"), width = 12, height = 7)
p

dev.off()


##########################check true and false percentage###################################

meta <- cs8@meta.data
meta_sub <- meta[,c("clusters","celltype","celltype_expression")]

meta_sub2 <- meta_sub[which(meta_sub$celltype_expression=="H.E.P"),]

meta_sub2$clusters[which(meta_sub2$clusters=="HEP")] <- "H.E.P"
meta_sub2$celltype[which(meta_sub2$celltype=="Hemogenic endothelial progenitors")] <- "H.E.P"
meta_sub2$celltype[which(meta_sub2$celltype=="Hemogenic endothelial progenitors")] <- "H.E.P"

# Calculate the percentages for each column compared to 'celltype_expression'
percentage_col1 <- length(which(meta_sub2[[1]] == meta_sub2$celltype_expression)) / nrow(meta_sub2) * 100
percentage_col2 <- length(which(meta_sub2[[2]] == meta_sub2$celltype_expression)) / nrow(meta_sub2) * 100

# Create a dataframe with the calculated percentages
percentage_df_HEP <- data.frame(
  Column = c(names(meta_sub2)[1], names(meta_sub2)[2] ),
  Percentage = c(percentage_col1, percentage_col2)
)
percentage_df_HEP$group <- "HEP"

#############################################################
meta_sub2 <- meta_sub[which(meta_sub$celltype_expression=="Erythroid"),]

meta_sub2$clusters[which(meta_sub2$clusters=="Ery")] <- "Erythroid"

# Calculate the percentages for each column compared to 'celltype_expression'
percentage_col1 <- length(which(meta_sub2[[1]] == meta_sub2$celltype_expression)) / nrow(meta_sub2) * 100
percentage_col2 <- length(which(meta_sub2[[2]] == meta_sub2$celltype_expression)) / nrow(meta_sub2) * 100

# Create a dataframe with the calculated percentages
percentage_df_ERY <- data.frame(
  Column = c(names(meta_sub2)[1], names(meta_sub2)[2]),
  Percentage = c(percentage_col1, percentage_col2)
)
percentage_df_ERY$group <- "Erythroid"

percentage_df <- rbind(percentage_df_HEP, percentage_df_ERY)


# Define the color palette you want to use
color_palette <- c("clusters" = "#DDDDDD", "celltype" = "#8595E1" )

percentage_df$Column <- factor(percentage_df$Column, levels = c("clusters","celltype"))

# Create the boxplot
p<- ggplot(percentage_df, aes(x = group, y = Percentage, fill = Column)) +
  geom_col(position = position_dodge()) +
  #facet_wrap(~ label, strip.position = "top") +  # Use 'label' for facetting
  scale_fill_manual(values = color_palette) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "percentage of coverage")

p

# Save the plot as a PDF
ggsave("spatial_percentage of coverage.pdf", plot = p, device = "pdf", width = 7, height = 4)


######################subset another slice############################################
# Subset cells based on x and y coordinates
cells_to_keep <- rownames(umap_coords[umap_coords[, 1] >= 11500 & umap_coords[, 1] <= 14000 &
                                        umap_coords[, 2] >= 8000 & umap_coords[, 2] <= 19800, ])

# Subset the Seurat object
cs8_subset <- subset(cs8, cells = cells_to_keep)


sort(table(cs8$celltype), decreasing = TRUE)


##################gene expression###########################################

p <-FeaturePlot(cs8, features = c("HOXA10"), reduction="spatial")

pdf(paste0(outputDir, "/human_CS8_RCTD_HOXA10.pdf"), width = 12, height = 10)
p

dev.off()

