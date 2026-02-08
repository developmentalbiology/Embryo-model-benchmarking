
setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3")
outputDir = getwd()

library(Seurat)
library(tidyverse)
library(spacexr)
library(Matrix)

# load data
cs8 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final data/public dataset/cs8_human_embryo.rds")
RCTD <- readRDS("RCTD_human_cs8_leiden_3_20250724_v3.rds")

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

#############################################################################
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


