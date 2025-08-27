
# ============================================================================
# huamn CS8 final RCTD_based_annotation_assignment
# ============================================================================

setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3")
outputDir = getwd()

# RCTD results
#reload dataset
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

# Convert xy matrix to a dimension reduction object
colnames(xy) <- c("spatial_1", "spatial_2")
# Create the spatial reduction object
spatial_reduction <- CreateDimReducObject(embeddings = xy, key = "spatial_", assay = DefaultAssay(cs8))
# Add the 'spatial' reduction to the Seurat object
cs8@reductions[["spatial"]] <- spatial_reduction

# convert weights matrix
norm_weights <- as.matrix(norm_weights)

source("multi_tier_integrated_final.R") # load custom script for cell type assignment


results <- multi_tier_assignment_integrated_final(
  weights_matrix = norm_weights,
  cell_type_names = cell_type_names,
  spatial_coords = spatialRNA@coords,
  
  # Original 5-tier parameters (fully preserved)
  excluded_cell_types = c('TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2', 'Hypoblast_1' ),
  protected_cell_types = NULL,
  max_proportion = 0.1,                       # proportion limit
  k_neighbors = 7,                            # number of neighbors
  spatial_weight = 0.8,                       # spatial weight
  protection_threshold = 0.4,                 # protection threshold
  dominance_threshold = 0.2,                  # dominance threshold
  
  # Signal balancing parameters, supress broadly observed signal
  suppressed_cell_types = c('Epiblast_3', 'Neural.ectoder.forebrain', 'Neural.ectoderm.midbrain'),    
  max_suppressed_proportion = 0.1,           
  suppression_strength = 0.8,                 
  
  # Weak signal protection parameters, protect weak but specific signal
  weak_signal_protected_types = c('Ectoderm', 'Paraxial.mesoderm','Lateral.plate.mesoderm_2', 'Primitive.streak', 'DE','YS.mesoderm_2', 'YS.endoderm','Primitive.megakaryocyte'),
  weak_signal_threshold = 0.02,               # Lower protection threshold
  weak_signal_dominance_threshold = 1.01,     # More lenient dominance requirement
  weak_signal_spatial_boost = 3.0,            # 30% spatial weight boost
  
  # activate parameters
  enable_signal_balancing = TRUE,             # Enable signal balancing
  enable_optimizations = TRUE,                # Maintain high performance
  verbose = TRUE
)


# Verify TE is completely gone
table(results$final_assignment)  # Should show 0 TE

############## check plots ###################################
# add results to the query data
cs8$final_anno <- results$final_assignment
cs8$final_anno_confidence <- results$assignment_confidence

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
final_anno_labels = c('TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
                      'Epiblast_1','Epiblast_2','Epiblast_3','Ectoderm',
                      'Amniontic.epi','Amniontic.ectoderm',
                      'PGC',
                      'Primitive.streak',
                      'Neuromesodermal.progenitor',
                      'Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 'Neural.ectoderm.midbrain','Spinal.cord',
                      'Paraxial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 'Rostral.mesoderm', 'Lateral.plate.mesoderm_1',
                      'Lateral.plate.mesoderm_2','Lateral.plate.mesoderm_3','Cardiac.mesoderm','Amniotic.mesoderm','Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2',
                      'Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm',
                      'DE','Gut',
                      'Notochord',
                      'Hemogenic.endothelial.progenitor','Endothelium','Erythroid','Primitive.megakaryocyte','Myeloid.progenitor'
)



# Reorder the 'first_type' column according to the custom list and remove those that don't exist
cs8$final_anno <- factor(cs8$final_anno, levels = final_anno_labels, ordered = TRUE)

# Optionally, you can drop unused levels
cs8$final_anno <- droplevels(cs8$final_anno)

# Get the unique clusters from your dataset
unique_clusters <- unique(cs8$final_anno)

# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)


# Generate the plot with the custom palette
p <- DimPlot(cs8, reduction = "spatial", group.by = "final_anno") +
  scale_color_manual(values = color_palette) +  # Apply the custom color palette
  ggtitle("UMAP Plot with Custom Color Palette") 

pdf(paste0(outputDir, "/human_CS8_RCTD_full_top_optimized.pdf"), width = 18, height = 13)
p

dev.off()

saveRDS(
  object = cs8,
  file = "human_cs8_reanno.Rds"
)



# ============================================================================
# huamn CS7 final RCTD_based_annotation_assignment
# ============================================================================

setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3")
outputDir = getwd()

# RCTD results
#reload dataset
cs7 <- readRDS("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/final data/public dataset/CS7_human_embryo.rds")
RCTD <- readRDS( "RCTD_human_cs7_leiden_3_20250724_v3.rds")

# change x and y coordination
tile_width <- 1500
tile_height <- 3000
tiles_per_row <- 10

sample_num <- as.numeric(sub("S", "", cs7$sample_final))

col_index <- (sample_num - 1) %% tiles_per_row
row_index <- (sample_num - 1) %/% tiles_per_row

cs7$x_new <- cs7$newx + col_index * tile_width
cs7$y_new <- cs7$newy + row_index * tile_height


colnames(cs7@meta.data)[8] <- "x"
colnames(cs7@meta.data)[9] <- "y"


#get xy coordinates
xy <- cs7@meta.data[,c("x","y")]
xy <- as.matrix(xy)

#RCTD results
results <- RCTD@results
norm_weights <- normalize_weights(results$weights)
cell_type_names <- RCTD@cell_type_info$info[[2]]
spatialRNA <- RCTD@spatialRNA
xy <- xy[row.names(spatialRNA@coords), ]
spatialRNA@coords <- as.data.frame(xy)

# Convert xy matrix to a dimension reduction object
colnames(xy) <- c("spatial_1", "spatial_2")
# Create the spatial reduction object
spatial_reduction <- CreateDimReducObject(embeddings = xy, key = "spatial_", assay = DefaultAssay(cs8))
# Add the 'spatial' reduction to the Seurat object
cs7@reductions[["spatial"]] <- spatial_reduction

# convert weights matrix
norm_weights <- as.matrix(norm_weights)

# cell type assignment
results <- multi_tier_assignment_integrated_final(
  weights_matrix = norm_weights,
  cell_type_names = cell_type_names,
  spatial_coords = spatialRNA@coords,
  
  # Original 5-tier parameters (fully preserved)
  excluded_cell_types = c('TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2', 'Hypoblast_1' ),
  protected_cell_types = NULL,
  max_proportion = 0.1,                       
  k_neighbors = 6,                            
  spatial_weight = 0.8,                       
  protection_threshold = 0.4,                 
  dominance_threshold = 0.2,                  
  
  # Signal balancing parameters, supress broadly observed signal
  suppressed_cell_types = c('Emergent.mesoderm', 'Neural.ectoderm.forebrain',  'Neural.ectoderm.midbrain'),        
  max_suppressed_proportion = 0.05,           
  suppression_strength = 0.6,                 
  
  # NEW: Weak signal protection parameters
  weak_signal_protected_types = c('Ectoderm', 'Paraxial.mesoderm','Lateral.plate.mesoderm_2', 'Primitive.streak', 'DE','YS.mesoderm_2', 'YS.endoderm','Primitive.megakaryocyte'),
  weak_signal_threshold = 0.02,               
  weak_signal_dominance_threshold = 1.01,     
  weak_signal_spatial_boost = 3.0,            
  
  # Control parameters
  enable_signal_balancing = TRUE,             
  enable_optimizations = TRUE,                
  verbose = TRUE
)


table(results$final_assignment)  # Check assignments
table(results$assignment_tier)   # Check quality tiers


############## check plots ###################################
# add results to the query data
cs7$final_anno <- results$final_assignment
cs7$final_anno_confidence <- results$assignment_confidence

#assign colors

# Reorder the 'first_type' column according to the custom list and remove those that don't exist
cs7$final_anno <- factor(cs7$final_anno, levels = final_anno_labels, ordered = TRUE)

# Optionally, you can drop unused levels
cs7$final_anno <- droplevels(cs7$final_anno)

# Get the unique clusters from your dataset
unique_clusters <- unique(cs7$final_anno)

# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)


# Generate the plot with the custom palette
p <- DimPlot(cs7, reduction = "spatial", group.by = "final_anno") +
  scale_color_manual(values = color_palette) +  # Apply the custom color palette
  ggtitle("UMAP Plot with Custom Color Palette") 

pdf(paste0(outputDir, "/human_CS7_RCTD_full_top_optimized.pdf"), width = 13, height = 10)
p

dev.off()

saveRDS(
  object = cs7,
  file = "human_cs7_reanno.Rds"
)

























