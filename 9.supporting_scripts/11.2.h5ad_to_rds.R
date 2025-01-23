library(Seurat)
library(Matrix)
library(readr)


setwd("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug")


mtx <- readMM("matrix.mtx")
bc <- read_tsv("barcodes_new.tsv",col_names = FALSE)
fs <- read_tsv("features_new.tsv",col_names = FALSE)
##check mtx
ncol(mtx) #ncells
nrow(mtx) #ngenes
mtx[1:5, 1:5]

colnames(mtx) <- bc$X1
rownames(mtx) <- fs$X1

##read meta
meta <- read.csv("obs_metadata.csv")

#check duplicated rownames
duplicates <- any(duplicated(meta$X))
if (duplicates) {
  cat("There are duplicate row names in the matrix.\n")
} else {
  cat("Row names are unique.\n")
}


# Identify all duplicated entries
duplicate_ids <- which(duplicated(meta$X) | duplicated(meta$X, fromLast = TRUE))

# Create a table to count occurrences of each duplicate
counts <- table(meta$X[duplicate_ids])

# Loop through and relabel duplicates
for (id in names(counts)) {
  # Get the indices of the duplicates for the current id
  indices <- which(meta$X == id)
  # Rename them with a suffix based on their order
  for (j in seq_along(indices)) {
    meta$X[indices[j]] <- paste0(meta$X[indices[j]], "_", j)
  }
}

# Loop through and relabel duplicates
for (id in names(counts)) {
  # Get the indices of the duplicates for the current id
  indices <- which(bc$X1 == id)
  # Rename them with a suffix based on their order
  for (j in seq_along(indices)) {
    bc$X1[indices[j]] <- paste0(bc$X1[indices[j]], "_", j)
  }
}

#check duplicated rownames
duplicates <- any(duplicated(meta$X))
if (duplicates) {
  cat("There are duplicate row names in the matrix.\n")
} else {
  cat("Row names are unique.\n")
}


colnames(mtx) <- bc$X1

#check duplicates
which(duplicated(colnames(mtx)))
which(duplicated(rownames(meta)))


##seurat object
embryo <- CreateSeuratObject(counts=mtx, project = "embryo",meta.data =meta) 
head(embryo)


# Read the UMAP coordinates
#umap_coords <- read.csv('umap_coordinates.csv', row.names = 1) #error: duplicate 'row.names' are not allowed
umap_coords <- read.csv('umap_coordinates.csv')

# Loop through and relabel duplicates
for (id in names(counts)) {
  # Get the indices of the duplicates for the current id
  indices <- which(umap_coords$X == id)
  # Rename them with a suffix based on their order
  for (j in seq_along(indices)) {
    umap_coords$X[indices[j]] <- paste0(umap_coords$X[indices[j]], "_", j)
  }
}

rownames(umap_coords) <- umap_coords$X
umap_coords <- umap_coords[,-1, drop=FALSE]


##compare cell names
# Extract cell names from both Seurat object and UMAP coordinates
seurat_cells <- colnames(embryo)
umap_cells <- rownames(umap_coords)

# Find cells present in Seurat object but not in UMAP coordinates
missing_in_umap <- setdiff(seurat_cells, umap_cells)

# Find cells present in UMAP coordinates but not in Seurat object
missing_in_seurat <- setdiff(umap_cells, seurat_cells)

# Print missing cells
if (length(missing_in_umap) > 0) {
  cat("Cells present in Seurat object but missing in UMAP coordinates:\n")
  print(missing_in_umap)
} else {
  cat("All cells in Seurat object are present in UMAP coordinates.\n")
}

if (length(missing_in_seurat) > 0) {
  cat("Cells present in UMAP coordinates but missing in Seurat object:\n")
  print(missing_in_seurat)
} else {
  cat("All cells in UMAP coordinates are present in Seurat object.\n")
}


# Ensure the row names of the UMAP coordinates match the Seurat object cell names
#rownames(umap_coords) <- gsub("\\.", "-", rownames(umap_coords))  # if necessary



# Add UMAP coordinates to the Seurat object
embryo[['umap']] <- CreateDimReducObject(embeddings = as.matrix(umap_coords), key = "UMAP_", assay = DefaultAssay(embryo))

# Verify that the UMAP coordinates have been added
DimPlot(embryo, reduction="umap", group.by = "reanno")


saveRDS(embryo, file = "human_monkey_ref.rds")