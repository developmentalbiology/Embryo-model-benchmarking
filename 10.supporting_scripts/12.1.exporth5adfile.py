import anndata
import scanpy as sc
import numpy as np
from anndata import AnnData
import numpy as np
import pandas as pd
import os

##load dataset
os.chdir("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug")
adata = anndata.read_h5ad("human_monkey_reanno.h5ad")
#adata = anndata.read_h5ad("human_reanno_clean.h5ad")
#adata = anndata.read_h5ad("monkey_reanno_clean.h5ad")

adata

# Extract the indices of highly variable genes
hvg = adata.var[adata.var['highly_variable']].index.tolist()

# Convert to DataFrame for easier export
hvg_df = pd.DataFrame(hvg, columns=['highly_variable_genes'])

# Export to CSV file
hvg_df.to_csv('highly_variable_genes2000.csv', index=False)

from scipy.io import mmwrite
#export adata for seurat object conversion
adata2 = adata.X
adata2_transposed = np.transpose(adata2)
mmwrite("matrix.mtx", adata2_transposed, field="integer")
gene_names = adata.var_names.tolist()
cell_names = adata.obs_names.tolist()
with open("barcodes_new.tsv", "wt") as f:
     f.write("\n".join(cell_names))

with open("features_new.tsv", "wt") as f:
    f.write("\n".join(gene_names))


# Export obs and var information
adata.obs.to_csv('obs_metadata.csv') 
adata.var.to_csv('var_metadata.csv')

# Export the UMAP coordinates
umap_coordinates = adata.obsm['X_umap']

# Convert to DataFrame for better handling and exporting
umap_df = pd.DataFrame(umap_coordinates, index=adata.obs.index, columns=['UMAP1', 'UMAP2'])

# Save the UMAP coordinates to a CSV file
umap_df.to_csv('umap_coordinates.csv')
