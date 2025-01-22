import os
import anndata
import scanpy as sc
import doubletdetection
import numpy as np

# Change to the directory where your h5ad files are stored
os.chdir("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/garfield/in_vitro_embryo_models")

# List all h5ad files in the directory
h5ad_files = [f for f in os.listdir() if f.endswith('.h5ad')]

# Loop through each h5ad file
for h5ad_file in h5ad_files:
    print(f"Processing {h5ad_file}...")
    
    try:
        # Load the data
        adata = anndata.read_h5ad(h5ad_file)
        
        # Replace NA values in the observation and variable dataframes
        adata.obs.fillna("NA", inplace=True)
        adata.var.fillna("NA", inplace=True)

        # Remove "empty" genes
        sc.pp.filter_genes(adata, min_cells=1)

        # Convert adata.X to a numpy array if it's a matrix or has 'todense' method
        if hasattr(adata.X, 'todense'):
            X_data = np.asarray(adata.X.todense())  # Convert to dense matrix and then to NumPy array
        elif hasattr(adata.X, 'toarray'):
            X_data = np.asarray(adata.X.toarray())  # Convert to dense matrix and then to NumPy array
        else:
            X_data = np.asarray(adata.X)  # Already dense, convert to NumPy array

        # Check for NaN values
        if np.isnan(X_data).any():
            print(f"NaN values found in {h5ad_file}, skipping...")
            continue

        # Detect doublets
        clf = doubletdetection.BoostClassifier(
            n_iters=40,
            clustering_algorithm="louvain",
            standard_scaling=True,
            pseudocount=0.1,
            n_jobs=-1,
        )
        doublets = clf.fit(X_data).predict(p_thresh=1e-16, voter_thresh=0.5)
        doublet_score = clf.doublet_score()

        # Add doublet information to the AnnData object
        adata.obs["doublet"] = doublets
        adata.obs["doublet_score"] = doublet_score

        # Filter out cells with doublet_score > 100
        adata = adata[adata.obs['doublet_score'] <= 100].copy()

        # Save the processed dataset
        adata.raw.var.rename(columns={'_index': 'index'}, inplace=True)
        
        output_filename = f"processed_{h5ad_file}"
        adata.write_h5ad(filename=output_filename)
        
        print(f"Saved processed data as {output_filename}")

    except Exception as e:
        print(f"An error occurred while processing {h5ad_file}: {e}")
        print(f"Skipping {h5ad_file}.")

print("All files processed.")
