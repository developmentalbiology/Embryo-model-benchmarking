#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgba
import seaborn as sns
import matplotlib as mpl

# Set up scanpy plotting settings
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 8))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

# Set random seed for reproducibility
sc.settings.seed = 42

# Directories
input_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/output'
figures_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/figures_visualized'

# Create output directory
os.makedirs(figures_folder, exist_ok=True)

# List all annotated h5ad files
h5ad_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('_scPoli_annotated.h5ad')]
print(f"Found {len(h5ad_files)} annotated files")

# Define the ordered labels
lineage_ordered_labels = [
    "Endoderm", "epi", "ExE_endo", "Exe_meso", "hemogenic",
    "mesoderm", "neural_ecto", "non_neuro_ecto", "Notochord",
    "PGC", "TE_TrB", "Gastru"
]

# Your ordered labels list for annotation visualization
ordered_labels = [ 
    'TE','CTBs_1', 'CTBs_2', 'CTBs_3', 'STBs_1', 'STBs_2', 'STBs_3', 'EVTs_1', 'EVTs_2', 'EVTs_3', 'EVTs_4',
    'Epi_1', 'Epi_2', 'Epi_3','Epi_4',
    'Allantois_1', 'Allantois_2', 'pre-YS.mesoderm', 'YS.mesoderm', 'Exe.endothelium', 
    'Amnion', 'Amniotic_epi', 'Ectoderm_1', 'Ectoderm_2',
    'Neural tube', 'Neural crest',
    'Primitive.streak', 'Nascent mesoderm','PGC',
    'Emergent mesoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 'Lateral plate mesoderm_1',
    'Lateral plate mesoderm_2', 'Lateral plate mesoderm_3', 'Lateral plate mesoderm_4',
    'Lateral plate mesoderm_5', 'pre-somatic mesoderm', 'Somite', 'Rostral mesoderm',
    'Cardiac myocyte', 
    'Notochord', 'DE', 'Gut',
    'Hypoblast', 'AVE', 'VE/YE', 'YS.Endoderm_1', 'YS.Endoderm_2', 
    'Hemogenic endothelial progenitors', 'Endothelium', 'Erythroid', 'Myeloid progenitor'
]

# Convert `ordered_labels` list to a set
ordered_labels_set = set(ordered_labels)

# Process each file
for h5ad_file in h5ad_files:
    file_name = os.path.basename(h5ad_file).replace("_scPoli_annotated.h5ad", "")
    print(f"Processing {file_name}...")
    
    try:
        # Load the annotated dataset
        adata = sc.read_h5ad(h5ad_file)
        print(f"Loaded data with shape: {adata.shape}")
        
        # List available columns
        print(f"Available columns in obs: {list(adata.obs.columns)}")
        
        # Process final_lineage_pred
        if 'final_lineage_pred' in adata.obs.columns:
            print(f"Processing lineage predictions, unique values: {adata.obs['final_lineage_pred'].unique()}")
            
            # Ensure the column is treated as categorical with the specified order
            adata.obs["final_lineage_pred"] = pd.Categorical(
            adata.obs["final_lineage_pred"],
            categories=lineage_ordered_labels,
            ordered=True
            )
            
            # Compute UMAP for lineage visualization
            print("Computing UMAP coordinates for lineage visualization")
            
            # If we have X_umap_lineage, use it, otherwise compute UMAP
            if 'X_umap_lineage' in adata.obsm:
                print("Using existing X_umap_lineage coordinates")
                if 'X_umap' not in adata.obsm:
                    adata.obsm['X_umap'] = adata.obsm['X_umap_lineage'].copy()
            else:
                # Normalize and log-transform if not already done
                if 'n_genes' not in adata.obs:
                    sc.pp.filter_cells(adata, min_genes=200)
                    sc.pp.filter_genes(adata, min_cells=3)
                    if 'counts' in adata.layers:
                        adata.X = adata.layers['counts'].copy()
                    sc.pp.normalize_total(adata, target_sum=1e4)
                    sc.pp.log1p(adata)
                
                print("Computing PCA")
                sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="cell_ranger")
                sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
                
                print("Computing neighbors and UMAP")
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)
            
            # Plot UMAP with lineage colors - MANUALLY SAVE FIGURE
            plt.figure(figsize=(12, 10))
            sc.pl.umap(
                adata, 
                color='final_lineage_pred',
                title=f"{file_name} - Lineage Prediction",
                show=False  # Don't show, we'll save manually
            )
            
            outfile = os.path.join(figures_folder, f"{file_name}_lineage_pred.pdf")
            plt.savefig(outfile, dpi=300)
            plt.close()
            print(f"Saved lineage UMAP plot to {outfile}")
        
            
        # Process final_anno_pred
        if 'final_anno_pred' in adata.obs.columns:
            print(f"Processing annotation predictions")
            
            
            # Convert to categorical
            adata.obs['final_anno_pred'] = pd.Categorical(adata.obs['final_anno_pred'], categories=ordered_labels, ordered=True)

            
            # Compute UMAP for annotation visualization
            print("Computing UMAP coordinates for annotation visualization")

            
            # Plot UMAP - MANUALLY SAVE FIGURE
            plt.figure(figsize=(14, 12))
            sc.pl.umap(
                adata, 
                color='final_anno_pred',
                title=f"{file_name} - Cell Type Annotation",
                show=False  # Don't show, we'll save manually
            )
            
            outfile = os.path.join(figures_folder, f"{file_name}_anno_pred.pdf")
            plt.savefig(outfile, dpi=300)
            plt.close()
            print(f"Saved annotation UMAP plot to {outfile}")
            
            
        else:
            print("Warning: 'final_anno_pred' column not found in the dataset")
            

        # Save the processed h5ad file, overwriting the original
        print(f"Saving processed h5ad file back to {h5ad_file}")
        adata.write_h5ad(h5ad_file)
            
        print(f"Completed visualization for {file_name}")
        
    except Exception as e:
        print(f"Error processing {file_name}: {str(e)}")
        import traceback
        print(traceback.format_exc())
        continue

print("All visualization completed!")