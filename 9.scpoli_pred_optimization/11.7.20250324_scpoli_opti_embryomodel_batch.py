#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
from scarches.models.scpoli import scPoli
from scarches.dataset.trvae.data_handling import remove_sparsity

import warnings
warnings.filterwarnings('ignore')

# Function to print with immediate flushing
def log_message(message):
    print(message, flush=True)
    # Also write to a separate log file as backup
    with open('scpoli_progress.log', 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

log_message("Script started")

# Set up Scanpy visualization parameters
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

# Set random seed for reproducibility
sc.settings.seed = 42

# Set directories
query_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/garfield/in_vitro_embryo_models/processed/data'
output_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/output'
figures_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/figures'

# Create output directories
os.makedirs(output_folder, exist_ok=True)
os.makedirs(figures_folder, exist_ok=True)

log_message(f"Output directories created: {output_folder}, {figures_folder}")

# List all h5ad files in the directory
h5ad_files = [os.path.join(query_folder, f) for f in os.listdir(query_folder) if f.startswith('corrected_processed_')]
log_message(f"Found {len(h5ad_files)} files to process: {[os.path.basename(f) for f in h5ad_files]}")

# Define models paths
lineage_model_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/enhanced_reference_model_lineage/'
anno_model_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/enhanced_reference_model_anno/'

log_message("Loading lineage reference model...")
# Load lineage reference model
adata_lineage_path = os.path.join(lineage_model_dir, "combined_adata.h5ad")
source_adata_lineage = sc.read_h5ad(adata_lineage_path)
log_message(f"Loaded lineage reference data: {source_adata_lineage.shape}")
enhanced_scpoli_lineage_model = scPoli.load(lineage_model_dir, adata=source_adata_lineage)
log_message("Lineage reference model loaded successfully")

log_message("Loading annotation reference model...")
# Load annotation reference model
adata_anno_path = os.path.join(anno_model_dir, "combined_adata.h5ad")
source_adata_anno = sc.read_h5ad(adata_anno_path)
log_message(f"Loaded annotation reference data: {source_adata_anno.shape}")
enhanced_scpoli_anno_model = scPoli.load(anno_model_dir, adata=source_adata_anno)
log_message("Annotation reference model loaded successfully")

# Common settings
condition_key = 'orig.ident'
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

# Function to extend query data with missing genes
def extend_query_data(query_adata, source_var_names):
    log_message(f"Extending query data to match reference genes ({len(source_var_names)} genes)")
    # Create a new AnnData object with genes in the same order as source_adata
    all_genes = source_var_names
    
    # Find common genes and their indices
    common_genes = query_adata.var_names.intersection(all_genes)
    log_message(f"Found {len(common_genes)} common genes between query and reference")
    
    query_gene_indices = [query_adata.var_names.get_loc(gene) for gene in common_genes]
    target_gene_indices = [all_genes.get_loc(gene) for gene in common_genes]
    
    # Create a new matrix with zeros, shape: (n_cells, n_genes_in_source)
    if sp.issparse(query_adata.X):
        log_message("Creating sparse matrix for extended data")
        # For sparse matrix
        X_new = sp.csr_matrix((query_adata.shape[0], len(all_genes)), dtype=np.float32)
        X_query = query_adata.X.tocsr()
        # Fill in values for common genes
        for i, (q_idx, t_idx) in enumerate(zip(query_gene_indices, target_gene_indices)):
            if i % 5000 == 0 and i > 0:
                log_message(f"  Processed {i}/{len(common_genes)} genes...")
            X_new[:, t_idx] = X_query[:, q_idx]
    else:
        log_message("Creating dense matrix for extended data")
        # For dense matrix
        X_new = np.zeros((query_adata.shape[0], len(all_genes)), dtype=np.float32)
        # Fill in values for common genes
        for i, (q_idx, t_idx) in enumerate(zip(query_gene_indices, target_gene_indices)):
            if i % 5000 == 0 and i > 0:
                log_message(f"  Processed {i}/{len(common_genes)} genes...")
            X_new[:, t_idx] = query_adata.X[:, q_idx]
    
    # Create a new AnnData object
    query_adata_extended = sc.AnnData(
        X=X_new,
        obs=query_adata.obs.copy(),
        var=pd.DataFrame(index=all_genes)
    )
    
    # Copy layers
    if 'counts' in query_adata.layers:
        log_message("Copying counts layer to extended data")
        if sp.issparse(query_adata.layers['counts']):
            counts_new = sp.csr_matrix((query_adata.shape[0], len(all_genes)), dtype=np.float32)
            counts_query = query_adata.layers['counts'].tocsr()
            for i, (q_idx, t_idx) in enumerate(zip(query_gene_indices, target_gene_indices)):
                if i % 5000 == 0 and i > 0:
                    log_message(f"  Processed {i}/{len(common_genes)} genes in counts layer...")
                counts_new[:, t_idx] = counts_query[:, q_idx]
        else:
            counts_new = np.zeros((query_adata.shape[0], len(all_genes)), dtype=np.float32)
            for i, (q_idx, t_idx) in enumerate(zip(query_gene_indices, target_gene_indices)):
                if i % 5000 == 0 and i > 0:
                    log_message(f"  Processed {i}/{len(common_genes)} genes in counts layer...")
                counts_new[:, t_idx] = query_adata.layers['counts'][:, q_idx]
        
        query_adata_extended.layers['counts'] = counts_new
    
    # Copy var attributes if needed
    if 'features' in query_adata.var.columns:
        log_message("Copying features column to extended data var")
        var_features = pd.Series(index=all_genes, data=all_genes)
        common_features = query_adata.var['features'][query_gene_indices].values
        var_features.iloc[target_gene_indices] = common_features
        query_adata_extended.var['features'] = var_features
    else:
        query_adata_extended.var['features'] = query_adata_extended.var.index
    
    log_message(f"Data extension completed. Extended shape: {query_adata_extended.shape}")
    return query_adata_extended

# Process each query dataset
for idx, h5ad_file in enumerate(h5ad_files):
    file_name = os.path.basename(h5ad_file).replace(".h5ad", "")
    log_message(f"Processing file {idx+1}/{len(h5ad_files)}: {file_name}")
    
    try:
        # Load query dataset
        log_message(f"Loading query dataset: {h5ad_file}")
        query_adata = sc.read_h5ad(h5ad_file)
        log_message(f"Query dataset loaded. Shape: {query_adata.shape}")
        
        # Preprocess query dataset if needed
        if 'counts' not in query_adata.layers:
            log_message("Adding counts layer to query data")
            query_adata.layers["counts"] = query_adata.X.copy()
        
        # LINEAGE LABEL TRANSFER
        log_message("Starting lineage label transfer")
        
        # Extend query data to match lineage reference
        query_adata_lineage_ext = extend_query_data(query_adata, source_adata_lineage.var_names)
        
        # Ensure orig.ident is set correctly
        if 'orig.ident' not in query_adata_lineage_ext.obs.columns:
            log_message("Setting orig.ident column to 'query'")
            query_adata_lineage_ext.obs['orig.ident'] = 'query'
        else:
            unique_values = query_adata_lineage_ext.obs['orig.ident'].unique()
            if len(unique_values) > 1:
                log_message(f"Warning: Multiple unique values found in 'orig.ident': {unique_values}. Overwriting all to 'query'.")
                query_adata_lineage_ext.obs['orig.ident'] = 'query'
        
        # Set X matrix to float32
        log_message("Converting data to float32")
        if sp.issparse(query_adata_lineage_ext.X):
            query_adata_lineage_ext.X = query_adata_lineage_ext.X.astype(np.float32)
        else:
            query_adata_lineage_ext.X = query_adata_lineage_ext.X.astype(np.float32)
        
        # Set final_lineage to Unknown
        log_message("Setting final_lineage to 'Unknown'")
        query_adata_lineage_ext.obs['final_lineage'] = 'Unknown'
        
        # Create scPoli query object for lineage
        log_message("Initializing scPoli query object for lineage prediction")
        scpoli_lineage_query = scPoli.load_query_data(
            adata=query_adata_lineage_ext,
            reference_model=enhanced_scpoli_lineage_model,
            labeled_indices=[]
        )
        
        # Train the query model
        log_message("Training scPoli lineage query model (this may take some time)...")
        scpoli_lineage_query.train(
            n_epochs=50,
            pretraining_epochs=40,
            eta=10
        )
        log_message("Lineage model training completed")
        
        # Get lineage predictions
        log_message("Getting lineage predictions")
        results_lineage = scpoli_lineage_query.classify(query_adata_lineage_ext, scale_uncertainties=True)
        log_message(f"Lineage predictions completed with {len(results_lineage['final_lineage']['preds'])} predictions")
        
        # Get latent representation
        log_message("Getting latent representation for lineage")
        scpoli_lineage_query.model.eval()
        data_latent_lineage = scpoli_lineage_query.get_latent(
            query_adata_lineage_ext,
            mean=True
        )
        
        # Create AnnData from latent representation
        adata_latent_lineage = sc.AnnData(data_latent_lineage)
        adata_latent_lineage.obs = query_adata_lineage_ext.obs.copy()
        
        # Add predictions to latent data
        log_message("Adding predictions to latent data")
        adata_latent_lineage.obs['final_lineage_pred'] = results_lineage['final_lineage']['preds'].tolist()
        adata_latent_lineage.obs['final_lineage_uncert'] = results_lineage['final_lineage']['uncert'].tolist()
        
        # Compute UMAP for visualization
        log_message("Computing UMAP for lineage visualization")
        sc.pp.neighbors(adata_latent_lineage)
        sc.tl.umap(adata_latent_lineage)
        
        # Transfer predictions to original adata
        log_message("Transferring lineage predictions to original adata")
        matching_indices = adata_latent_lineage.obs.index.intersection(query_adata.obs.index)
        if "final_lineage_pred" not in query_adata.obs.columns:
            query_adata.obs["final_lineage_pred"] = np.nan
        if "final_lineage_uncert" not in query_adata.obs.columns:
            query_adata.obs["final_lineage_uncert"] = np.nan
            
        query_adata.obs.loc[matching_indices, "final_lineage_pred"] = adata_latent_lineage.obs.loc[matching_indices, "final_lineage_pred"]
        query_adata.obs.loc[matching_indices, "final_lineage_uncert"] = adata_latent_lineage.obs.loc[matching_indices, "final_lineage_uncert"]
        
        # Store UMAP coordinates
        log_message("Storing lineage UMAP coordinates")
        query_adata.obsm['X_umap_lineage'] = adata_latent_lineage.obsm['X_umap']
        
        # ANNOTATION LABEL TRANSFER
        log_message("Starting annotation label transfer")
        
        # Extend query data to match annotation reference
        query_adata_anno_ext = extend_query_data(query_adata, source_adata_anno.var_names)
        
        # Ensure orig.ident is set correctly
        if 'orig.ident' not in query_adata_anno_ext.obs.columns:
            query_adata_anno_ext.obs['orig.ident'] = 'query'
        else:
            unique_values = query_adata_anno_ext.obs['orig.ident'].unique()
            if len(unique_values) > 1:
                log_message(f"Warning: Multiple unique values found in 'orig.ident': {unique_values}. Overwriting all to 'query'.")
                query_adata_anno_ext.obs['orig.ident'] = 'query'
        
        # Set X matrix to float32
        if sp.issparse(query_adata_anno_ext.X):
            query_adata_anno_ext.X = query_adata_anno_ext.X.astype(np.float32)
        else:
            query_adata_anno_ext.X = query_adata_anno_ext.X.astype(np.float32)
        
        # Set final_anno to Unknown
        log_message("Setting final_anno to 'Unknown'")
        query_adata_anno_ext.obs['final_anno'] = 'Unknown'
        
        # Create scPoli query object for annotation
        log_message("Initializing scPoli query object for annotation prediction")
        scpoli_anno_query = scPoli.load_query_data(
            adata=query_adata_anno_ext,
            reference_model=enhanced_scpoli_anno_model,
            labeled_indices=[]
        )
        
        # Train the query model
        log_message("Training scPoli annotation query model (this may take some time)...")
        scpoli_anno_query.train(
            n_epochs=50,
            pretraining_epochs=40,
            eta=10
        )
        log_message("Annotation model training completed")
        
        # Get annotation predictions
        log_message("Getting annotation predictions")
        results_anno = scpoli_anno_query.classify(query_adata_anno_ext, scale_uncertainties=True)
        log_message(f"Annotation predictions completed with {len(results_anno['final_anno']['preds'])} predictions")
        
        # Get latent representation
        log_message("Getting latent representation for annotation")
        scpoli_anno_query.model.eval()
        data_latent_anno = scpoli_anno_query.get_latent(
            query_adata_anno_ext,
            mean=True
        )
        
        # Create AnnData from latent representation
        adata_latent_anno = sc.AnnData(data_latent_anno)
        adata_latent_anno.obs = query_adata_anno_ext.obs.copy()
        
        # Add predictions to latent data
        log_message("Adding predictions to latent data")
        adata_latent_anno.obs['final_anno_pred'] = results_anno['final_anno']['preds'].tolist()
        adata_latent_anno.obs['final_anno_uncert'] = results_anno['final_anno']['uncert'].tolist()
        
        # Compute UMAP for visualization
        log_message("Computing UMAP for annotation visualization")
        sc.pp.neighbors(adata_latent_anno)
        sc.tl.umap(adata_latent_anno)
        
        # Transfer predictions to original adata
        log_message("Transferring annotation predictions to original adata")
        matching_indices = adata_latent_anno.obs.index.intersection(query_adata.obs.index)
        if "final_anno_pred" not in query_adata.obs.columns:
            query_adata.obs["final_anno_pred"] = np.nan
        if "final_anno_uncert" not in query_adata.obs.columns:
            query_adata.obs["final_anno_uncert"] = np.nan
            
        query_adata.obs.loc[matching_indices, "final_anno_pred"] = adata_latent_anno.obs.loc[matching_indices, "final_anno_pred"]
        query_adata.obs.loc[matching_indices, "final_anno_uncert"] = adata_latent_anno.obs.loc[matching_indices, "final_anno_uncert"]
        
        # Store UMAP coordinates
        log_message("Storing annotation UMAP coordinates")
        query_adata.obsm['X_umap_anno'] = adata_latent_anno.obsm['X_umap']
        
        # Save UMAP visualizations
        log_message("Saving UMAP visualizations")
        sc.pl.umap(adata_latent_lineage, color=['final_lineage_pred', 'final_lineage_uncert'], 
                  save=f'_{file_name}_lineage.pdf', show=False)
        sc.pl.umap(adata_latent_anno, color=['final_anno_pred', 'final_anno_uncert'], 
                  save=f'_{file_name}_anno.pdf', show=False)
        
        # Save the annotated dataset
        output_file = os.path.join(output_folder, f"{file_name}_scPoli_annotated.h5ad")
        log_message(f"Saving annotated dataset to {output_file}")
        query_adata.write_h5ad(filename=output_file)
        
        # Export metadata as CSV
        metadata_file = os.path.join(output_folder, f"{file_name}_metadata.csv")
        log_message(f"Exporting metadata to {metadata_file}")
        query_adata.obs.to_csv(metadata_file)
        
        log_message(f"Completed processing {file_name}. Saved to {output_file}")
        
    except Exception as e:
        log_message(f"ERROR processing {file_name}: {str(e)}")
        # Print the full error traceback for debugging
        import traceback
        log_message(traceback.format_exc())
        continue

log_message("All datasets processed successfully!")