#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
from scarches.dataset.trvae.data_handling import remove_sparsity

import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

# Load dataset
os.chdir("/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250729_scpoli_optimization_v3")
adata = sc.read('/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/data/human_clustering_20250724_v3.h5ad')
adata

# use 2000hvg
#sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="cell_ranger", batch_key="orig.ident")

# Reorganize adata with highly variable genes (HVGs)
adata_hvg = adata[:, adata.var.highly_variable].copy()

print("Original obs:", adata_hvg.obs.columns)
print("Original var:", adata_hvg.var.columns)
print("Original uns:", adata_hvg.uns.keys())
print("Original obsm:", adata_hvg.obsm.keys())
print("Original varm:", adata_hvg.varm.keys())
print("Original layers:", adata_hvg.layers.keys())

obs_to_keep = ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'stage','species', 'percent.mt', 'platform', 'reanno', 'lineage']
obs_columns_to_remove = [col for col in adata_hvg.obs.columns if col not in obs_to_keep]
adata_hvg.obs.drop(columns=obs_columns_to_remove, inplace=True)

features_to_keep = ['features']
var_columns_to_remove = [col for col in adata_hvg.var.columns if col not in features_to_keep]
adata_hvg.var.drop(columns=var_columns_to_remove, inplace=True)

adata_hvg.uns = {}
adata_hvg.obsm = {}
adata_hvg.varm = {}

layers_to_keep = ['counts']
layers_keys_to_remove = [key for key in adata_hvg.layers.keys() if key not in layers_to_keep]
for key in layers_keys_to_remove:
    del adata_hvg.layers[key]

print("Modified obs:", adata_hvg.obs.columns)
print("Modified var:", adata_hvg.var.columns)
print("Modified uns:", adata_hvg.uns.keys())
print("Modified obsm:", adata_hvg.obsm.keys())
print("Modified varm:", adata_hvg.varm.keys())
print("Modified layers:", adata_hvg.layers.keys())

counts_matrix = adata_hvg.layers["counts"].toarray()
adata = sc.AnnData(
    X=counts_matrix,
    obs=adata_hvg.obs.copy(),
    var=adata_hvg.var.copy(),
    uns=adata_hvg.uns.copy(),
    obsm=adata_hvg.obsm.copy(),
    varm=adata_hvg.varm.copy(),
    layers={'counts': counts_matrix}
)

adata = remove_sparsity(adata)

# Use adata as the reference
source_adata = adata.copy()

# Train reference model for "final_lineage" (only once)
condition_key = 'orig.ident'
cell_type_key = "lineage"
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

reference_model_dir = './lineage_model_hvg2000_dim50_reseed/'
os.makedirs(reference_model_dir, exist_ok=True)  # Create directory for saving the model
log_dir = os.path.join(reference_model_dir, "logs") 

# Check if the reference model already exists
if not os.path.exists(os.path.join(reference_model_dir, 'model.pt')):
    print("Training reference model...")
    scpoli_model = scPoli(
        adata=source_adata,
        condition_keys=condition_key,
        cell_type_keys=cell_type_key,
        embedding_dims=50, # default=5
        recon_loss='nb',
    )
    scpoli_model.train(
        n_epochs=200, # default=50
        pretraining_epochs=40,
        early_stopping_kwargs=early_stopping_kwargs,
        eta=5,  

    )
    # Save the trained reference model to the directory
    scpoli_model.save(reference_model_dir, overwrite=True, save_anndata=True)
else:
    print("Loading pre-trained reference model...")
    scpoli_model = scPoli.load(reference_model_dir)
    
#get latent representation of reference data
scpoli_model.model.eval()
data_latent_source = scpoli_model.get_latent(
    source_adata,
    mean=True
)

# obtain latent
source_adata.obsm["X_scpoli"] = data_latent_source

# calculate UMAP
sc.pp.neighbors(source_adata, use_rep="X_scpoli")
sc.tl.umap(source_adata)

# 1. plot cell type
plt.figure()
sc.pl.umap(source_adata, color=cell_type_key, title=f'UMAP by {cell_type_key}', frameon=False, show=True, save="lineage_hvg2000_dim50_reseed.pdf")

# 2. plot batch
plt.figure()
sc.pl.umap(source_adata, color=condition_key, title=f'UMAP by {condition_key} (Batch)', frameon=False, show=True,save="lineage_batch_hvg2000_dim50_reseed.pdf")