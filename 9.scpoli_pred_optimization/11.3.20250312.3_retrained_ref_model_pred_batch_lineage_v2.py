#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli

import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

# Change the working directory to the Garfield folder (if needed)
os.chdir('/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim')
os.getcwd()

# load reference data
embryomodel = ad.read_h5ad("/storage/liuxiaodongLab/fanxueying/mayanalysis/embryomodel_integration/embryo_model_integration_scPoli.h5ad")
embryomodel

# filter certain predictions

# Step 1: Identify indices of high-confidence cells
high_confidence_indices = np.where(
    (embryomodel.obs["human_ref_final_lineage_uncert"] < 0.1) | 
    (embryomodel.obs["human_ref_final_anno_pred"] =="PGC")
)[0]
# Step 2: Extract high-confidence cells
high_confidence_cells = embryomodel[high_confidence_indices]

# Step 3: Extract high-confidence labels (cell type annotations)
high_confidence_labels = embryomodel.obs["human_ref_final_lineage_pred"].iloc[high_confidence_indices]


print(high_confidence_cells.X)


from scarches.models.scpoli import scPoli
from scarches.dataset.trvae.data_handling import remove_sparsity

#load reference data
# Load dataset
adata = sc.read('/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/human_reanno_20250108.h5ad')
adata

# Reorganize adata with highly variable genes (HVGs)
adata_hvg = adata[:, adata.var.highly_variable].copy()

print("Original obs:", adata_hvg.obs.columns)
print("Original var:", adata_hvg.var.columns)
print("Original uns:", adata_hvg.uns.keys())
print("Original obsm:", adata_hvg.obsm.keys())
print("Original varm:", adata_hvg.varm.keys())
print("Original layers:", adata_hvg.layers.keys())

obs_to_keep = ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'stage', 'percent.mt', 'species', 'embryo', 'platform', 'ann_level_2', 'ann_level_3', 'ann_level_1', 'final_anno', 'final_lineage']
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
cell_type_key = "final_lineage"
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

reference_model_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/reference_model_lineage/'
os.makedirs(reference_model_dir, exist_ok=True)  # Create directory for saving the model

# Check if the reference model already exists
if not os.path.exists(os.path.join(reference_model_dir, 'model.pt')):
    print("Training reference model...")
    scpoli_model = scPoli(
        adata=source_adata,
        condition_keys=condition_key,
        cell_type_keys=cell_type_key,
        embedding_dims=5,
        recon_loss='nb',
    )
    scpoli_model.train(
        n_epochs=50,
        pretraining_epochs=40,
        early_stopping_kwargs=early_stopping_kwargs,
        eta=5,
    )
    # Save the trained reference model to the directory
    scpoli_model.save(reference_model_dir, overwrite=True)
else:
    print("Loading pre-trained reference model...")
    scpoli_model = scPoli.load(reference_model_dir)


##############################################################################################

## use high_confidence_cells to retrain the reference model
#reorganize query dataset
query_adata = high_confidence_cells

query_adata.layers["counts"] = query_adata.X.copy()
counts_matrix = query_adata.layers["counts"].toarray()
query_adata = sc.AnnData(
            X=counts_matrix,
            obs=query_adata.obs.copy(),
            var=query_adata.var.copy(),
            layers={'counts': counts_matrix}
        )
query_adata = remove_sparsity(query_adata)

# Reorganize query dataset to match genes in the reference dataset
all_genes = source_adata.var_names
missing_genes = all_genes.difference(query_adata.var_names)
missing_data = np.zeros((query_adata.shape[0], len(missing_genes)))
query_adata_df = pd.DataFrame(query_adata.X, columns=query_adata.var_names, index=query_adata.obs_names)
missing_df = pd.DataFrame(missing_data, columns=missing_genes, index=query_adata.obs_names)
query_adata_combined_df = pd.concat([query_adata_df, missing_df], axis=1)[all_genes]
query_adata_extended = sc.AnnData(
            X=query_adata_combined_df.values,
            obs=query_adata.obs,
            var=pd.DataFrame(index=all_genes),
            layers={'counts': query_adata_combined_df.values}
        )

query_adata_extended

source_adata

# Rename the column in query_adata_extended.obs
query_adata_extended.obs.rename(columns={'human_ref_final_lineage_pred': 'final_lineage'}, inplace=True)

# Verify the renaming
print(query_adata_extended.obs.columns)

# Find common columns between query_adata_extended.obs and source_adata.obs
common_columns = query_adata_extended.obs.columns.intersection(source_adata.obs.columns)

# Print the common columns
print("Common columns:")
print(common_columns)

# Retain only the common columns in query_adata_extended.obs
query_adata_extended.obs = query_adata_extended.obs[common_columns]
query_adata_extended

# Step 1: Extract gene names/IDs
query_genes = query_adata_extended.var.index  # Genes from query_adata
source_genes = source_adata.var.index  # Genes from source_adata

# Step 2: Find overlapping genes
common_genes = query_genes.intersection(source_genes)
print(f"Number of common genes: {len(common_genes)}")

# Step 3: Subset both datasets to common genes
query_adata_common = query_adata_extended[:, common_genes].copy()
source_adata_common = source_adata[:, source_adata.var.index.isin(common_genes)].copy()

# Identify all unique columns across both datasets
all_columns = source_adata.obs.columns.union(query_adata_common.obs.columns)

# Add missing columns to query_adata_extended.obs and fill with NaN
for col in all_columns:
    if col not in query_adata_common.obs.columns:
        query_adata_common.obs[col] = np.nan

# Verify the updated columns in query_adata_extended.obs
print("Updated query_adata_common.obs columns:")
print(query_adata_common.obs.columns)

# Concatenate the datasets
combined_adata = sc.concat([source_adata_common, query_adata_common], axis=0, join="outer")
combined_adata

print(combined_adata.obs['final_lineage'])

# Retrain the scPoli model with the combined dataset
enhanced_model_dir='/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/enhanced_reference_model_lineage/'
os.makedirs(enhanced_model_dir, exist_ok=True)  # Create directory for saving the model

print("Retraining scPoli model with enhanced dataset...")
enhanced_scpoli_model = scPoli(
    adata=combined_adata,
    condition_keys=condition_key,
    cell_type_keys=cell_type_key,
    embedding_dims=5,
    recon_loss='nb',
)
enhanced_scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)

# Save the enhanced model and AnnData object
try:
    print("Saving the enhanced model...")
    enhanced_scpoli_model.save(enhanced_model_dir, overwrite=True)
    print(f"Enhanced model saved to {enhanced_model_dir}")
    
    # Save the combined AnnData object
    combined_adata_file = os.path.join(enhanced_model_dir, 'combined_adata.h5ad')
    combined_adata.write(combined_adata_file)
    print(f"Combined AnnData object saved to {combined_adata_file}")
except Exception as e:
    print(f"Error while saving the enhanced model or AnnData: {e}")


##############################################################################################

# Process query datasets from folder
query_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/garfield/in_vitro_embryo_models/processed/data'
figures_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/figures_final_lineage'
os.makedirs(figures_folder, exist_ok=True)  # Create the figures directory if it doesn't exist

# List all h5ad files in the directory
h5ad_files = [os.path.join(query_folder, f) for f in os.listdir(query_folder) if f.startswith('corrected_processed_')]

for h5ad_file in h5ad_files:
    file_name = os.path.basename(h5ad_file)  # Extract just the file name
    print(f"Processing {file_name}...")  # Print only the file name

    try:
        # Load the query dataset
        query_adata = sc.read_h5ad(h5ad_file)

        # Create a set of unique values
        unique_values = set(query_adata.obs["orig.ident"])
        # Ensure there is only one unique value
        if len(unique_values) == 1:
            query_name = unique_values.pop()  # Extract the single value from the set
        else:
            raise ValueError("Expected only one unique value in 'orig.ident', but found multiple or none.")
            
        print(f"Query name: {query_name}")

        # Preprocess query dataset
        sc.settings.seed = 42
        query_adata.layers["counts"] = query_adata.X.copy()
        sc.pp.normalize_total(query_adata, target_sum=1e4)
        sc.pp.log1p(query_adata)
        query_adata.layers["logcounts"] = query_adata.X.copy()
        sc.pp.highly_variable_genes(query_adata, n_top_genes=2000, flavor="cell_ranger", batch_key="orig.ident")
        sc.tl.pca(query_adata, n_comps=30, use_highly_variable=True)
        counts_matrix = query_adata.layers["counts"].toarray()

        query_adata = sc.AnnData(
            X=counts_matrix,
            obs=query_adata.obs.copy(),
            var=query_adata.var.copy(),
            layers={'counts': counts_matrix}
        )

        query_adata = remove_sparsity(query_adata)

        # Reorganize query dataset to match genes in the reference dataset
        all_genes = source_adata.var_names
        missing_genes = all_genes.difference(query_adata.var_names)
        missing_data = np.zeros((query_adata.shape[0], len(missing_genes)))
        query_adata_df = pd.DataFrame(query_adata.X, columns=query_adata.var_names, index=query_adata.obs_names)
        missing_df = pd.DataFrame(missing_data, columns=missing_genes, index=query_adata.obs_names)
        query_adata_combined_df = pd.concat([query_adata_df, missing_df], axis=1)[all_genes]
        query_adata_extended = sc.AnnData(
            X=query_adata_combined_df.values,
            obs=query_adata.obs,
            var=pd.DataFrame(index=all_genes),
            layers={'counts': query_adata_combined_df.values}
        )
        query_adata_extended.var['features'] = query_adata.var.reindex(all_genes)['features']

        # Label transfer to query dataset
        query_adata_extended.obs['final_lineage'] = 'Unknown'
        scpoli_query = scPoli.load_query_data(
            adata=query_adata_extended,
            reference_model=enhanced_scpoli_model,
            labeled_indices=[],
        )
        scpoli_query.train(
            n_epochs=50,
            pretraining_epochs=40,
            eta=10
        )
        query_adata_extended.X = query_adata_extended.X.astype(np.float32)

        #Label transfer from reference to query
        results_dict = scpoli_query.classify(query_adata_extended, scale_uncertainties=True)

        # check the label transfer performance achieved
        for i in range(len(cell_type_key)):
            preds = results_dict[cell_type_key]["preds"]
            results_dict[cell_type_key]["uncert"]
            classification_df = pd.DataFrame(
                classification_report(
                    y_true=query_adata_extended.obs[cell_type_key],
                    y_pred=preds,
                    output_dict=True,
                )
            ).transpose()
        print(classification_df)

        #get latent representation of reference data
        scpoli_query.model.eval()
        data_latent_source = scpoli_query.get_latent(
            source_adata,
            mean=True
        )
        adata_latent_source = sc.AnnData(data_latent_source)
        adata_latent_source.obs = source_adata.obs.copy()

        # Get latent representation of query data
        data_latent = scpoli_query.get_latent(query_adata_extended, mean=True)
        adata_latent = sc.AnnData(data_latent)
        adata_latent.obs = query_adata_extended.obs.copy()

        adata_latent.obs['final_lineage_pred'] = results_dict['final_lineage']['preds'].tolist()
        adata_latent.obs['final_lineage_uncert'] = results_dict['final_lineage']['uncert'].tolist()
        adata_latent.obs['classifier_outcome'] = (
            adata_latent.obs['final_lineage_pred'] == adata_latent.obs['final_lineage']
        )

        # Get prototypes
        labeled_prototypes = scpoli_query.get_prototypes_info()
        labeled_prototypes.obs['study'] = 'labeled prototype'
        unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
        unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

        # Join AnnDatas
        adata_latent_full = adata_latent_source.concatenate(
            [adata_latent, labeled_prototypes, unlabeled_prototypes],
            batch_key='query'
        )
        adata_latent_full.obs['final_lineage_pred'][adata_latent_full.obs['query'].isin(['0'])] = np.nan
        sc.pp.neighbors(adata_latent_full, n_neighbors=15)
        sc.tl.umap(adata_latent_full)

        # Get AnnData without prototypes
        adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]

        # Plot UMAP
        sc.pl.umap(
            adata_no_prototypes,
            color='final_lineage_pred',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_final_lineage.png'
        )
        sc.pl.umap(
            adata_no_prototypes,
            color='orig.ident',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_lineage_dataset.png'
        )
        sc.pl.umap(
            adata_no_prototypes,
            color='final_lineage_uncert',
            show=False,
            frameon=False,
            cmap='magma',
            vmax=1,
            save=f'_{file_name}_scPoli_lineage_uncert.png'
        )

        sc.pp.neighbors(adata_latent)
        sc.tl.leiden(adata_latent)
        sc.tl.umap(adata_latent)
        sc.pl.umap(
            adata_latent,
            color='final_lineage_pred',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_lineage_query.png'
        )

        adata_latent_full
        full_latent = adata_latent_full.obs[adata_latent_full.obs["orig.ident"] == query_name]
        full_latent

        # Save results
        output_csv = os.path.join(figures_folder, f'{file_name}_scPoli_query.csv')
        # Save the .obs DataFrame (metadata) to CSV
        full_latent.to_csv(output_csv, index=True)
        print(f"Saved results to {output_csv}")

    except Exception as e:
        print(f"Error processing {file_name}: {e}")
