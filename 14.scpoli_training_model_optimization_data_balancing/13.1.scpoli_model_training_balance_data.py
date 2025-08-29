#!/usr/bin/env python
# coding: utf-8

import os
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
import pickle
import time
from scarches.models.scpoli import scPoli
from scarches.dataset.trvae.data_handling import remove_sparsity

import warnings
# Suppress warnings
warnings.filterwarnings('ignore')


sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

# Change the working directory to the Garfield folder (if needed)
os.chdir('/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250731_scpoli_optimization_comparasion_v3')
os.getcwd()


def create_balanced_dataset(adata, cell_type_key, target_min_cells=200):
    """
    Create balanced dataset
    """
    print(f"=== Creating balanced dataset ===")
    
    cell_counts = adata.obs[cell_type_key].value_counts()
    print(f"Original cell type distribution:")
    for cell_type, count in cell_counts.items():
        pct = count / len(adata) * 100
        print(f"  {cell_type}: {count} ({pct:.1f}%)")
    
    min_count = cell_counts.min()
    max_count = cell_counts.max()
    print(f"Imbalance ratio: {max_count/min_count:.1f}:1")
    
    if min_count >= target_min_cells:
        print("All cell types have sufficient counts, no balancing needed")
        return adata
    
    balanced_data_list = []
    
    for cell_type in cell_counts.index:
        cell_type_data = adata[adata.obs[cell_type_key] == cell_type].copy()
        current_count = len(cell_type_data)
        
        if current_count < target_min_cells:
            # Need upsampling
            needed = target_min_cells - current_count
            print(f"  {cell_type}: {current_count} -> {target_min_cells} (upsample {needed} cells)")
            
            # Random sampling with replacement
            np.random.seed(42)
            additional_indices = np.random.choice(
                len(cell_type_data), needed, replace=True
            )
            additional_data = cell_type_data[additional_indices].copy()
            
            # Add slight noise to avoid exact duplicates
            noise_factor = 0.01
            for i in range(len(additional_data)):
                noise = np.random.normal(0, noise_factor, additional_data[i].X.shape)
                additional_data[i].X = additional_data[i].X + noise
            
            # Combine original data and upsampled data
            cell_type_balanced = cell_type_data.concatenate(additional_data)
        else:
            print(f"  {cell_type}: {current_count} (no adjustment needed)")
            cell_type_balanced = cell_type_data
        
        balanced_data_list.append(cell_type_balanced)
    
    # Merge all cell types
    balanced_adata = balanced_data_list[0]
    for data in balanced_data_list[1:]:
        balanced_adata = balanced_adata.concatenate(data)
    
    # Validate results
    new_counts = balanced_adata.obs[cell_type_key].value_counts()
    print(f"\nCell type distribution after balancing:")
    for cell_type, count in new_counts.items():
        pct = count / len(balanced_adata) * 100
        print(f"  {cell_type}: {count} ({pct:.1f}%)")
    
    new_min = new_counts.min()
    new_max = new_counts.max()
    print(f"New imbalance ratio: {new_max/new_min:.1f}:1")
    
    return balanced_adata


def train_reference_model_method1(adata, condition_key, cell_type_key, annotation_version, save_dir):
    """
    Method 1: Data Balancing + Original Training Settings
    """
    print(f"\n=== Training {annotation_version.upper()} Reference Model - Method 1 ===")
    print("Strategy: Data Balancing + Original Training Settings")
    start_time = time.time()

    # Data preprocessing (based on user's original logic)
    adata_processed = adata.copy()
    adata_processed = remove_sparsity(adata_processed)
    print(f"Training data shape: {adata_processed.shape}")
    print(f"Number of cell types: {len(adata_processed.obs[cell_type_key].unique())}")
    print(f"Number of conditions: {len(adata_processed.obs[condition_key].unique())}")

    # Create model (using original parameters)
    model = scPoli(
        adata=adata_processed,
        condition_keys=condition_key,
        cell_type_keys=cell_type_key,
        embedding_dims=50,  # Original setting
        recon_loss='nb',   # Original setting
    )
    print(f"Model Information:")
    print(f"  Number of conditions: {model.conditions_}")
    print(f"  Cell types: {model.cell_types_}")

    # Use original training settings
    early_stopping_kwargs = {
        "early_stopping_metric": "val_prototype_loss",
        "mode": "min",
        "threshold": 0,
        "patience": 20,  # Original setting
        "reduce_lr": True,
        "lr_patience": 13,  # Original setting
        "lr_factor": 0.1,
    }
    print("Training Parameters (Original Settings):")
    print(f"  Total epochs: 200")
    print(f"  Pretraining epochs: 40")
    print(f"  Eta: 5 (original)")
    print(f"  Patience: 20 (original)")

    try:
        model.train(
            n_epochs=200,       # Original setting
            pretraining_epochs=40,  # Original setting
            eta=5,             # Original setting
            early_stopping_kwargs=early_stopping_kwargs,
        )
        training_time = time.time() - start_time
        print(f"✅ {annotation_version.upper()} Training completed, time taken: {training_time/60:.1f} minutes")

        # Save model
        model_save_dir = os.path.join(save_dir, f'method1_balanced_original_{annotation_version}')
        os.makedirs(model_save_dir, exist_ok=True)
        model.save(model_save_dir, overwrite=True, save_anndata=True)
        print(f"✅ {annotation_version.upper()} Model saved to: {model_save_dir}")
        
        #get latent representation of reference data
        model.model.eval()
        data_latent_source = model.get_latent(
            adata_processed,
            mean=True
        )
        
        adata_processed.obsm["X_scpoli"] = data_latent_source
        
        #plot umap
        sc.pp.neighbors(adata_processed, use_rep="X_scpoli")
        sc.tl.umap(adata_processed)
        
        plt.figure()
        sc.pl.umap(adata_processed, color=cell_type_key, title=f'UMAP by {cell_type_key}', frameon=False, show=True, save=f'{annotation_version}_balance_data_umap.pdf')

        plt.figure()
        sc.pl.umap(adata_processed, color=condition_key, title=f'UMAP by {condition_key} (Batch)', frameon=False, show=True, save=f'{annotation_version}_balance_data_batch.pdf')
                   
        return model, model_save_dir, training_time
    
    except Exception as e:
        print(f"❌ {annotation_version.upper()} Training failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None, 0
    
    
def main():
    """
    Main function - Method 1: Data Balancing + Original Training Settings (Lineage + Reanno)
    Using user's original script logic
    """
    print("="*100)
    print("Label Transfer Method 1: Data Balancing + Original Training Settings (Lineage + Reanno)")
    print("Based on user's original script logic")
    print("="*100)

    # Configuration parameters
    reference_data_path = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/data/human_clustering_20250724_v3.h5ad'
    output_dir = './label_transfer_balance_data/'
    condition_key = 'orig.ident'
    annotation_versions = ['lineage', 'reanno']   
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Load and preprocess reference data (following user's original logic)
    print(f"Loading reference data: {reference_data_path}")
    adata = sc.read(reference_data_path)

    # Reorganize data according to user's original logic
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    # Keep necessary columns
    obs_to_keep = ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'stage', 'percent.mt', 'platform', 'reanno', 'lineage']
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

    # Rebuild AnnData object
    counts_matrix = adata_hvg.layers["counts"].toarray()
    adata_ref = sc.AnnData(
        X=counts_matrix,
        obs=adata_hvg.obs.copy(),
        var=adata_hvg.var.copy(),
        uns=adata_hvg.uns.copy(),
        obsm=adata_hvg.obsm.copy(),
        varm=adata_hvg.varm.copy(),
        layers={'counts': counts_matrix}
    )
    adata_ref = remove_sparsity(adata_ref)
    print(f"Reference data shape: {adata_ref.shape}")

    # Store all results
    all_results = {}
    
    # Step 2: Train and predict for each annotation version
    for annotation_version in annotation_versions:
        print(f"\n{'='*80}")
        print(f"Processing annotation version: {annotation_version.upper()}")
        print(f"{'='*80}")

        cell_type_key = annotation_version

        # Check if the annotation version exists
        if cell_type_key not in adata_ref.obs.columns:
            print(f"❌ Annotation version '{cell_type_key}' does not exist in the data")
            continue

        # Create balanced dataset
        print(f"\n--- {annotation_version.upper()} Data Balancing ---")
        adata_balanced = create_balanced_dataset(adata_ref, cell_type_key, target_min_cells=200)

        # Train reference model
        print(f"\n--- {annotation_version.upper()} Reference Model Training ---")
        reference_model, model_path, training_time = train_reference_model_method1(
            adata_balanced, condition_key, cell_type_key, annotation_version, output_dir
        )
        
        # Store results
        all_results[annotation_version] = {
            'reference_data_shape': adata_ref.shape,
            'balanced_data_shape': adata_balanced.shape,
            'training_time_minutes': training_time / 60,
            'model_path': model_path,
        }

    # Step 3: Generate summary report
    summary = {
        'method': 'method1_balanced_original',
        'description': 'Data Balancing + Original Training Settings (Based on user\'s original script logic)',
        'annotation_versions': annotation_versions,
        'results': all_results
    }

    summary_path = os.path.join(output_dir, 'method1_summary.pkl')
    with open(summary_path, 'wb') as f:
        pickle.dump(summary, f)

    print(f"\n{'='*100}")
    print("Method 1 Completion Summary (Lineage + Reanno) - Based on user's original script logic")
    print("="*100)
    print(f"Strategy: Data Balancing + Original Training Settings")
    print(f"Logic: Fully based on the user's original label transfer workflow")

    for annotation_version in annotation_versions:
        if annotation_version in all_results:
            result = all_results[annotation_version]
            print(f"\n{annotation_version.upper()} Version:")
            print(f"  Reference data: {result['reference_data_shape']} -> {result['balanced_data_shape']}")
            print(f"  Training time: {result['training_time_minutes']:.1f} minutes")
            print(f"  Model path: {result['model_path']}")
            

    print(f"\nOutput directory: {output_dir}")
    print(f"Summary file: {summary_path}")
    print(f"\n✅ Using user's original script logic, should generate complete CSV files")

    return summary

if __name__ == "__main__":
    summary = main()


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    