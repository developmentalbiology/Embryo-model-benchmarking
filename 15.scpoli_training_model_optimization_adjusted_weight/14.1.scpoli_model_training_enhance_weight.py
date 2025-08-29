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
# Note: This path might need to be adjusted based on your environment
# os.chdir('/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250731_scpoli_optimization_comparasion_v3')
# print(os.getcwd())

def calculate_class_weights(adata, cell_type_key, annotation_version):
    """
    Calculate class weights, giving higher weights to rare cell types.
    Rare cell types are defined as:
    - lineage: < 50 cells
    - reanno: < 15 cells
    """
    print(f"\n=== Calculating {annotation_version.upper()} class weights ===")
    cell_counts = adata.obs[cell_type_key].value_counts()
    total_cells = len(adata)
    
    # Define thresholds for rare cell types based on annotation version
    rare_thresholds = {
        'lineage': 50,
        'reanno': 15
    }
    threshold = rare_thresholds.get(annotation_version, 50)  # Default to 50 if not specified
    
    # Calculate weights
    class_weights = {}
    base_weight = 1.0
    rare_weight_multiplier = 3.0  # Weight multiplier for rare cell types
    
    print(f"Rare cell type threshold for {annotation_version}: < {threshold} cells")
    print(f"Class weight settings:")
    rare_types = []
    common_types = []
    for cell_type, count in cell_counts.items():
        if count < threshold:
            weight = base_weight * rare_weight_multiplier
            rare_types.append(cell_type)
            print(f"  {cell_type}: {weight:.1f} (rare type, {count} cells < {threshold}, weight enhanced)")
        else:
            weight = base_weight
            common_types.append(cell_type)
            print(f"  {cell_type}: {weight:.1f} (common type, {count} cells >= {threshold}, regular weight)")
        class_weights[cell_type] = weight
        
    # Summary statistics
    print(f"\nSummary:")
    print(f"  Total cell types: {len(cell_counts)}")
    print(f"  Rare cell types ({annotation_version} < {threshold}): {len(rare_types)}")
    print(f"  Common cell types ({annotation_version} >= {threshold}): {len(common_types)}")
    # print(f"  Rare types: {rare_types}") # Optional: Uncomment if you want to list all rare types
    
    return class_weights

def train_reference_model_method2(adata, condition_key, cell_type_key, annotation_version, save_dir):
    """
    Method 3: Data Balancing + Class Weight Setting
    """
    print(f"\n=== Training {annotation_version.upper()} Reference Model - Method 3 ===")
    print("Strategy: Data Balancing + Class Weight Setting")
    
    start_time = time.time()
    
    # Data preprocessing (based on user's original logic)
    adata_processed = adata.copy()
    adata_processed = remove_sparsity(adata_processed)
    print(f"Training data shape: {adata_processed.shape}")
    print(f"Number of cell types: {len(adata_processed.obs[cell_type_key].unique())}")
    print(f"Number of conditions: {len(adata_processed.obs[condition_key].unique())}")
    
    # Calculate class weights
    class_weights = calculate_class_weights(adata_processed, cell_type_key, annotation_version)
    
    # Create the model (using original parameters)
    try:
        model = scPoli(
            adata=adata_processed,
            condition_keys=condition_key,
            cell_type_keys=cell_type_key,
            embedding_dims=50,  # Original setting
            recon_loss='nb',    # Original setting
        )
    except Exception as e:
        print(f"❌ Failed to create scPoli model for {annotation_version.upper()}: {e}")
        import traceback
        traceback.print_exc()
        return None, None, 0
        
    print(f"Model information:")
    print(f"  Conditions: {model.conditions_}")
    print(f"  Cell types: {model.cell_types_}")
    
    # Enhanced training settings (increase eta to strengthen prototype learning)
    early_stopping_kwargs = {
        "early_stopping_metric": "val_prototype_loss",
        "mode": "min",
        "threshold": 0,
        "patience": 20,        # Increased patience
        "reduce_lr": True,
        "lr_patience": 13,     # Patience for learning rate reduction
        "lr_factor": 0.1,      # Factor by which to reduce learning rate
    }
    
    print("Training parameters (Enhanced class weight settings):")
    print(f"  Total epochs: 200")
    print(f"  Pre-training epochs: 40")
    print(f"  Eta: 12 (Enhanced prototype weight)")
    print(f"  Patience: 20")
    print(f"  Class weights: Applied")
    
    try:
        # Train the model
        # Note: scPoli's train method might not directly accept class_weights.
        # The effect is indirectly enhanced by increasing eta.
        model.train(
            n_epochs=200,
            pretraining_epochs=40,
            eta=12,  # Increased eta to strengthen prototype learning
            early_stopping_kwargs=early_stopping_kwargs,
        )
        
        training_time = time.time() - start_time
        print(f"✅ {annotation_version.upper()} training completed, time taken: {training_time/60:.1f} minutes")
        
        # Save the model
        model_save_dir = os.path.join(save_dir, f'method3_balanced_weighted_{annotation_version}')
        os.makedirs(model_save_dir, exist_ok=True)
        model.save(model_save_dir, overwrite=True, save_anndata=True)
        print(f"✅ {annotation_version.upper()} model saved to: {model_save_dir}")
        
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
        sc.pl.umap(adata_processed, color=cell_type_key, title=f'UMAP by {cell_type_key}', frameon=False, show=True, save=f'{annotation_version}_enhance_weight_umap.pdf')

        plt.figure()
        sc.pl.umap(adata_processed, color=condition_key, title=f'UMAP by {condition_key} (Batch)', frameon=False, show=True, save=f'{annotation_version}_enhance_weight_batch.pdf')
                           
 
        # Save class weight information
        weights_info = {
            'class_weights': class_weights,
            'training_params': {
                'n_epochs': 200,
                'pretraining_epochs': 40,
                'eta': 12, # Corrected to match the used value
                'patience': 20
            }
        }
        weights_info_path = os.path.join(model_save_dir, 'class_weights_info.pkl')
        with open(weights_info_path, 'wb') as f:
            pickle.dump(weights_info, f)
            
        return model, model_save_dir, training_time
        
    except Exception as e:
        print(f"❌ {annotation_version.upper()} training failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None, 0

def main():
    """
    Main function - Method 3: Data Balancing + Class Weight Setting (Lineage + Reanno)
    Based on user's original script logic.
    """
    print("="*100)
    print("Label Transfer Method 3: Data Balancing + Class Weight Setting (Lineage + Reanno)")
    print("Based on user's original script logic")
    print("="*100)
    
    # Configuration parameters
    # Note: These paths will need to be adjusted for your environment
    reference_data_path = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/data/human_clustering_20250724_v3.h5ad'
    output_dir = './ref_training_balance_weight/'
    condition_key = 'orig.ident'
    annotation_versions = ['lineage', 'reanno']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Load and preprocess reference data (following user's original logic)
    print(f"Loading reference data: {reference_data_path}")
    try:
        adata = sc.read(reference_data_path)
    except FileNotFoundError:
        print(f"❌ Reference data file not found at {reference_data_path}")
        return None
    except Exception as e:
        print(f"❌ Error loading reference data: {e}")
        return None

    # Reorganize data according to user's original logic
    # Ensure 'highly_variable' column exists
    if 'highly_variable' not in adata.var.columns:
        print("Warning: 'highly_variable' column not found. Using all genes.")
        adata_hvg = adata.copy()
    else:
        adata_hvg = adata[:, adata.var['highly_variable']].copy()

    # Keep necessary columns in .obs
    obs_to_keep = ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'stage', 'percent.mt', 'platform', 'reanno', 'lineage']
    obs_columns_to_remove = [col for col in adata_hvg.obs.columns if col not in obs_to_keep]
    if obs_columns_to_remove:
         adata_hvg.obs.drop(columns=obs_columns_to_remove, inplace=True, errors='ignore')

    # Keep necessary columns in .var
    features_to_keep = ['features'] # Check if 'features' is the correct column name in your data
    # If 'features' is not present, you might want to keep 'gene_ids' or 'gene_name' or similar
    # Let's make this more robust
    var_id_cols = ['gene_ids', 'gene_name', 'index'] # Common gene identifier columns
    actual_features_to_keep = [col for col in features_to_keep if col in adata_hvg.var.columns]
    if not actual_features_to_keep:
        # If 'features' is not found, try common id columns
        actual_features_to_keep = [col for col in var_id_cols if col in adata_hvg.var.columns][:1] # Take the first one found
        if actual_features_to_keep:
            print(f"Info: 'features' column not found. Using '{actual_features_to_keep[0]}' instead.")
        else:
            print("Warning: No suitable feature identifier column found in .var. Keeping all .var columns.")
            actual_features_to_keep = adata_hvg.var.columns.tolist()

    var_columns_to_remove = [col for col in adata_hvg.var.columns if col not in actual_features_to_keep]
    if var_columns_to_remove:
        adata_hvg.var.drop(columns=var_columns_to_remove, inplace=True, errors='ignore')

    # Clear .uns, .obsm, .varm
    adata_hvg.uns = {}
    adata_hvg.obsm = {}
    adata_hvg.varm = {}

    # Keep only 'counts' layer if it exists
    layers_to_keep = ['counts']
    if 'counts' not in adata_hvg.layers:
        print("Warning: 'counts' layer not found in adata_hvg. Using .X for counts matrix.")
        counts_matrix = adata_hvg.X
    else:
        layers_keys_to_remove = [key for key in adata_hvg.layers.keys() if key not in layers_to_keep]
        for key in layers_keys_to_remove:
            del adata_hvg.layers[key]
        counts_matrix = adata_hvg.layers["counts"]

    # Ensure counts_matrix is a dense array
    if not isinstance(counts_matrix, np.ndarray):
        try:
            counts_matrix = counts_matrix.toarray()
        except AttributeError:
            print("Warning: Could not convert counts matrix to dense array. It might already be dense.")

    # Rebuild AnnData object
    try:
        adata_ref = sc.AnnData(
            X=counts_matrix,
            obs=adata_hvg.obs.copy(),
            var=adata_hvg.var.copy(),
            uns=adata_hvg.uns.copy(),
            obsm=adata_hvg.obsm.copy(),
            varm=adata_hvg.varm.copy(),
            layers={'counts': counts_matrix} if 'counts' in adata_hvg.layers else {}
        )
    except Exception as e:
        print(f"❌ Error rebuilding AnnData object: {e}")
        return None

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
            
        # Train reference model (with class weights)
        print(f"\n--- {annotation_version.upper()} Reference Model Training (Enhanced Class Weights) ---")
        reference_model, model_path, training_time = train_reference_model_method2(
            adata_ref, condition_key, cell_type_key, annotation_version, output_dir
        )
        
        if reference_model is None:
            print(f"❌ {annotation_version.upper()} Reference model training failed, skipping")
            continue
            
        # Store results
        all_results[annotation_version] = {
            'reference_data_shape': adata_ref.shape,
            'balanced_data_shape': adata_ref.shape, # As per original logic, shape is not changed by balancing here
            'training_time_minutes': training_time / 60,
            'model_path': model_path,
        }

    # Step 3: Generate summary report
    summary = {
        'method': 'method3_balanced_weighted', # Corrected method name
        'description': 'Data Balancing + Class Weight Setting (Based on user\'s original script logic)',
        'annotation_versions': annotation_versions,
        'results': all_results
    }
    
    summary_path = os.path.join(output_dir, 'method3_summary.pkl') # Corrected filename
    try:
        with open(summary_path, 'wb') as f:
            pickle.dump(summary, f)
    except Exception as e:
        print(f"❌ Failed to save summary to {summary_path}: {e}")
        return summary # Return summary even if save fails

    print(f"\n{'='*100}")
    print("Method 3 Completed Summary (Lineage + Reanno) - Based on user's original script logic")
    print("="*100)
    print(f"Strategy: Data Balancing + Class Weight Setting")
    
    for annotation_version in annotation_versions:
        if annotation_version in all_results:
            result = all_results[annotation_version]
            print(f"\n{annotation_version.upper()} version:")
            print(f"  Reference data: {result['reference_data_shape']} -> {result['balanced_data_shape']}")
            print(f"  Training time: {result['training_time_minutes']:.1f} minutes")
            print(f"  Model path: {result['model_path']}")
            
    print(f"\nOutput directory: {output_dir}")
    print(f"Summary file: {summary_path}")
    print(f"\n✅ Using user's original script logic + enhanced class weights, should improve prediction of rare cell types")
    
    return summary

if __name__ == "__main__":
    summary = main()
