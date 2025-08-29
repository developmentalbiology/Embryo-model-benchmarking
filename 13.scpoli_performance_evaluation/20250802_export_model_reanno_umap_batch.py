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

# Set up Scanpy visualization parameters
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

def log_message(message):
    """Log message with timestamp"""
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}")

# Your ordered labels list for annotation visualization
ordered_labels = [
    'TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
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
]

def setup_model_configs():
    """
    Setup configuration for all your models
    """
    
    base_path = "/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/final_model"
    
    model_configs = [
        # Lineage models
        {
            'name': 'lineage_hvg2000_dim50_reseed',
            'model_dir': f'{base_path}/lineage_model_hvg2000_dim50_reseed',
            'adata_path': f'{base_path}/lineage_model_hvg2000_dim50_reseed/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_lineage',
            'model_dir': f'{base_path}/enhanced_reference_model_lineage',
            'adata_path': f'{base_path}/enhanced_reference_model_lineage/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_lineage_2ndround',
            'model_dir': f'{base_path}/enhanced_reference_model_lineage_2ndround',
            'adata_path': f'{base_path}/enhanced_reference_model_lineage_2ndround/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method1_balanced_original_lineage',
            'model_dir': f'{base_path}/method1_balanced_original_lineage',
            'adata_path': f'{base_path}/method1_balanced_original_lineage/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method3_balanced_weighted_lineage',
            'model_dir': f'{base_path}/method3_balanced_weighted_lineage',
            'adata_path': f'{base_path}/method3_balanced_weighted_lineage/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        
        # Reanno models
        {
            'name': 'reanno_hvg4000_dim50',
            'model_dir': f'{base_path}/reanno_model_hvg4000_dim50',
            'adata_path': f'{base_path}/reanno_model_hvg4000_dim50/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_reanno',
            'model_dir': f'{base_path}/enhanced_reference_model_reanno',
            'adata_path': f'{base_path}/enhanced_reference_model_reanno/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_reanno_2ndround',
            'model_dir': f'{base_path}/enhanced_reference_model_reanno_2ndround',
            'adata_path': f'{base_path}/enhanced_reference_model_reanno_2ndround/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method1_balanced_original_reanno',
            'model_dir': f'{base_path}/method1_balanced_original_reanno',
            'adata_path': f'{base_path}/method1_balanced_original_reanno/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method3_balanced_weighted_reanno',
            'model_dir': f'{base_path}/method3_balanced_weighted_reanno',
            'adata_path': f'{base_path}/method3_balanced_weighted_reanno/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        }
    ]
    
    return model_configs

def get_unified_label_order(all_labels):
    """
    Get unified label order based on ordered_labels, let scanpy handle colors
    """
    # Create ordered list: first from ordered_labels, then remaining
    final_order = []
    for label in ordered_labels:
        if label in all_labels:
            final_order.append(label)
    
    # Add any labels not in ordered_labels
    for label in all_labels:
        if label not in final_order:
            final_order.append(label)
    
    return final_order

def standardize_label_categories(adata, label_key, unified_order):
    """
    Standardize categorical labels to ensure consistent ordering, let scanpy handle colors
    """
    # Get unique labels in the data
    unique_labels = list(adata.obs[label_key].unique())
    
    # Filter unified_order to only include labels present in this data
    final_order = [label for label in unified_order if label in unique_labels]
    
    # Add any labels in this data that weren't in unified_order
    for label in unique_labels:
        if label not in final_order:
            final_order.append(label)
    
    # Convert to categorical with specified order
    adata.obs[label_key] = pd.Categorical(
        adata.obs[label_key], 
        categories=final_order, 
        ordered=True
    )
    
    log_message(f"Standardized {label_key} categories: {len(final_order)} categories")
    return final_order

def load_and_process_model(config, unified_label_order):
    """
    Load a single model and process it for UMAP visualization
    """
    log_message(f"Processing model: {config['name']}")
    
    try:
        # Load the reference data
        if not os.path.exists(config['adata_path']):
            log_message(f"Warning: Data file not found: {config['adata_path']}")
            return None
            
        source_adata = sc.read_h5ad(config['adata_path'])
        log_message(f"Loaded reference data: {source_adata.shape}")
        
        # Load the model
        if not os.path.exists(config['model_dir']):
            log_message(f"Warning: Model directory not found: {config['model_dir']}")
            return None
            
        scpoli_model = scPoli.load(config['model_dir'], adata=source_adata)
        log_message(f"Model loaded successfully: {config['name']}")
        
        # Get latent representation
        scpoli_model.model.eval()
        data_latent = scpoli_model.get_latent(source_adata, mean=True)
        
        # Add latent representation to adata
        source_adata.obsm["X_scpoli"] = data_latent
        
        # Standardize label categories using unified order
        final_order = standardize_label_categories(
            source_adata, 
            config['label_key'], 
            unified_label_order
        )
        
        # Compute neighbors and UMAP
        sc.pp.neighbors(source_adata, use_rep="X_scpoli")
        sc.tl.umap(source_adata)
        
        return {
            'name': config['name'],
            'adata': source_adata,
            'label_key': config['label_key'],
            'batch_key': config['batch_key']
        }
        
    except Exception as e:
        log_message(f"Error processing model {config['name']}: {str(e)}")
        return None

def generate_umap_plots(processed_models, output_dir="umap_outputs"):
    """
    Generate UMAP plots for all processed models
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    for model_data in processed_models:
        if model_data is None:
            continue
            
        log_message(f"Generating UMAP for {model_data['name']}")
        
        try:
            # Plot UMAP colored by cell type
            plt.figure(figsize=(8, 6))
            sc.pl.umap(
                model_data['adata'], 
                color=model_data['label_key'], 
                title=f"{model_data['name']} - {model_data['label_key']}", 
                frameon=False, 
                show=False,
                legend_loc='right margin',
                save=f"{model_data['name']}_{model_data['label_key']}_umap.pdf"
            )
            plt.close()
            
            # Plot UMAP colored by batch
            plt.figure(figsize=(8, 6))
            sc.pl.umap(
                model_data['adata'], 
                color=model_data['batch_key'], 
                title=f"{model_data['name']} - {model_data['batch_key']}", 
                frameon=False, 
                show=False,
                legend_loc='right margin',
                save=f"{model_data['name']}_{model_data['batch_key']}_umap.pdf"
            )
            plt.close()
            
            log_message(f"UMAP plots saved for {model_data['name']}")
            
        except Exception as e:
            log_message(f"Error generating plots for {model_data['name']}: {str(e)}")

def main():
    """
    Main function to process all models
    """
    log_message("Starting unified model processing and UMAP generation")
    
    # Setup model configurations
    model_configs = setup_model_configs()
    log_message(f"Found {len(model_configs)} models to process")
    
    # First pass: collect all unique labels to create unified order
    all_labels = set()
    for config in model_configs:
        try:
            if os.path.exists(config['adata_path']):
                temp_adata = sc.read_h5ad(config['adata_path'])
                all_labels.update(temp_adata.obs[config['label_key']].unique())
        except Exception as e:
            log_message(f"Warning: Could not read {config['adata_path']} for label collection: {str(e)}")
    
    # Create unified label order (scanpy will handle colors automatically)
    unified_label_order = get_unified_label_order(list(all_labels))
    log_message(f"Created unified label order for {len(all_labels)} unique labels")
    
    # Process all models
    processed_models = []
    for config in model_configs:
        processed_model = load_and_process_model(config, unified_label_order)
        if processed_model:
            processed_models.append(processed_model)
    
    log_message(f"Successfully processed {len(processed_models)} models")
    
    # Generate UMAP plots
    generate_umap_plots(processed_models)
    
    log_message("All processing completed!")
    
    return processed_models

if __name__ == "__main__":
    # Run the main function
    processed_models = main()
    
    # Optional: Print summary
    print("\n" + "="*50)
    print("PROCESSING SUMMARY")
    print("="*50)
    for model_data in processed_models:
        if model_data:
            print(f"✓ {model_data['name']}: {model_data['adata'].shape[0]} cells, {model_data['adata'].shape[1]} genes")
    print("="*50)