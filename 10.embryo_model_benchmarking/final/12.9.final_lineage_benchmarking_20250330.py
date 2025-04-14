#!/usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp
from collections import defaultdict
from glob import glob

# ---------- CELL TYPE TO LINEAGE MAPPING ----------

CELL_TYPE_TO_LINEAGE = {
    # TE_TrB lineage
    'TE': 'TE_TrB',
    'CTBs_1': 'TE_TrB', 
    'CTBs_2': 'TE_TrB', 
    'CTBs_3': 'TE_TrB', 
    'STBs_1': 'TE_TrB', 
    'STBs_2': 'TE_TrB', 
    'STBs_3': 'TE_TrB', 
    'EVTs_1': 'TE_TrB', 
    'EVTs_2': 'TE_TrB', 
    'EVTs_3': 'TE_TrB', 
    'EVTs_4': 'TE_TrB',
    
    # epi lineage
    'Epi_1': 'epi', 
    'Epi_2': 'epi', 
    'Epi_3': 'epi',
    'Epi_4': 'epi',
    
    # Exe_meso lineage
    'Allantois_1': 'Exe_meso', 
    'Allantois_2': 'Exe_meso', 
    'pre-YS.mesoderm': 'Exe_meso', 
    'YS.mesoderm': 'Exe_meso', 
    'Exe.endothelium': 'Exe_meso',
    
    # non_neuro_ecto lineage
    'Amnion': 'non_neuro_ecto', 
    'Amniotic_epi': 'non_neuro_ecto', 
    'Ectoderm_1': 'non_neuro_ecto', 
    'Ectoderm_2': 'non_neuro_ecto',
    
    # neural_ecto lineage
    'Neural tube': 'neural_ecto', 
    'Neural crest': 'neural_ecto',
    
    # Gastru lineage
    'Primitive.streak': 'Gastru', 
    'Nascent mesoderm': 'Gastru',
    
    # PGC lineage
    'PGC': 'PGC',
    
    # mesoderm lineage
    'Emergent mesoderm': 'mesoderm', 
    'Paraxial mesoderm': 'mesoderm', 
    'Intermediate mesoderm': 'mesoderm', 
    'Lateral plate mesoderm_1': 'mesoderm',
    'Lateral plate mesoderm_2': 'mesoderm', 
    'Lateral plate mesoderm_3': 'mesoderm', 
    'Lateral plate mesoderm_4': 'mesoderm',
    'Lateral plate mesoderm_5': 'mesoderm', 
    'pre-somatic mesoderm': 'mesoderm', 
    'Somite': 'mesoderm', 
    'Rostral mesoderm': 'mesoderm',
    'Cardiac myocyte': 'mesoderm',
    
    # Notochord lineage
    'Notochord': 'Notochord',
    
    # Endoderm lineage
    'DE': 'Endoderm', 
    'Gut': 'Endoderm',
    
    # ExE_endo lineage
    'Hypoblast': 'ExE_endo', 
    'AVE': 'ExE_endo', 
    'VE/YE': 'ExE_endo', 
    'YS.Endoderm_1': 'ExE_endo', 
    'YS.Endoderm_2': 'ExE_endo',
    
    # hemogenic lineage
    'Hemogenic endothelial progenitors': 'hemogenic', 
    'Endothelium': 'hemogenic', 
    'Erythroid': 'hemogenic', 
    'Myeloid progenitor': 'hemogenic'
}

# ---------- PATHS ----------

# Define paths to data
DATA_DIR = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/output'
MARKER_OVERLAP_DIR = "/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/marker_overlap_metrics"
EXPRESSION_CORR_DIR = "/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/expression_correlation_metrics"
REFERENCE_FILE = '/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/human_reanno_20250108.h5ad'
OUTPUT_DIR = "./lineage_comparison_results"

# ---------- UTILITY FUNCTIONS ----------

def log_message(message, log_file='lineage_comparison.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")


def map_cell_type_to_lineage(cell_type):
    """Map a cell type to its corresponding lineage."""
    return CELL_TYPE_TO_LINEAGE.get(cell_type, "Unknown")


def preprocess_data(adata):
    """
    Add final_anno_pred_lineage column to the data
    """
    # Make a copy to avoid modifying the original
    adata_copy = adata.copy()
    
    # Add derived lineage column if cell type prediction exists
    if 'final_anno_pred' in adata_copy.obs:
        adata_copy.obs['final_anno_pred_lineage'] = adata_copy.obs['final_anno_pred'].apply(map_cell_type_to_lineage)
        log_message(f"Added final_anno_pred_lineage based on cell type mapping")
    
    return adata_copy


def get_certain_cells_mask(adata, pred_type='final_lineage_pred'):
    """
    Get a boolean mask for certain cells based on prediction type.
    """
    # Default threshold for certainty
    certainty_threshold = 0.2  # Cells with uncertainty < 0.2 are considered certain
    
    if pred_type == 'final_lineage_pred':
        # Check for lineage certainty column
        if 'is_lineage_certain' in adata.obs:
            certain_mask = adata.obs['is_lineage_certain'].values.astype(bool)
            log_message(f"Using 'is_lineage_certain': {np.sum(certain_mask)}/{len(certain_mask)} cells are certain")
        elif 'final_lineage_uncert' in adata.obs:
            certain_mask = adata.obs['final_lineage_uncert'].values < certainty_threshold
            log_message(f"Using 'final_lineage_uncert' threshold {certainty_threshold}: {np.sum(certain_mask)}/{len(certain_mask)} cells are certain")
        else:
            certain_mask = np.ones(adata.n_obs, dtype=bool)
            log_message(f"No lineage certainty found, using all {np.sum(certain_mask)} cells")
            
    else:  # pred_type == 'final_anno_pred_lineage'
        # Check for cell type certainty column
        if 'is_celltype_certain' in adata.obs:
            certain_mask = adata.obs['is_celltype_certain'].values.astype(bool)
            log_message(f"Using 'is_celltype_certain': {np.sum(certain_mask)}/{len(certain_mask)} cells are certain")
        elif 'final_anno_uncert' in adata.obs:
            certain_mask = adata.obs['final_anno_uncert'].values < certainty_threshold
            log_message(f"Using 'final_anno_uncert' threshold {certainty_threshold}: {np.sum(certain_mask)}/{len(certain_mask)} cells are certain")
        else:
            certain_mask = np.ones(adata.n_obs, dtype=bool)
            log_message(f"No cell type certainty found, using all {np.sum(certain_mask)} cells")
    
    return certain_mask

# ---------- METRICS CALCULATION FUNCTIONS ----------

def calculate_lineage_consistency(adata):
    """
    Calculate consistency between direct lineage predictions and 
    lineages derived from cell type predictions.
    """
    # Ensure necessary columns exist
    required_cols = ['final_lineage_pred', 'final_anno_pred_lineage']
    if not all(col in adata.obs.columns for col in required_cols):
        log_message(f"Missing required columns for consistency calculation: {[col for col in required_cols if col not in adata.obs.columns]}")
        return {'error': 'Missing columns for calculation'}
    
    # Calculate overall consistency
    matching_cells = (adata.obs['final_lineage_pred'] == adata.obs['final_anno_pred_lineage']).sum()
    total_cells = adata.n_obs
    overall_consistency = matching_cells / total_cells if total_cells > 0 else 0
    
    # Calculate consistency for each lineage
    lineage_consistency = {}
    all_lineages = set(CELL_TYPE_TO_LINEAGE.values())
    
    for lineage in all_lineages:
        # Get cells with this lineage in either prediction
        lineage_mask = ((adata.obs['final_lineage_pred'] == lineage) | 
                     (adata.obs['final_anno_pred_lineage'] == lineage))
        
        lineage_cells = np.sum(lineage_mask)
        
        if lineage_cells > 0:
            matching = np.sum((adata.obs.loc[lineage_mask, 'final_lineage_pred'] == 
                           adata.obs.loc[lineage_mask, 'final_anno_pred_lineage']))
            lineage_consistency[lineage] = matching / lineage_cells
        else:
            lineage_consistency[lineage] = np.nan
    
    return {
        'overall_consistency': overall_consistency,
        'lineage_specific_consistency': lineage_consistency
    }


def calculate_lineage_abundance(adata, lineage, pred_type='final_lineage_pred'):
    """
    Calculate abundance (proportion of cells) for a lineage using only certain cells.
    """
    # Get certain cells mask
    certain_mask = get_certain_cells_mask(adata, pred_type)
    
    # Get adata with only certain cells
    adata_certain = adata[certain_mask]
    
    # Total number of certain cells
    total_certain_cells = adata_certain.n_obs
    
    if total_certain_cells == 0:
        return 0
    
    # Get lineage mask
    if pred_type == 'final_lineage_pred':
        lineage_mask = adata_certain.obs['final_lineage_pred'] == lineage
    else:  # final_anno_pred_lineage
        lineage_mask = adata_certain.obs['final_anno_pred_lineage'] == lineage
    
    # Calculate abundance
    lineage_cells = np.sum(lineage_mask)
    abundance = lineage_cells / total_certain_cells
    
    return abundance


def calculate_percentage_of_certain_cells(adata, lineage, pred_type='final_lineage_pred'):
    """
    Calculate percentage of cells with high prediction certainty for a lineage.
    """
    # Get certainty threshold 
    certainty_threshold = 0.2
    
    if pred_type == 'final_lineage_pred':
        # Get cells belonging to this lineage
        lineage_mask = adata.obs['final_lineage_pred'] == lineage
        
        # Lineage cells
        lineage_cells = np.sum(lineage_mask)
        if lineage_cells == 0:
            return 0
        
        # Check for certainty column
        if 'final_lineage_uncert' in adata.obs:
            uncertainty_values = adata.obs.loc[lineage_mask, 'final_lineage_uncert'].values
            certain_cells = np.sum(uncertainty_values < certainty_threshold)
            return certain_cells / lineage_cells
        
        elif 'is_lineage_certain' in adata.obs:
            certainty_values = adata.obs.loc[lineage_mask, 'is_lineage_certain'].values
            certain_cells = np.sum(certainty_values)
            return certain_cells / lineage_cells
    
    else:  # final_anno_pred_lineage
        # For cell type-derived lineage
        lineage_mask = adata.obs['final_anno_pred_lineage'] == lineage
        
        # Lineage cells
        lineage_cells = np.sum(lineage_mask)
        if lineage_cells == 0:
            return 0
            
        if 'final_anno_uncert' in adata.obs:
            uncertainty_values = adata.obs.loc[lineage_mask, 'final_anno_uncert'].values
            certain_cells = np.sum(uncertainty_values < certainty_threshold)
            return certain_cells / lineage_cells
            
        elif 'is_celltype_certain' in adata.obs:
            certainty_values = adata.obs.loc[lineage_mask, 'is_celltype_certain'].values
            certain_cells = np.sum(certainty_values)
            return certain_cells / lineage_cells
    
    return 0  # Default if we couldn't calculate


def calculate_mean_certainty(adata, lineage, pred_type='final_lineage_pred'):
    """
    Calculate mean prediction certainty for a lineage using only certain cells.
    """
    # Get certainty mask
    certain_mask = get_certain_cells_mask(adata, pred_type)
    
    if pred_type == 'final_lineage_pred':
        # Combine lineage and certainty masks
        lineage_mask = adata.obs['final_lineage_pred'] == lineage
        combined_mask = lineage_mask & certain_mask
        
        if np.sum(combined_mask) == 0:
            return 0
        
        # Calculate mean certainty for certain cells
        if 'final_lineage_uncert' in adata.obs:
            uncertainty_values = adata.obs.loc[combined_mask, 'final_lineage_uncert'].values
            return 1 - np.mean(uncertainty_values)  # Convert uncertainty to certainty
    else:
        # For cell type-derived lineage
        lineage_mask = adata.obs['final_anno_pred_lineage'] == lineage
        combined_mask = lineage_mask & certain_mask
        
        if np.sum(combined_mask) == 0:
            return 0
            
        if 'final_anno_uncert' in adata.obs:
            uncertainty_values = adata.obs.loc[combined_mask, 'final_anno_uncert'].values
            return 1 - np.mean(uncertainty_values)
    
    return 0  # Default if we couldn't calculate


def calculate_cell_type_presence(adata, lineage, pred_type='final_lineage_pred', ref_adata=None):
    """
    Calculate cell type presence score for a lineage using only certain cells.
    """
    # Get certainty mask
    certain_mask = get_certain_cells_mask(adata, pred_type)
    
    # Get adata with only certain cells
    adata_certain = adata[certain_mask]
    
    # Create reverse mapping from lineage to expected cell types
    lineage_to_cell_types = defaultdict(list)
    for cell_type, lin in CELL_TYPE_TO_LINEAGE.items():
        lineage_to_cell_types[lin].append(cell_type)
    
    # Check if we should use reference data for mapping
    if ref_adata is not None and 'final_anno' in ref_adata.obs.columns:
        # Build mapping from reference data
        ref_cell_types = sorted(ref_adata.obs['final_anno'].unique().tolist())
        lineage_to_cell_types = defaultdict(list)
        for cell_type in ref_cell_types:
            lin = map_cell_type_to_lineage(cell_type)
            if lin != "Unknown":
                lineage_to_cell_types[lin].append(cell_type)
    
    # Get expected cell types for this lineage
    expected_cell_types = lineage_to_cell_types.get(lineage, [])
    
    # Skip if no expected cell types
    if not expected_cell_types:
        return 0
    
    # Get cells of this lineage
    if pred_type == 'final_lineage_pred':
        lineage_mask = adata_certain.obs['final_lineage_pred'] == lineage
    else:  # final_anno_pred_lineage
        lineage_mask = adata_certain.obs['final_anno_pred_lineage'] == lineage
    
    # Skip if no cells of this lineage
    if np.sum(lineage_mask) == 0:
        return 0
    
    # Get unique cell types for this lineage
    observed_types = set(adata_certain.obs.loc[lineage_mask, 'final_anno_pred'].unique())
    
    # Calculate presence score
    present_types = sum(1 for ct in expected_cell_types if ct in observed_types)
    presence_score = present_types / len(expected_cell_types)
    
    return presence_score

# ---------- LOADING PRE-CALCULATED METRICS ----------

def extract_model_name(filename):
    """Extract model name from filename"""
    # Remove file extension
    base_name = os.path.basename(filename)
    
    # Remove the suffix part
    if "_celltype_marker_overlap.csv" in base_name:
        model_name = base_name.replace("_celltype_marker_overlap.csv", "")
    elif "_lineage_marker_overlap.csv" in base_name:
        model_name = base_name.replace("_lineage_marker_overlap.csv", "")
    elif "_cell_type_pearson_correlation.csv" in base_name:
        model_name = base_name.replace("_cell_type_pearson_correlation.csv", "")
    elif "_lineage_pearson_correlation.csv" in base_name:
        model_name = base_name.replace("_lineage_pearson_correlation.csv", "")
    elif "_scPoli_annotated.h5ad" in base_name:
        model_name = base_name.replace(".h5ad", "")
    else:
        model_name = base_name
    
    return model_name


def load_marker_overlap_files(directory):
    """
    Load marker overlap files (both cell type and lineage) from directory
    
    Returns:
    --------
    Dict with model names as keys, each containing 'celltype' and 'lineage' DataFrames
    """
    marker_files = {}
    
    # Find all cell type marker overlap files
    celltype_pattern = "*_celltype_marker_overlap.csv"
    celltype_files = glob(os.path.join(directory, celltype_pattern))
    
    for file_path in celltype_files:
        model_name = extract_model_name(file_path)
        try:
            df = pd.read_csv(file_path)
            
            if model_name not in marker_files:
                marker_files[model_name] = {}
                
            marker_files[model_name]['celltype'] = df
            log_message(f"Loaded cell type marker data for {model_name}: {len(df)} cell types")
        except Exception as e:
            log_message(f"Error loading {file_path}: {str(e)}")
    
    # Find all lineage marker overlap files
    lineage_pattern = "*_lineage_marker_overlap.csv"
    lineage_files = glob(os.path.join(directory, lineage_pattern))
    
    for file_path in lineage_files:
        model_name = extract_model_name(file_path)
        try:
            df = pd.read_csv(file_path)
            
            if model_name not in marker_files:
                marker_files[model_name] = {}
                
            marker_files[model_name]['lineage'] = df
            log_message(f"Loaded lineage marker data for {model_name}: {len(df)} lineages")
        except Exception as e:
            log_message(f"Error loading {file_path}: {str(e)}")
            
    return marker_files


def load_expression_correlation_files(directory):
    """
    Load expression correlation files (both cell type and lineage) from directory
    
    Returns:
    --------
    Dict with model names as keys, each containing 'celltype' and 'lineage' DataFrames
    """
    correlation_files = {}
    
    # Find all cell type correlation files
    celltype_pattern = "*_cell_type_pearson_correlation.csv"
    celltype_files = glob(os.path.join(directory, celltype_pattern))
    
    for file_path in celltype_files:
        model_name = extract_model_name(file_path)
        try:
            df = pd.read_csv(file_path)
            
            if model_name not in correlation_files:
                correlation_files[model_name] = {}
                
            correlation_files[model_name]['celltype'] = df
            log_message(f"Loaded cell type correlation data for {model_name}: {len(df)} cell types")
        except Exception as e:
            log_message(f"Error loading {file_path}: {str(e)}")
    
    # Find all lineage correlation files
    lineage_pattern = "*_lineage_pearson_correlation.csv"
    lineage_files = glob(os.path.join(directory, lineage_pattern))
    
    for file_path in lineage_files:
        model_name = extract_model_name(file_path)
        try:
            df = pd.read_csv(file_path)
            
            if model_name not in correlation_files:
                correlation_files[model_name] = {}
                
            correlation_files[model_name]['lineage'] = df
            log_message(f"Loaded lineage correlation data for {model_name}: {len(df)} lineages")
        except Exception as e:
            log_message(f"Error loading {file_path}: {str(e)}")
            
    return correlation_files


def group_celltype_metrics_by_lineage(celltype_df):
    """
    Group cell type metrics by their derived lineage
    """
    if celltype_df.empty:
        return pd.DataFrame()
        
    # Create a copy of cell type dataframe
    df = celltype_df.copy()
    
    # Map cell types to lineages
    df['Lineage'] = df['Group'].apply(map_cell_type_to_lineage)
    
    # Group by lineage and calculate mean metrics
    lineage_metrics = []
    for lineage, group in df.groupby('Lineage'):
        if lineage == "Unknown":
            continue
            
        # Get all metrics
        metrics = {
            'Lineage': lineage,
            'Cell_Types': len(group),
            'Mean_Precision': group['Precision'].mean(),
            'Mean_Recall': group['Recall'].mean(),
            'Mean_F1_Score': group['F1_Score'].mean(),
            'Mean_Jaccard': group['Jaccard'].mean(),
            'Total_Overlap': group['Overlap'].sum(),
            'Total_Query_Markers': group['Query_Markers'].sum(),
            'Total_Ref_Markers': group['Ref_Markers'].sum()
        }
        
        lineage_metrics.append(metrics)
    
    # Create lineage metrics dataframe
    lineage_df = pd.DataFrame(lineage_metrics)
    
    # Calculate overall precision and recall
    if not lineage_df.empty and 'Total_Query_Markers' in lineage_df.columns:
        lineage_df['Overall_Precision'] = lineage_df['Total_Overlap'] / lineage_df['Total_Query_Markers']
        lineage_df['Overall_Recall'] = lineage_df['Total_Overlap'] / lineage_df['Total_Ref_Markers']
    
    return lineage_df


def group_celltype_correlation_by_lineage(celltype_corr_df):
    """
    Group cell type correlation by their derived lineage
    """
    if celltype_corr_df.empty:
        return pd.DataFrame()
        
    # Create a copy of cell type dataframe
    df = celltype_corr_df.copy()
    
    # Map cell types to lineages
    df['Lineage'] = df['Label'].apply(map_cell_type_to_lineage)
    
    # Group by lineage and calculate mean correlation
    lineage_metrics = []
    for lineage, group in df.groupby('Lineage'):
        if lineage == "Unknown":
            continue
            
        # Get mean correlation
        metrics = {
            'Lineage': lineage,
            'Cell_Types': len(group),
            'Mean_Pearson_r': group['Pearson_r'].mean(),
            'Mean_Pearson_p': group['Pearson_p'].mean(),
            'Total_Query_Cells': group['Query_Cells'].sum()
        }
        
        lineage_metrics.append(metrics)
    
    # Create lineage metrics dataframe
    lineage_df = pd.DataFrame(lineage_metrics)
    
    return lineage_df


def calculate_composite_score(metrics_dict):
    """
    Calculate a composite score from multiple metrics.
    
    Parameters:
    -----------
    metrics_dict : dict
        Dictionary with metric names and values
    """
    # Define weights for different metrics
    weights = {
        'Precision': 0.10, 
        'Recall': 0.10,
        'F1_Score': 0.10,
        'Jaccard': 0.10,
        'Pearson_r': 0.10,
        'consistency': 0.10,
        'abundance': 0.10,
        'presence': 0.10,
        'percentage_certain_cells': 0.10,
        'mean_certainty': 0.10
    }
    
    # Calculate weighted average
    total_weight = 0
    weighted_sum = 0
    
    for metric, weight in weights.items():
        if metric in metrics_dict and not pd.isna(metrics_dict[metric]):
            value = metrics_dict[metric]
            # Ensure value is valid
            if isinstance(value, (int, float)) and np.isfinite(value):
                weighted_sum += value * weight
                total_weight += weight
    
    # Return composite score
    if total_weight == 0:
        return 0
    return weighted_sum / total_weight


def create_radar_plot(direct_metrics, derived_metrics, lineage, model_name, output_dir):
    """
    Create a radar plot comparing direct lineage metrics with derived lineage metrics
    """
    # Create output directory
    os.makedirs(os.path.join(output_dir, 'radar_plots'), exist_ok=True)
    
    # Convert to dictionaries if Series
    if isinstance(direct_metrics, pd.Series):
        direct_metrics = direct_metrics.to_dict()
    if isinstance(derived_metrics, pd.Series):
        derived_metrics = derived_metrics.to_dict()
    
    # Define metrics to include in radar plot
    metrics_to_plot = [
        'Precision', 'Recall', 'F1_Score', 'Jaccard', 'Pearson_r',
        'consistency', 'abundance', 'presence', 
        'percentage_certain_cells', 'mean_certainty'
    ]
    
    # Define labels for plotting
    metric_labels = [m.replace('_', ' ').capitalize() for m in metrics_to_plot]
    
    # Get metrics mapping between direct and derived
    direct_to_derived = {
        'Precision': 'Mean_Precision', 
        'Recall': 'Mean_Recall', 
        'F1_Score': 'Mean_F1_Score', 
        'Jaccard': 'Mean_Jaccard',
        'Pearson_r': 'Mean_Pearson_r',
        'consistency': 'consistency',
        'abundance': 'abundance',
        'presence': 'presence',
        'percentage_certain_cells': 'percentage_certain_cells',
        'mean_certainty': 'mean_certainty'
    }
    
    # Set up the plot
    plt.figure(figsize=(12, 10))
    ax = plt.subplot(111, polar=True)
    
    # Set angles for each metric
    angles = np.linspace(0, 2 * np.pi, len(metrics_to_plot), endpoint=False).tolist()
    angles += angles[:1]  # Close the loop
    
    # Get values for direct metrics
    direct_values = []
    for metric in metrics_to_plot:
        if metric in direct_metrics:
            val = direct_metrics[metric]
        else:
            val = 0
        direct_values.append(val)
    direct_values += direct_values[:1]  # Close the loop
    
    # Get values for derived metrics
    derived_values = []
    for metric in metrics_to_plot:
        derived_key = direct_to_derived.get(metric, metric)
        if derived_key in derived_metrics:
            val = derived_metrics[derived_key]
        else:
            val = 0
        derived_values.append(val)
    derived_values += derived_values[:1]  # Close the loop
    
    # Plot the direct metrics
    ax.plot(angles, direct_values, 'o-', linewidth=2, 
          label=f"Direct (final_lineage_pred)", color='royalblue')
    ax.fill(angles, direct_values, alpha=0.1, color='royalblue')
    
    # Plot the derived metrics
    ax.plot(angles, derived_values, 'o-', linewidth=2, 
          label=f"Derived (cell types by lineage)", color='darkgreen')
    ax.fill(angles, derived_values, alpha=0.1, color='darkgreen')
    
    # Calculate composite scores
    direct_composite = calculate_composite_score(direct_metrics)
    derived_composite = calculate_composite_score({
        k: derived_metrics.get(direct_to_derived.get(k, k), 0) 
        for k in metrics_to_plot
    })
    
    # Set chart properties
    plt.xticks(angles[:-1], metric_labels)
    ax.set_ylim(0, 1)
    plt.yticks([0.2, 0.4, 0.6, 0.8], ["0.2", "0.4", "0.6", "0.8"], fontsize=8)
    
    # Add title with composite scores
    plt.title(f"{lineage} - All Metrics Comparison\n{model_name}\nComposite: Direct={direct_composite:.2f}, Derived={derived_composite:.2f}")
    
    # Add number of cell types to legend if available
    cell_types_count = derived_metrics.get('Cell_Types', 0)
    if cell_types_count > 0:
        derived_label = f"Derived (avg of {cell_types_count} cell types)"
    else:
        derived_label = "Derived (cell types by lineage)"
        
    # Add custom legend
    plt.legend([
        f"Direct (final_lineage_pred)",
        derived_label
    ], loc='upper right', bbox_to_anchor=(0.1, 0.1))
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'radar_plots', f"{model_name}_{lineage}_radar.png"), 
              dpi=300, bbox_inches='tight')
    plt.close()


def combine_all_lineage_metrics(direct_lineage_df, derived_lineage_df, lineage_corr_df, 
                              celltype_corr_by_lineage_df, direct_metrics=None, derived_metrics=None):
    """
    Combine direct and derived lineage metrics from all sources into a comprehensive DataFrame
    
    Parameters:
    -----------
    direct_lineage_df : pd.DataFrame
        DataFrame with direct lineage marker metrics
    derived_lineage_df : pd.DataFrame
        DataFrame with derived lineage marker metrics
    lineage_corr_df : pd.DataFrame
        DataFrame with direct lineage correlation metrics
    celltype_corr_by_lineage_df : pd.DataFrame
        DataFrame with derived lineage correlation metrics
    direct_metrics : dict, optional
        Dictionary with additional direct lineage metrics (consistency, abundance, etc.)
    derived_metrics : dict, optional
        Dictionary with additional derived lineage metrics
        
    Returns:
    --------
    pd.DataFrame
        Combined metrics for all lineages
    """
    # Start with all lineages from all sources
    all_lineages = set()
    
    # Add lineages from marker overlap data
    if not direct_lineage_df.empty and 'Group' in direct_lineage_df.columns:
        all_lineages.update(direct_lineage_df['Group'])
    
    if not derived_lineage_df.empty and 'Lineage' in derived_lineage_df.columns:
        all_lineages.update(derived_lineage_df['Lineage'])
        
    # Add lineages from additional direct metrics
    if direct_metrics is not None:
        all_lineages.update(direct_metrics.keys())
    
    # Create combined metrics
    combined_metrics = []
    
    for lineage in all_lineages:
        # Direct lineage metrics
        direct_metrics_row = {
            'Lineage': lineage,
            'Type': 'Direct'
        }
        
        # Get marker metrics for direct lineage
        if not direct_lineage_df.empty and 'Group' in direct_lineage_df.columns:
            direct_marker = direct_lineage_df[direct_lineage_df['Group'] == lineage]
            if not direct_marker.empty:
                direct_row = direct_marker.iloc[0]
                direct_metrics_row.update({
                    'Precision': direct_row['Precision'],
                    'Recall': direct_row['Recall'],
                    'F1_Score': direct_row['F1_Score'],
                    'Jaccard': direct_row['Jaccard'],
                    'Query_Markers': direct_row['Query_Markers'],
                    'Ref_Markers': direct_row['Ref_Markers'],
                    'Overlap': direct_row['Overlap']
                })
        
        # Get correlation metrics for direct lineage
        if not lineage_corr_df.empty and 'Label' in lineage_corr_df.columns:
            direct_corr = lineage_corr_df[lineage_corr_df['Label'] == lineage]
            if not direct_corr.empty:
                direct_metrics_row['Pearson_r'] = direct_corr.iloc[0]['Pearson_r']
                direct_metrics_row['Pearson_p'] = direct_corr.iloc[0]['Pearson_p']
                direct_metrics_row['Query_Cells'] = direct_corr.iloc[0]['Query_Cells']
                
        # Add additional direct metrics if available
        if direct_metrics is not None and lineage in direct_metrics:
            lineage_direct_metrics = direct_metrics[lineage]
            for key, value in lineage_direct_metrics.items():
                direct_metrics_row[key] = value
        
        # Calculate composite score
        direct_metrics_row['Composite_Score'] = calculate_composite_score(direct_metrics_row)
        
        # Add direct metrics to results
        combined_metrics.append(direct_metrics_row)
        
        # Derived lineage metrics
        derived_metrics_row = {
            'Lineage': lineage,
            'Type': 'Derived'
        }
        
        # Get marker metrics for derived lineage (cell types grouped by lineage)
        if not derived_lineage_df.empty and 'Lineage' in derived_lineage_df.columns:
            derived_marker = derived_lineage_df[derived_lineage_df['Lineage'] == lineage]
            if not derived_marker.empty:
                derived_row = derived_marker.iloc[0]
                derived_metrics_row.update({
                    'Mean_Precision': derived_row['Mean_Precision'],
                    'Mean_Recall': derived_row['Mean_Recall'],
                    'Mean_F1_Score': derived_row['Mean_F1_Score'],
                    'Mean_Jaccard': derived_row['Mean_Jaccard'],
                    'Cell_Types': derived_row['Cell_Types'],
                    'Total_Query_Markers': derived_row['Total_Query_Markers'],
                    'Total_Ref_Markers': derived_row['Total_Ref_Markers'],
                    'Total_Overlap': derived_row['Total_Overlap']
                })
        
        # Get correlation metrics for derived lineage (cell types grouped by lineage)
        if not celltype_corr_by_lineage_df.empty and 'Lineage' in celltype_corr_by_lineage_df.columns:
            derived_corr = celltype_corr_by_lineage_df[celltype_corr_by_lineage_df['Lineage'] == lineage]
            if not derived_corr.empty:
                derived_metrics_row['Mean_Pearson_r'] = derived_corr.iloc[0]['Mean_Pearson_r']
                derived_metrics_row['Mean_Pearson_p'] = derived_corr.iloc[0]['Mean_Pearson_p']
                derived_metrics_row['Cell_Types_With_Corr'] = derived_corr.iloc[0]['Cell_Types']
        
        # Add additional derived metrics if available
        if derived_metrics is not None and lineage in derived_metrics:
            lineage_derived_metrics = derived_metrics[lineage]
            for key, value in lineage_derived_metrics.items():
                derived_metrics_row[key] = value
                
        # Calculate composite score based on mean metrics
        derived_composite_inputs = {
            'Precision': derived_metrics_row.get('Mean_Precision', 0),
            'Recall': derived_metrics_row.get('Mean_Recall', 0),
            'F1_Score': derived_metrics_row.get('Mean_F1_Score', 0),
            'Jaccard': derived_metrics_row.get('Mean_Jaccard', 0),
            'Pearson_r': derived_metrics_row.get('Mean_Pearson_r', 0),
            'consistency': derived_metrics_row.get('consistency', 0),
            'abundance': derived_metrics_row.get('abundance', 0),
            'presence': derived_metrics_row.get('presence', 0),
            'percentage_certain_cells': derived_metrics_row.get('percentage_certain_cells', 0),
            'mean_certainty': derived_metrics_row.get('mean_certainty', 0)
        }
        derived_metrics_row['Composite_Score'] = calculate_composite_score(derived_composite_inputs)
        
        # Add derived metrics to results
        combined_metrics.append(derived_metrics_row)
        
    return pd.DataFrame(combined_metrics)


def calculate_additional_metrics(adata, ref_adata):
    """
    Calculate additional metrics (consistency, abundance, etc.) for all lineages
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with annotations
    ref_adata : AnnData, optional
        Reference AnnData object
        
    Returns:
    --------
    tuple
        (direct_metrics, derived_metrics) dictionaries with metrics for each lineage
    """
    # Ensure we have derived lineage column
    if 'final_anno_pred_lineage' not in adata.obs:
        adata = preprocess_data(adata)
    
    # Calculate overall consistency
    consistency_results = calculate_lineage_consistency(adata)
    
    # Get all lineages
    all_lineages = set(CELL_TYPE_TO_LINEAGE.values())
    
    # Calculate metrics for each lineage
    direct_metrics = {}
    derived_metrics = {}
    
    for lineage in all_lineages:
        log_message(f"Calculating additional metrics for lineage: {lineage}")
        
        # Calculate direct lineage metrics
        direct_metrics[lineage] = {
            'consistency': consistency_results['lineage_specific_consistency'].get(lineage, np.nan),
            'abundance': calculate_lineage_abundance(adata, lineage, 'final_lineage_pred'),
            'percentage_certain_cells': calculate_percentage_of_certain_cells(adata, lineage, 'final_lineage_pred'),
            'mean_certainty': calculate_mean_certainty(adata, lineage, 'final_lineage_pred'),
            'presence': calculate_cell_type_presence(adata, lineage, 'final_lineage_pred', ref_adata)
        }
        
        # Calculate derived lineage metrics
        derived_metrics[lineage] = {
            'consistency': consistency_results['lineage_specific_consistency'].get(lineage, np.nan),
            'abundance': calculate_lineage_abundance(adata, lineage, 'final_anno_pred_lineage'),
            'percentage_certain_cells': calculate_percentage_of_certain_cells(adata, lineage, 'final_anno_pred_lineage'),
            'mean_certainty': calculate_mean_certainty(adata, lineage, 'final_anno_pred_lineage'),
            'presence': calculate_cell_type_presence(adata, lineage, 'final_anno_pred_lineage', ref_adata)
        }
    
    return direct_metrics, derived_metrics

# ---------- MAIN PROCESSING FUNCTIONS ----------

def analyze_model(model_name, file_path, reference_file, marker_files, correlation_files, output_dir):
    """
    Analyze a single model using both pre-calculated metrics and on-the-fly calculation for other metrics
    
    Parameters:
    -----------
    model_name : str
        Name of the model
    file_path : str
        Path to the model's h5ad file
    reference_file : str
        Path to reference h5ad file
    marker_files : dict
        Dictionary with marker data for the model
    correlation_files : dict
        Dictionary with correlation data for the model
    output_dir : str
        Directory to save results
    
    Returns:
    --------
    dict
        Combined metrics for the model
    """
    log_message(f"Analyzing model: {model_name}")
    
    # Create output directories
    os.makedirs(os.path.join(output_dir, 'data'), exist_ok=True)
    
    # Load h5ad files for calculating additional metrics
    log_message(f"Loading {model_name} h5ad file for calculating additional metrics")
    try:
        adata = sc.read_h5ad(file_path)
        log_message(f"Loaded adata with {adata.n_obs} cells and {adata.n_vars} genes")
        
        # Preprocess to add derived lineage column
        adata = preprocess_data(adata)
        
        # Load reference data
        ref_adata = sc.read_h5ad(reference_file)
        log_message(f"Loaded reference data with {ref_adata.n_obs} cells and {ref_adata.n_vars} genes")
        
        # Calculate additional metrics
        log_message(f"Calculating additional metrics")
        direct_additional_metrics, derived_additional_metrics = calculate_additional_metrics(adata, ref_adata)
    except Exception as e:
        log_message(f"Error loading h5ad files or calculating additional metrics: {str(e)}")
        direct_additional_metrics, derived_additional_metrics = {}, {}
    
    # Get pre-calculated marker and correlation data
    direct_lineage_df = marker_files.get(model_name, {}).get('lineage', pd.DataFrame())
    celltype_df = marker_files.get(model_name, {}).get('celltype', pd.DataFrame())
    lineage_corr_df = correlation_files.get(model_name, {}).get('lineage', pd.DataFrame())
    celltype_corr_df = correlation_files.get(model_name, {}).get('celltype', pd.DataFrame())

    log_message(f"Cell type marker data for {model_name}: {celltype_df.head()}")
    
    # Group cell type metrics by lineage
    derived_lineage_df = group_celltype_metrics_by_lineage(celltype_df)
    celltype_corr_by_lineage_df = group_celltype_correlation_by_lineage(celltype_corr_df)

    
    # Save grouped metrics
    if not derived_lineage_df.empty:
        derived_lineage_df.to_csv(os.path.join(output_dir, 'data', f"{model_name}_derived_lineage_metrics.csv"), index=False)
        
    if not celltype_corr_by_lineage_df.empty:
        celltype_corr_by_lineage_df.to_csv(os.path.join(output_dir, 'data', f"{model_name}_derived_lineage_correlation.csv"), index=False)
        
    # Combine all metrics
    combined_df = combine_all_lineage_metrics(
        direct_lineage_df, 
        derived_lineage_df, 
        lineage_corr_df, 
        celltype_corr_by_lineage_df,
        direct_additional_metrics,
        derived_additional_metrics
    )
    
    # Save combined metrics
    if not combined_df.empty:
        combined_df.to_csv(os.path.join(output_dir, 'data', f"{model_name}_combined_lineage_metrics.csv"), index=False)
        
    # Create radar plots for each lineage
    if not combined_df.empty:
        for lineage in combined_df['Lineage'].unique():
            # Get direct metrics
            direct_row = combined_df[(combined_df['Lineage'] == lineage) & (combined_df['Type'] == 'Direct')]
            # Get derived metrics
            derived_row = combined_df[(combined_df['Lineage'] == lineage) & (combined_df['Type'] == 'Derived')]
            
            if not direct_row.empty and not derived_row.empty:
                create_radar_plot(
                    direct_row.iloc[0], 
                    derived_row.iloc[0], 
                    lineage, 
                    model_name, 
                    output_dir
                )
                
    # Return combined metrics
    return {
        'combined_df': combined_df,
        'derived_lineage_df': derived_lineage_df,
        'celltype_corr_by_lineage_df': celltype_corr_by_lineage_df
    }


def compare_models(all_metrics, output_dir):
    """
    Create comparison visualizations across models
    """
    log_message("Creating model comparison visualizations")
    
    # Create output directory
    os.makedirs(os.path.join(output_dir, 'comparisons'), exist_ok=True)
    
    # Combine all data into a single DataFrame
    all_data = []
    
    for model_name, metrics in all_metrics.items():
        if 'combined_df' in metrics and not metrics['combined_df'].empty:
            df = metrics['combined_df'].copy()
            df['Model'] = model_name
            all_data.append(df)
            
    if not all_data:
        log_message("No data available for comparison")
        return
        
    # Combine DataFrames
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Save to file
    combined_df.to_csv(os.path.join(output_dir, 'all_models_lineage_metrics.csv'), index=False)
    
    # Create comparison plots for each lineage
    for lineage in combined_df['Lineage'].unique():
        lineage_df = combined_df[combined_df['Lineage'] == lineage]
        
        # Skip lineages with insufficient data
        if len(lineage_df) < 2:
            continue
            
        # Create composite score comparison
        plt.figure(figsize=(12, 6))
        sns.barplot(data=lineage_df, x='Model', y='Composite_Score', hue='Type')
        plt.title(f"{lineage} - Composite Score Comparison")
        plt.ylabel('Composite Score')
        plt.ylim(0, 1)
        plt.xticks(rotation=45, ha='right')
        plt.legend(title="Prediction Type")
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'comparisons', f'{lineage}_composite_score_comparison.png'), dpi=300)
        plt.close()
        
        # Create key metrics comparison (using facets)
        key_metrics_to_plot = [
            {'direct': 'Precision', 'derived': 'Mean_Precision', 'title': 'Precision'},
            {'direct': 'F1_Score', 'derived': 'Mean_F1_Score', 'title': 'F1 Score'},
            {'direct': 'consistency', 'derived': 'consistency', 'title': 'Consistency'},
            {'direct': 'presence', 'derived': 'presence', 'title': 'Presence'}
        ]
        
        # Prepare data for plotting
        plot_data = []
        
        for _, row in lineage_df.iterrows():
            for metric_dict in key_metrics_to_plot:
                if row['Type'] == 'Direct':
                    metric_key = metric_dict['direct']
                    if metric_key in row and not pd.isna(row[metric_key]):
                        plot_data.append({
                            'Model': row['Model'],
                            'Type': 'Direct',
                            'Metric': metric_dict['title'],
                            'Value': row[metric_key]
                        })
                else:  # Derived
                    metric_key = metric_dict['derived']
                    if metric_key in row and not pd.isna(row[metric_key]):
                        plot_data.append({
                            'Model': row['Model'],
                            'Type': 'Derived',
                            'Metric': metric_dict['title'],
                            'Value': row[metric_key]
                        })
        
        # Create plot if we have data
        if plot_data:
            plot_df = pd.DataFrame(plot_data)
            
            # Create facet grid
            g = sns.FacetGrid(plot_df, col='Metric', height=4, aspect=1.2, sharey=True)
            g.map_dataframe(sns.barplot, x='Model', y='Value', hue='Type')
            g.add_legend(title='Prediction Type')
            g.set_axis_labels("Model", "Value")
            g.set_titles(col_template="{col_name}")
            
            # Rotate x-axis labels
            for ax in g.axes.flat:
                plt.sca(ax)
                plt.xticks(rotation=45, ha='right')
                plt.ylim(0, 1)
                plt.title(ax.get_title(), fontsize=12, pad=10)
                
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'comparisons', f'{lineage}_key_metrics_comparison.png'), dpi=300)
            plt.close()

# Add these functions to your existing code

def create_multi_model_radar_plots(all_metrics, output_dir, plot_type="direct"):
    """
    Create radar plots comparing all models for each lineage, either direct or derived lineages
    """
    log_message(f"Creating radar plots comparing all models for each lineage ({plot_type})")
    
    # Create output directory
    radar_dir = os.path.join(output_dir, f'radar_{plot_type}_lineages')
    os.makedirs(radar_dir, exist_ok=True)
    
    # Get all lineages and models
    all_lineages = set()
    all_model_names = set()
    
    # Collect all lineages and models
    for model_name, metrics in all_metrics.items():
        if 'combined_df' in metrics and not metrics['combined_df'].empty:
            all_model_names.add(model_name)
            df = metrics['combined_df']
            all_lineages.update(df['Lineage'].unique())
    
    # Manually add missing lineages if needed
    all_lineages = sorted(all_lineages | set(CELL_TYPE_TO_LINEAGE.values()))
    
    log_message(f"Found {len(all_lineages)} lineages and {len(all_model_names)} models")
    
    # Define metrics to include in radar plot based on your CSV column names
    metrics_to_plot = [
        'Precision', 'Recall', 'F1_Score', 'Jaccard', 'Pearson_r',
        'consistency', 'abundance', 'presence', 
        'percentage_certain_cells', 'mean_certainty'
    ]
    
    # Alternative column names that might be in the CSV
    alt_column_names = {
        'Precision': ['precision', 'Mean_Precision', 'mean_precision'],
        'Recall': ['recall', 'Mean_Recall', 'mean_recall'],
        'F1_Score': ['f1_score', 'Mean_F1_Score', 'mean_f1_score', 'F1 Score', 'f1'],
        'Jaccard': ['jaccard', 'Mean_Jaccard', 'mean_jaccard', 'jaccard_index', 'Jaccard Index'],
        'Pearson_r': ['Pearson_r', 'Mean_Pearson_r', 'pearson_r', 'correlation', 'Correlation'],
        'consistency': ['consistency', 'Consistency'],
        'abundance': ['abundance', 'Abundance'],
        'presence': ['presence', 'Presence'],
        'percentage_certain_cells': ['percentage_certain_cells', 'percentage_certain', 'percent_certain'],
        'mean_certainty': ['mean_certainty', 'certainty']
    }
    
    # For derived metrics, we need to adjust the column names
    if plot_type == "derived":
        primary_mapping = {
            'Precision': 'Mean_Precision',
            'Recall': 'Mean_Recall',
            'F1_Score': 'Mean_F1_Score',
            'Jaccard': 'Mean_Jaccard',
            'Pearson_r': 'Mean_Pearson_r'
        }
    else:
        primary_mapping = {m: m for m in metrics_to_plot}
    
    # Define labels for plotting
    metric_labels = [m.replace('_', ' ').capitalize() for m in metrics_to_plot]
    
    # Create a radar plot for each lineage
    for lineage in all_lineages:
        # Skip "Unknown" lineage
        if lineage == "Unknown":
            continue
            
        # Collect data for this lineage from all models
        lineage_data = []
        
        for model_name in all_model_names:
            if model_name not in all_metrics or 'combined_df' not in all_metrics[model_name]:
                continue
                
            combined_df = all_metrics[model_name]['combined_df']
            
            # Filter for this lineage and type
            type_value = "Derived" if plot_type == "derived" else "Direct"
            lineage_rows = combined_df[(combined_df['Lineage'] == lineage) & (combined_df['Type'] == type_value)]
            
            if lineage_rows.empty:
                continue
                
            # Get metrics for this model
            model_metrics = lineage_rows.iloc[0].to_dict()
            model_metrics['model_name'] = model_name
            
            lineage_data.append(model_metrics)
        
        # Skip if no data for this lineage
        if not lineage_data:
            log_message(f"No {plot_type} data for lineage {lineage}, skipping radar plot")
            continue
            
        # Create the radar plot
        plt.figure(figsize=(12, 10))
        ax = plt.subplot(111, polar=True)
        
        # Set angles for each metric
        angles = np.linspace(0, 2 * np.pi, len(metrics_to_plot), endpoint=False).tolist()
        angles += angles[:1]  # Close the loop
        
        # Set up the radar chart
        plt.xticks(angles[:-1], metric_labels, size=12)
        ax.set_rlabel_position(0)
        plt.yticks([0.2, 0.4, 0.6, 0.8], ["0.2", "0.4", "0.6", "0.8"], size=10)
        plt.ylim(0, 1)
        
        # Use different colors for each model
        colormap = plt.cm.get_cmap('tab10', len(lineage_data))
        
        # Plot each model
        for i, model_data in enumerate(lineage_data):
            # Extract values for each metric
            values = []
            for metric in metrics_to_plot:
                # First try the primary mapping
                mapped_metric = primary_mapping.get(metric, metric)
                value = 0
                
                # Try to get the value using the primary mapping
                if mapped_metric in model_data and not pd.isna(model_data[mapped_metric]):
                    value = model_data[mapped_metric]
                else:
                    # Try alternative column names
                    for alt_name in alt_column_names.get(metric, []):
                        if alt_name in model_data and not pd.isna(model_data[alt_name]):
                            value = model_data[alt_name]
                            break
                
                values.append(value)
            
            # Close the loop for plotting
            values += values[:1]
            
            # Calculate composite score
            if 'Composite_Score' in model_data and not pd.isna(model_data['Composite_Score']):
                composite_score = model_data['Composite_Score']
            elif 'composite_score' in model_data and not pd.isna(model_data['composite_score']):
                composite_score = model_data['composite_score']
            else:
                # Calculate it
                composite_inputs = {}
                for idx, metric in enumerate(metrics_to_plot):
                    composite_inputs[metric] = values[idx]
                composite_score = calculate_composite_score(composite_inputs)
            
            # Plot this model's data
            color = colormap(i)
            model_name = model_data['model_name']
            label = f"{model_name} (Score: {composite_score:.2f})"
            
            ax.plot(angles, values, 'o-', linewidth=2, label=label, color=color)
            ax.fill(angles, values, alpha=0.1, color=color)
        
        # Add legend and title
        plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
        plt.title(f"{lineage} Lineage - {plot_type.title()} Prediction - All Models Comparison")
        
        # Save the figure
        plt.tight_layout()
        plt.savefig(os.path.join(radar_dir, f"{lineage}_all_models_{plot_type}.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        log_message(f"Created radar plot for {lineage} ({plot_type}) with {len(lineage_data)} models")


def create_mean_score_radar_plots(all_metrics, output_dir):
    """
    Create radar plots showing mean score between direct and derived lineages for all models
    """
    log_message("Creating radar plots with mean scores between direct and derived lineages")
    
    # Create output directory
    radar_dir = os.path.join(output_dir, 'radar_mean_scores')
    os.makedirs(radar_dir, exist_ok=True)
    
    # Get all lineages and models
    all_lineages = set()
    all_model_names = list(all_metrics.keys())
    
    # Collect all lineages
    for model_name, metrics in all_metrics.items():
        if 'combined_df' in metrics and not metrics['combined_df'].empty:
            df = metrics['combined_df']
            all_lineages.update(df['Lineage'].unique())
    
    # Add all possible lineages
    all_lineages = sorted(all_lineages | set(CELL_TYPE_TO_LINEAGE.values()))
    
    log_message(f"Found {len(all_lineages)} lineages and {len(all_model_names)} models")
    
    # Define metrics to include in radar plot
    metrics_to_plot = [
        'Precision', 'Recall', 'F1_Score', 'Jaccard', 'Pearson_r',
        'consistency', 'abundance', 'presence', 
        'percentage_certain_cells', 'mean_certainty'
    ]
    
    # Alternative column names that might be in the CSV
    alt_column_names = {
        'Precision': ['precision', 'Mean_Precision', 'mean_precision'],
        'Recall': ['recall', 'Mean_Recall', 'mean_recall'],
        'F1_Score': ['f1_score', 'Mean_F1_Score', 'mean_f1_score', 'F1 Score', 'f1'],
        'Jaccard': ['jaccard', 'Mean_Jaccard', 'mean_jaccard', 'jaccard_index', 'Jaccard Index'],
        'Pearson_r': ['Pearson_r', 'Mean_Pearson_r', 'pearson_r', 'correlation', 'Correlation'],
        'consistency': ['consistency', 'Consistency'],
        'abundance': ['abundance', 'Abundance'],
        'presence': ['presence', 'Presence'],
        'percentage_certain_cells': ['percentage_certain_cells', 'percentage_certain', 'percent_certain'],
        'mean_certainty': ['mean_certainty', 'certainty']
    }
    
    # Mapping for derived metrics
    derived_mapping = {
        'Precision': 'Mean_Precision',
        'Recall': 'Mean_Recall',
        'F1_Score': 'Mean_F1_Score',
        'Jaccard': 'Mean_Jaccard',
        'Pearson_r': 'Mean_Pearson_r'
    }
    
    # Define labels for plotting
    metric_labels = [m.replace('_', ' ').capitalize() for m in metrics_to_plot]
    
    # Create a radar plot for each lineage
    for lineage in all_lineages:
        # Skip "Unknown" lineage
        if lineage == "Unknown":
            continue
            
        # Collect data for this lineage from all models
        lineage_data = []
        
        for model_name in all_model_names:
            if model_name not in all_metrics or 'combined_df' not in all_metrics[model_name]:
                continue
                
            combined_df = all_metrics[model_name]['combined_df']
            
            # Get direct and derived metrics for this lineage
            direct_rows = combined_df[(combined_df['Lineage'] == lineage) & (combined_df['Type'] == "Direct")]
            derived_rows = combined_df[(combined_df['Lineage'] == lineage) & (combined_df['Type'] == "Derived")]
            
            if direct_rows.empty or derived_rows.empty:
                log_message(f"Missing direct or derived data for {lineage} in model {model_name}, skipping")
                continue
                
            direct_metrics = direct_rows.iloc[0].to_dict()
            derived_metrics = derived_rows.iloc[0].to_dict()
            
            # Calculate mean metrics
            mean_metrics = {
                'model_name': model_name,
                'lineage': lineage
            }
            
            # For each metric, calculate mean between direct and derived
            for metric in metrics_to_plot:
                # Get direct value
                direct_val = 0
                if metric in direct_metrics and not pd.isna(direct_metrics[metric]):
                    direct_val = direct_metrics[metric]
                else:
                    # Try alternative column names
                    for alt_name in alt_column_names.get(metric, []):
                        if alt_name in direct_metrics and not pd.isna(direct_metrics[alt_name]):
                            direct_val = direct_metrics[alt_name]
                            break
                
                # Get derived value
                derived_val = 0
                derived_metric = derived_mapping.get(metric, metric)
                if derived_metric in derived_metrics and not pd.isna(derived_metrics[derived_metric]):
                    derived_val = derived_metrics[derived_metric]
                else:
                    # Try alternative column names
                    for alt_name in alt_column_names.get(metric, []):
                        if alt_name in derived_metrics and not pd.isna(derived_metrics[alt_name]):
                            derived_val = derived_metrics[alt_name]
                            break
                
                # Calculate mean
                mean_metrics[metric] = (direct_val + derived_val) / 2
            
            # Calculate composite score on mean values
            mean_metrics['composite_score'] = calculate_composite_score(mean_metrics)
            
            lineage_data.append(mean_metrics)
        
        # Skip if no data for this lineage
        if not lineage_data:
            log_message(f"No mean score data for lineage {lineage}, skipping radar plot")
            continue
            
        # Create the radar plot
        plt.figure(figsize=(12, 10))
        ax = plt.subplot(111, polar=True)
        
        # Set angles for each metric
        angles = np.linspace(0, 2 * np.pi, len(metrics_to_plot), endpoint=False).tolist()
        angles += angles[:1]  # Close the loop
        
        # Set up the radar chart
        plt.xticks(angles[:-1], metric_labels, size=12)
        ax.set_rlabel_position(0)
        plt.yticks([0.2, 0.4, 0.6, 0.8], ["0.2", "0.4", "0.6", "0.8"], size=10)
        plt.ylim(0, 1)
        
        # Use different colors for each model
        colormap = plt.cm.get_cmap('tab10', len(lineage_data))
        
        # Plot each model
        for i, model_data in enumerate(lineage_data):
            # Extract values for each metric
            values = []
            for metric in metrics_to_plot:
                if metric in model_data and not pd.isna(model_data[metric]):
                    values.append(model_data[metric])
                else:
                    values.append(0)
            
            # Close the loop for plotting
            values += values[:1]
            
            # Plot this model's data
            color = colormap(i)
            model_name = model_data['model_name']
            composite_score = model_data.get('composite_score', 0)
            label = f"{model_name} (Score: {composite_score:.2f})"
            
            ax.plot(angles, values, 'o-', linewidth=2, label=label, color=color)
            ax.fill(angles, values, alpha=0.1, color=color)
        
        # Add legend and title
        plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
        plt.title(f"{lineage} Lineage - Mean Between Direct and Derived - All Models")
        
        # Save the figure
        plt.tight_layout()
        plt.savefig(os.path.join(radar_dir, f"{lineage}_mean_score_all_models.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        log_message(f"Created mean score radar plot for {lineage} with {len(lineage_data)} models")

# ---------- MAIN FUNCTION ----------

def main():
    """Main function"""
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load pre-calculated metrics from CSV files
    log_message("Loading marker overlap files")
    marker_files = load_marker_overlap_files(MARKER_OVERLAP_DIR)
    log_message(f"Loaded marker data for {len(marker_files)} models")
    
    log_message("Loading expression correlation files")
    correlation_files = load_expression_correlation_files(EXPRESSION_CORR_DIR)
    log_message(f"Loaded correlation data for {len(correlation_files)} models")
    
    # Find h5ad files
    h5ad_files = {}
    file_pattern = "corrected_processed_*_scPoli_annotated.h5ad"
    files = glob(os.path.join(DATA_DIR, file_pattern))
    
    for file_path in files:
        model_name = extract_model_name(file_path)
        h5ad_files[model_name] = file_path
        log_message(f"Found h5ad file for model: {model_name}")
    
    # Get all model names
    all_models = set(marker_files.keys()) | set(correlation_files.keys()) | set(h5ad_files.keys())
    log_message(f"Found {len(all_models)} models to analyze")
    
    # Analyze each model
    all_metrics = {}
    
    for model_name in all_models:
        try:
            # Get h5ad file path if available
            file_path = h5ad_files.get(model_name)
            
            if file_path is None:
                log_message(f"Warning: No h5ad file found for {model_name}, skipping additional metrics calculation")
                continue
                
            model_metrics = analyze_model(
                model_name, 
                file_path, 
                REFERENCE_FILE, 
                marker_files, 
                correlation_files, 
                OUTPUT_DIR
            )
            all_metrics[model_name] = model_metrics
            log_message(f"Completed analysis for {model_name}")
        except Exception as e:
            log_message(f"Error analyzing {model_name}: {str(e)}")
            import traceback
            log_message(traceback.format_exc())
    

    # Compare metrics across models
    if len(all_metrics) > 1:
        try:
            compare_models(all_metrics, OUTPUT_DIR)
            log_message("Completed model comparison")
            
            # Generate multi-model radar plots for direct lineages
            create_multi_model_radar_plots(all_metrics, OUTPUT_DIR, plot_type="direct")
            
            # Generate multi-model radar plots for derived lineages
            create_multi_model_radar_plots(all_metrics, OUTPUT_DIR, plot_type="derived")
            
            # Generate mean score radar plots
            create_mean_score_radar_plots(all_metrics, OUTPUT_DIR)
            
        except Exception as e:
            log_message(f"Error in model comparison or radar plots: {str(e)}")
            import traceback
            log_message(traceback.format_exc())
    
    log_message(f"Analysis complete! Results saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()