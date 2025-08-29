#!/usr/bin/env python
# coding: utf-8
"""
Cell Type Evaluation Pipeline for Single-Cell Models

This script evaluates the performance of multiple single-cell annotation models (e.g., scPoli)
by comparing predicted cell types against a reference. It follows a simplified logic:

1.  Compute consistent cells (where reanno_pred == celltype_pred if available)
2.  For each reanno_pred cell type, compute consistency (%)
3.  Using consistent cells, compute mean certainty per cell type
4.  Using consistent cells, compute presence percentage (fraction of expected cell types detected)
5.  Using consistent cells, compute abundance per cell type
6.  Load precomputed cell-type expression correlations directly
7.  Load precomputed marker overlap metrics directly
8.  Compute a composite score combining all metrics

Outputs:
- Per-model cell-type-level metrics CSV
- Radar plots for each model-celltype
- Bar plots comparing composite scores across models
"""

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from collections import defaultdict

# ---------- PATHS ----------
DATA_DIR = '/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output'
MARKER_OVERLAP_DIR = "/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250806_embryo_model_benchmarking/marker_overlap_metrics"
EXPRESSION_CORR_DIR = "/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250806_embryo_model_benchmarking/expression_correlation_metrics"
REFERENCE_FILE = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/data/human_clustering_20250724_v3.h5ad'
OUTPUT_DIR = "./celltype_comparison_results"

# ---------- UTILITY FUNCTIONS ----------
def log_message(message, log_file='celltype_comparison.log'):
    """Log message to console and file with timestamp."""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

def get_consistent_cells_mask(adata):
    """
    Get a boolean mask for cells where reanno_pred_lineage == lineage_pred.
    Falls back to uncertainty threshold if consistency columns are missing.
    """
    if 'reanno_pred_lineage' in adata.obs and 'lineage_pred' in adata.obs:
        reanno_str = adata.obs['reanno_pred_lineage'].astype(str)
        lineage_str = adata.obs['lineage_pred'].astype(str)
        mask = (reanno_str == lineage_str).values
        log_message(f"Using consistent cells (reanno_pred_lineage == lineage_pred): {np.sum(mask)}/{len(mask)}")
        return mask

    # Fallback: Use certainty if available
    certainty_threshold = 0.2
    if 'is_celltype_certain' in adata.obs:
        mask = adata.obs['is_celltype_certain'].astype(bool).values
        log_message(f"Fallback: Using 'is_celltype_certain': {np.sum(mask)}/{len(mask)} certain cells")
        return mask
    elif 'reanno_uncert' in adata.obs:
        mask = adata.obs['reanno_uncert'].values < certainty_threshold
        log_message(f"Fallback: Using 'reanno_uncert' < {certainty_threshold}: {np.sum(mask)}/{len(mask)} certain")
        return mask
    else:
        mask = np.ones(adata.n_obs, dtype=bool)
        log_message("No consistency data; using all cells")
        return mask

# ---------- METRICS CALCULATION FUNCTIONS ----------
def calculate_celltype_consistency(adata):
    """
    Compute overall and per-celltype consistency between reanno_pred and celltype_pred.
    """
    required_cols = ['reanno_pred_lineage', 'lineage_pred']
    if not all(col in adata.obs for col in required_cols):
        log_message(f"Missing columns: {required_cols}")
        return {'error': 'Missing columns'}

    reanno_str = adata.obs['reanno_pred_lineage'].astype(str)
    lineage_str = adata.obs['lineage_pred'].astype(str)
    overall = (reanno_str == lineage_str).mean()

    # Calculate consistency per cell type based on lineage consistency
    celltype_consistency = {}
    if 'reanno_pred' in adata.obs:
        all_celltypes = adata.obs['reanno_pred'].dropna().unique()
        for celltype in all_celltypes:
            celltype_mask = adata.obs['reanno_pred'].astype(str) == str(celltype)
            if celltype_mask.sum() > 0:
                celltype_consistency[celltype] = (reanno_str[celltype_mask] == lineage_str[celltype_mask]).mean()
            else:
                celltype_consistency[celltype] = np.nan

    return {'overall_consistency': overall, 'celltype_specific_consistency': celltype_consistency}

def calculate_celltype_consistency_percentage(adata, celltype):
    """
    Step 2: Compute fraction of cells in a celltype that are consistent based on lineage consistency.
    """
    celltype_mask = adata.obs['reanno_pred'].astype(str) == str(celltype)
    if celltype_mask.sum() == 0:
        return 0.0

    if 'reanno_pred_lineage' in adata.obs and 'lineage_pred' in adata.obs:
        reanno_str = adata.obs['reanno_pred_lineage'].astype(str)
        lineage_str = adata.obs['lineage_pred'].astype(str)
        consistent_mask = (reanno_str == lineage_str).values
        return (celltype_mask & consistent_mask).sum() / celltype_mask.sum()
    return 1.0  # Assume consistent if no data

def calculate_celltype_mean_certainty_consistent_cells(adata, celltype, consistent_mask):
    """
    Step 3: Compute mean certainty for cells in a celltype among consistent cells.
    """
    celltype_mask = adata.obs['reanno_pred'].astype(str) == str(celltype)
    combined_mask = celltype_mask & consistent_mask
    if not combined_mask.any():
        return 0.0

    obs = adata.obs[combined_mask]
    if 'reanno_uncert' in obs:
        return (1 - obs['reanno_uncert']).mean()
    elif 'is_celltype_certain' in obs:
        return obs['is_celltype_certain'].astype(float).mean()
    else:
        return 1.0

def calculate_celltype_abundance_consistent_cells(adata, celltype, consistent_mask):
    """
    Step 5: Compute abundance of a celltype among consistent cells.
    """
    celltype_mask = adata.obs['reanno_pred'].astype(str) == str(celltype)
    celltype_count = (celltype_mask & consistent_mask).sum()
    total_consistent = consistent_mask.sum()
    return celltype_count / total_consistent if total_consistent > 0 else 0.0

def calculate_celltype_presence(adata, celltype, pred_type='reanno_pred', ref_adata=None):
    """
    Step 4: Compute if this celltype is present in the query data.
    For cell types, this is binary: 1 if present, 0 if not.
    """
    consistent_mask = get_consistent_cells_mask(adata)
    adata_cons = adata[consistent_mask]

    if pred_type not in adata_cons.obs:
        return 0.0
    
    observed_celltypes = set(adata_cons.obs[pred_type].dropna().astype(str).unique())
    return 1.0 if str(celltype) in observed_celltypes else 0.0

# ---------- LOADING PRE-CALCULATED METRICS ----------
def extract_model_name(filename):
    """Extract model name from filename"""
    # Remove file extension
    base_name = os.path.basename(filename)
    
    # Remove the suffix part
    if "_celltype_marker_overlap.csv" in base_name:
        model_name = base_name.replace("_celltype_marker_overlap.csv", "")
    elif "_cell_type_pearson_correlation.csv" in base_name:
        model_name = base_name.replace("_cell_type_pearson_correlation.csv", "")
    elif "_scPoli_annotated.h5ad" in base_name:
        model_name = base_name.replace(".h5ad", "")
    else:
        model_name = base_name
    
    return model_name

def load_marker_overlap_files(directory):
    """
    Load cell type marker overlap files from directory
    
    Returns:
    --------
    Dict with model names as keys, each containing 'celltype' DataFrames
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
            
    return marker_files

def load_expression_correlation_files(directory):
    """
    Load cell type expression correlation files from directory
    
    Returns:
    --------
    Dict with model names as keys, each containing 'celltype' DataFrames
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
            
    return correlation_files

def get_all_reference_celltypes(ref_adata):
    """
    Get all cell types from reference data.
    """
    if ref_adata is None or 'reanno' not in ref_adata.obs:
        log_message("Reference data missing or no 'reanno' column")
        return set()
    
    ref_celltypes = set(ref_adata.obs['reanno'].dropna().astype(str).unique())
    log_message(f"Found {len(ref_celltypes)} reference cell types")
    return ref_celltypes
def calculate_additional_metrics(adata, ref_adata):
    """Steps 1-5: Compute basic celltype metrics."""
    log_message("Starting celltype metrics calculation")
    consistent_mask = get_consistent_cells_mask(adata)
    log_message(f"Found {consistent_mask.sum()}/{len(consistent_mask)} consistent cells")

    if 'reanno_pred' not in adata.obs:
        log_message("Error: 'reanno_pred' not found")
        return {}

    # Get all cell types from reference
    all_celltypes = get_all_reference_celltypes(ref_adata)
    if not all_celltypes:
        # Fallback to query cell types if no reference
        all_celltypes = set(adata.obs['reanno_pred'].dropna().astype(str).unique())
    
    log_message(f"Processing {len(all_celltypes)} cell types")

    results = {}
    for celltype in all_celltypes:
        celltype_str = str(celltype)
        log_message(f"Processing cell type: {celltype_str}")

        # Step 2
        consistency = calculate_celltype_consistency_percentage(adata, celltype_str)
        # Step 3
        certainty = calculate_celltype_mean_certainty_consistent_cells(adata, celltype_str, consistent_mask)
        # Step 4
        presence = calculate_celltype_presence(adata, celltype_str, 'reanno_pred', ref_adata)
        # Step 5
        abundance = calculate_celltype_abundance_consistent_cells(adata, celltype_str, consistent_mask)

        results[celltype_str] = {
            'consistency': consistency,
            'mean_certainty': certainty,
            'presence': presence,
            'abundance': abundance
        }
    log_message("Completed basic celltype metrics calculation")
    return results

def calculate_composite_score(metrics_dict):
    """Step 8: Compute weighted average of all metrics."""
    weights = {
        'Precision': 0.05, 'Recall': 0.05, 'F1_Score': 0.05,
        'Jaccard': 0.05, 'Pearson_r': 0.2, 'consistency': 0.05,
        'abundance': 0.2, 'presence': 0.2, 'mean_certainty': 0.15
    }
    total_weight = 0.0
    weighted_sum = 0.0
    for metric, weight in weights.items():
        val = metrics_dict.get(metric, 0.0)
        if isinstance(val, (int, float)) and np.isfinite(val):
            weighted_sum += val * weight
            total_weight += weight
    return weighted_sum / total_weight if total_weight > 0 else 0.0

def combine_all_celltype_metrics(adata, ref_adata, marker_overlap_data, expression_correlation_data, model_name):
    """Combine all celltype metrics."""
    log_message(f"Combining all celltype metrics for model: {model_name}")
    
    # Get all cell types from reference data
    all_celltypes = get_all_reference_celltypes(ref_adata)
    if not all_celltypes:
        log_message("No reference cell types found, using query cell types")
        all_celltypes = set(adata.obs['reanno_pred'].dropna().astype(str).unique()) if 'reanno_pred' in adata.obs else set()
    
    log_message(f"Processing {len(all_celltypes)} cell types")
    
    # Calculate basic metrics (consistency, certainty, presence, abundance)
    additional_metrics = calculate_additional_metrics(adata, ref_adata)

    # Load correlation metrics
    celltype_correlation_metrics = {}
    if model_name in expression_correlation_data and 'celltype' in expression_correlation_data[model_name]:
        celltype_corr_df = expression_correlation_data[model_name]['celltype']
        for _, row in celltype_corr_df.iterrows():
            celltype_str = str(row['Label'])
            celltype_correlation_metrics[celltype_str] = {
                'Pearson_r': row.get('Pearson_r', 0),
                'Pearson_p': row.get('Pearson_p', 1),
                'Query_Cells': row.get('Query_Cells', 0)
            }

    # Load marker metrics
    celltype_marker_metrics = {}
    if model_name in marker_overlap_data and 'celltype' in marker_overlap_data[model_name]:
        celltype_marker_df = marker_overlap_data[model_name]['celltype']
        for _, row in celltype_marker_df.iterrows():
            celltype_str = str(row['Group'])
            celltype_marker_metrics[celltype_str] = {
                'Precision': row.get('Precision', 0),
                'Recall': row.get('Recall', 0),
                'F1_Score': row.get('F1_Score', 0),
                'Jaccard': row.get('Jaccard', 0),
                'Overlap': row.get('Overlap', 0),
                'Query_Markers': row.get('Query_Markers', 0),
                'Ref_Markers': row.get('Ref_Markers', 0)
            }

    combined_metrics = []

    # Process all reference cell types
    for celltype in all_celltypes:
        celltype_str = str(celltype)
        log_message(f"Processing cell type for model {model_name}: {celltype_str}")

        # Start with a base row of zeros
        metrics_row = {
            'Model': model_name,
            'CellType': celltype_str,
            'consistency': 0.0,
            'mean_certainty': 0.0,
            'presence': 0.0,
            'abundance': 0.0,
            'Pearson_r': 0.0,
            'Precision': 0.0,
            'Recall': 0.0,
            'F1_Score': 0.0,
            'Jaccard': 0.0
        }

        # Update with basic metrics if they exist
        if celltype_str in additional_metrics:
            metrics_row.update({
                'consistency': additional_metrics[celltype_str].get('consistency', 0.0),
                'mean_certainty': additional_metrics[celltype_str].get('mean_certainty', 0.0),
                'presence': additional_metrics[celltype_str].get('presence', 0.0),
                'abundance': additional_metrics[celltype_str].get('abundance', 0.0)
            })

        # Update with correlation metrics if they exist
        if celltype_str in celltype_correlation_metrics:
            corr_data = celltype_correlation_metrics[celltype_str]
            metrics_row['Pearson_r'] = corr_data.get('Pearson_r', 0.0)

        # Update with marker metrics if they exist
        if celltype_str in celltype_marker_metrics:
            marker_data = celltype_marker_metrics[celltype_str]
            metrics_row.update({
                'Precision': marker_data.get('Precision', 0.0),
                'Recall': marker_data.get('Recall', 0.0),
                'F1_Score': marker_data.get('F1_Score', 0.0),
                'Jaccard': marker_data.get('Jaccard', 0.0)
            })

        # Calculate composite score for this row
        metrics_row['Composite_Score'] = calculate_composite_score(metrics_row)
        combined_metrics.append(metrics_row)

    return pd.DataFrame(combined_metrics)

# ---------- PLOTTING ----------
def create_simplified_radar_plot(metrics_row, celltype, model_name, output_dir):
    """Create a radar plot for one model-celltype pair."""
    metrics = ['Precision', 'Recall', 'F1_Score', 'Jaccard', 'Pearson_r',
               'consistency', 'abundance', 'presence', 'mean_certainty']
    values = [metrics_row.get(m, 0) * 100 for m in metrics]
    labels = [m.replace('_', ' ').title() for m in metrics]

    N = len(metrics)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]
    values += values[:1]

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    ax.plot(angles, values, 'o-', linewidth=2)
    ax.fill(angles, values, alpha=0.25)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 100)
    ax.grid(True)
    
    # Clean celltype name for filename
    clean_celltype = "".join(c if c.isalnum() or c in ".-_" else "_" for c in celltype)
    plt.title(f"{model_name} - {celltype}\nComposite Score: {metrics_row.get('Composite_Score', 0):.3f}",
              size=14, pad=20)
    plt.tight_layout()
    out_path = os.path.join(output_dir, 'radar_plots', f"{model_name}_{clean_celltype}_radar.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_radar_plots_for_model(df, model_name, output_dir):
    """Create radar plots for all celltypes in a model."""
    for _, row in df.iterrows():
        try:
            create_simplified_radar_plot(row, row['CellType'], model_name, output_dir)
        except Exception as e:
            log_message(f"Failed to create radar plot for {model_name}-{row['CellType']}: {e}")

def compare_models(all_metrics, output_dir):
    """Create comparison bar plots across models."""
    data = [df.assign(Model=name) for name, df in all_metrics.items() if isinstance(df, pd.DataFrame) and not df.empty]
    if not data:
        return
    combined = pd.concat(data, ignore_index=True)
    combined.to_csv(os.path.join(output_dir, 'all_models_celltype_metrics.csv'), index=False)

    # Create comparison plots for top cell types by abundance
    os.makedirs(os.path.join(output_dir, 'comparisons'), exist_ok=True)
    
    # Get top cell types by average abundance across all models
    top_celltypes = combined.groupby('CellType')['abundance'].mean().sort_values(ascending=False).head(20).index
    
    for celltype in top_celltypes:
        df = combined[combined['CellType'] == celltype]
        if len(df) < 2:
            continue
        plt.figure(figsize=(10, 6))
        sns.barplot(data=df, x='Model', y='Composite_Score')
        clean_celltype = "".join(c if c.isalnum() or c in ".-_" else "_" for c in celltype)
        plt.title(f"{celltype} - Composite Score Comparison")
        plt.ylim(0, 1)
        plt.xticks(rotation=45)
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'comparisons', f"{clean_celltype}_composite.png"), dpi=300)
        plt.close()

# ---------- MAIN ----------
def analyze_model(model_name, file_path, ref_file, marker_files, corr_files, output_dir):
    """Analyze one model."""
    log_message(f"\n{'='*60}\nAnalyzing model: {model_name}\n{'='*60}")
    os.makedirs(os.path.join(output_dir, 'data'), exist_ok=True)

    try:
        adata = sc.read_h5ad(file_path)
        ref_adata = sc.read_h5ad(ref_file) if ref_file else None
        log_message(f"Loaded query: {adata.n_obs} cells")
        if ref_adata:
            log_message(f"Loaded reference: {ref_adata.n_obs} cells")

        marker_data = {model_name: marker_files.get(model_name, {})}
        corr_data = {model_name: corr_files.get(model_name, {})}

        df = combine_all_celltype_metrics(adata, ref_adata, marker_data, corr_data, model_name)
        if df.empty:
            return df

        out_file = os.path.join(output_dir, 'data', f"{model_name}_celltype_metrics.csv")
        df.to_csv(out_file, index=False)
        log_message(f"Saved to {out_file}")

        create_radar_plots_for_model(df, model_name, output_dir)
        return df
    except Exception as e:
        log_message(f"Error analyzing {model_name}: {e}")
        return pd.DataFrame()

def main():
    """Main entry point."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    log_message("Loading marker overlap files")
    marker_files = load_marker_overlap_files(MARKER_OVERLAP_DIR)
    log_message(f"Loaded marker data for {len(marker_files)} models")

    log_message("Loading expression correlation files")
    corr_files = load_expression_correlation_files(EXPRESSION_CORR_DIR)
    log_message(f"Loaded correlation data for {len(corr_files)} models")

    h5ad_files = {}
    pattern = "corrected_processed_*_scPoli_annotated.h5ad"
    for path in glob(os.path.join(DATA_DIR, pattern)):
        name = extract_model_name(path)
        h5ad_files[name] = path
        log_message(f"Found h5ad file for {name}")

    all_models = set(h5ad_files.keys()) | set(marker_files.keys()) | set(corr_files.keys())
    log_message(f"Found {len(all_models)} models to analyze")

    all_metrics = {}
    for model in all_models:
        if model not in h5ad_files:
            log_message(f"Skipping {model}: no h5ad file")
            continue
        df = analyze_model(model, h5ad_files[model], REFERENCE_FILE, marker_files, corr_files, OUTPUT_DIR)
        all_metrics[model] = df

    if len(all_metrics) > 1:
        compare_models(all_metrics, OUTPUT_DIR)
        log_message("Model comparison completed.")

    log_message(f"Analysis finished. Results in {OUTPUT_DIR}")

if __name__ == "__main__":
    main()