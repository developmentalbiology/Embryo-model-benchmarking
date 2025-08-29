#!/usr/bin/env python
# coding: utf-8
"""
Lineage Evaluation Pipeline for Single-Cell Models

This script evaluates the performance of multiple single-cell annotation models (e.g., scPoli)
by comparing predicted lineages against a reference. It follows an 8-step logic:

1.  Compute consistent cells (where reanno_pred_lineage == lineage_pred)
2.  For each reanno_pred_lineage, compute consistency (%)
3.  Using consistent cells, compute mean certainty per lineage
4.  Using consistent cells, compute presence percentage (fraction of expected cell types detected)
5.  Using consistent cells, compute abundance per lineage
6.  Using precomputed cell-type expression correlations, compute mean correlation per lineage
7.  Using precomputed marker overlap metrics, compute mean precision/recall/F1/Jaccard per lineage
8.  Compute a composite score combining all metrics

Outputs:
- Per-model lineage-level metrics CSV
- Radar plots for each model-lineage
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
OUTPUT_DIR = "./lineage_comparison_results"

# ---------- UTILITY FUNCTIONS ----------
def log_message(message, log_file='lineage_comparison.log'):
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
    if 'is_lineage_certain' in adata.obs:
        mask = adata.obs['is_lineage_certain'].astype(bool).values
        log_message(f"Fallback: Using 'is_lineage_certain': {np.sum(mask)}/{len(mask)} certain cells")
        return mask
    elif 'lineage_uncert' in adata.obs:
        mask = adata.obs['lineage_uncert'].values < certainty_threshold
        log_message(f"Fallback: Using 'lineage_uncert' < {certainty_threshold}: {np.sum(mask)}/{len(mask)} certain")
        return mask
    else:
        mask = np.ones(adata.n_obs, dtype=bool)
        log_message("No consistency data; using all cells")
        return mask

# ---------- METRICS CALCULATION FUNCTIONS ----------
def calculate_lineage_consistency(adata):
    """
    Compute overall and per-lineage consistency between reanno_pred_lineage and lineage_pred.
    """
    required_cols = ['lineage_pred', 'reanno_pred_lineage']
    if not all(col in adata.obs for col in required_cols):
        log_message(f"Missing columns: {required_cols}")
        return {'error': 'Missing columns'}

    reanno_str = adata.obs['reanno_pred_lineage'].astype(str)
    lineage_str = adata.obs['lineage_pred'].astype(str)
    overall = (reanno_str == lineage_str).mean()

    lineage_consistency = {}
    all_lineages = set(reanno_str.unique()) | set(lineage_str.unique())
    for lineage in all_lineages:
        mask = (reanno_str == lineage) | (lineage_str == lineage)
        if mask.sum() > 0:
            lineage_consistency[lineage] = (reanno_str[mask] == lineage_str[mask]).mean()
        else:
            lineage_consistency[lineage] = np.nan

    return {'overall_consistency': overall, 'lineage_specific_consistency': lineage_consistency}

def calculate_lineage_abundance(adata, lineage, pred_type='lineage_pred'):
    """
    Compute proportion of consistent cells assigned to a lineage.
    """
    consistent_mask = get_consistent_cells_mask(adata)
    adata_cons = adata[consistent_mask]
    total = adata_cons.n_obs
    if total == 0:
        return 0.0

    col = 'lineage_pred' if pred_type == 'lineage_pred' else 'reanno_pred_lineage'
    mask = adata_cons.obs[col].astype(str) == lineage
    return np.sum(mask) / total

def calculate_lineage_consistency_percentage(adata, lineage):
    """
    Step 2: Compute fraction of cells in a lineage that are consistent.
    Uses reanno_pred_lineage as the primary grouping.
    """
    lineage_mask = adata.obs['reanno_pred_lineage'].astype(str) == str(lineage)
    if lineage_mask.sum() == 0:
        return 0.0

    if 'reanno_pred_lineage' in adata.obs and 'lineage_pred' in adata.obs:
        reanno_str = adata.obs['reanno_pred_lineage'].astype(str)
        lineage_str = adata.obs['lineage_pred'].astype(str)
        consistent_mask = (reanno_str == lineage_str).values
        return (lineage_mask & consistent_mask).sum() / lineage_mask.sum()
    return 1.0  # Assume consistent if no data

def calculate_lineage_mean_certainty_consistent_cells(adata, lineage, consistent_mask):
    """
    Step 3: Compute mean certainty for cells in a lineage among consistent cells.
    """
    lineage_mask = adata.obs['reanno_pred_lineage'].astype(str) == str(lineage)
    combined_mask = lineage_mask & consistent_mask
    if not combined_mask.any():
        return 0.0

    obs = adata.obs[combined_mask]
    if 'reanno_uncert' in obs:
        return (1 - obs['reanno_uncert']).mean()
    elif 'is_celltype_certain' in obs:
        return obs['is_celltype_certain'].astype(float).mean()
    else:
        return 1.0

def calculate_lineage_abundance_consistent_cells(adata, lineage, consistent_mask):
    """
    Step 5: Compute abundance of a lineage among consistent cells.
    """
    lineage_mask = adata.obs['reanno_pred_lineage'].astype(str) == str(lineage)
    lineage_count = (lineage_mask & consistent_mask).sum()
    total_consistent = consistent_mask.sum()
    return lineage_count / total_consistent if total_consistent > 0 else 0.0

def calculate_cell_type_presence(adata, lineage, pred_type='reanno_pred_lineage', ref_adata=None):
    """
    Step 4: Compute fraction of expected cell types (from ref) detected in query.
    """
    consistent_mask = get_consistent_cells_mask(adata)
    adata_cons = adata[consistent_mask]

    # Get expected cell types from reference
    if ref_adata is None or 'lineage' not in ref_adata.obs or 'reanno' not in ref_adata.obs:
        log_message("Reference missing lineage/reanno; skipping presence")
        return 0.0

    ref_lineage_mask = ref_adata.obs['lineage'].astype(str) == str(lineage)
    expected_types = set(ref_adata.obs.loc[ref_lineage_mask, 'reanno'].dropna().unique())
    if not expected_types:
        return 0.0

    # Get observed cell types in query
    lineage_mask = adata_cons.obs[pred_type].astype(str) == str(lineage)
    if 'reanno_pred' not in adata_cons.obs:
        return 0.0
    observed_types = set(adata_cons.obs.loc[lineage_mask, 'reanno_pred'].dropna().unique())

    present = [ct for ct in expected_types if ct in observed_types]
    score = len(present) / len(expected_types)

    log_message(f"Lineage '{lineage}' presence: {score:.3f} ({len(present)}/{len(expected_types)})")
    return score

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
    """Group cell-type metrics by 'Lineage' column."""
    if celltype_df.empty or 'Lineage' not in celltype_df:
        log_message("Cannot group celltype metrics: missing 'Lineage' column")
        return pd.DataFrame()
    return celltype_df.groupby('Lineage').mean(numeric_only=True).reset_index()

def group_celltype_correlation_by_lineage(celltype_corr_df):
    """Group cell-type correlation metrics by 'Lineage' column."""
    if celltype_corr_df.empty or 'Lineage' not in celltype_corr_df:
        log_message("Cannot group correlation data: missing 'Lineage' column")
        return pd.DataFrame()
    return celltype_corr_df.groupby('Lineage')['Pearson_r'].mean().reset_index()

# ---------- 8-STEP LINEAGE METRICS CALCULATION ----------
def calculate_additional_metrics(adata, ref_adata):
    """Steps 1-5: Compute basic lineage metrics."""
    log_message("Starting 8-step lineage metrics calculation")
    consistent_mask = get_consistent_cells_mask(adata)
    log_message(f"Found {consistent_mask.sum()}/{len(consistent_mask)} consistent cells")

    if 'reanno_pred_lineage' not in adata.obs:
        log_message("Error: 'reanno_pred_lineage' not found")
        return {}

    lineages = adata.obs['reanno_pred_lineage'].dropna().unique()
    log_message(f"Processing {len(lineages)} lineages: {list(lineages)}")

    results = {}
    for lineage in lineages:
        log_message(f"Processing lineage: {lineage}")

        # Step 2
        consistency = calculate_lineage_consistency_percentage(adata, lineage)
        # Step 3
        certainty = calculate_lineage_mean_certainty_consistent_cells(adata, lineage, consistent_mask)
        # Step 4
        presence = calculate_cell_type_presence(adata, lineage, 'reanno_pred_lineage', ref_adata)
        # Step 5
        abundance = calculate_lineage_abundance_consistent_cells(adata, lineage, consistent_mask)

        results[str(lineage)] = {
            'consistency': consistency,
            'mean_certainty': certainty,
            'presence': presence,
            'abundance': abundance
        }
    log_message("Completed steps 1-5")
    return results

def aggregate_celltype_metrics_to_lineage_by_reanno_mapping(adata, consistent_mask, celltype_metrics_dict, metric_name):
    """
    Steps 6-7: Aggregate cell type metrics to lineage level using reanno and reanno_pred_lineage mapping.
    """
    lineage_aggregated_metrics = {}
    all_lineages = adata.obs['reanno_pred_lineage'].unique()
    for lineage in all_lineages:
        log_message(f"Aggregating {metric_name} for lineage: {lineage}")
        lineage_mask = adata.obs['reanno_pred_lineage'].astype(str) == str(lineage)
        lineage_consistent_mask = lineage_mask & consistent_mask
        if not np.any(lineage_consistent_mask):
            lineage_aggregated_metrics[str(lineage)] = {}
            continue

        lineage_consistent_obs = adata.obs[lineage_consistent_mask]
        if 'reanno_pred' not in lineage_consistent_obs.columns:
            lineage_aggregated_metrics[str(lineage)] = {}
            continue

        lineage_celltypes = lineage_consistent_obs['reanno_pred'].unique()
        aggregated_metrics = {}
        valid_celltypes = []
        for celltype in lineage_celltypes:
            celltype_str = str(celltype)
            if celltype_str in celltype_metrics_dict:
                valid_celltypes.append(celltype_str)

        if not valid_celltypes:
            lineage_aggregated_metrics[str(lineage)] = {}
            continue

        celltype_data = [celltype_metrics_dict[ct] for ct in valid_celltypes]
        if celltype_data:
            metric_columns = celltype_data[0].keys()
            for metric_col in metric_columns:
                values = []
                for ct_data in celltype_data:
                    if metric_col in ct_data and not pd.isna(ct_data[metric_col]):
                        values.append(ct_data[metric_col])
                if values:
                    aggregated_metrics[metric_col] = np.mean(values)
                else:
                    aggregated_metrics[metric_col] = 0.0

        lineage_aggregated_metrics[str(lineage)] = aggregated_metrics
        log_message(f"Lineage '{lineage}' aggregated from {len(valid_celltypes)} cell types")
    return lineage_aggregated_metrics

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

# define all lineages
ALL_LINEAGES = [
    'Primitive.streak', 'Amniotic_ecto', 'epi', 'meso_Exe.meso', 
    'neural_ecto', 'Notochord', 'hemogenic', 'TE_TrB', 'PGC',
    'ExE_endo', 'NMP', 'Endoderm'
]

def combine_all_lineage_metrics(adata, ref_adata, marker_overlap_data, expression_correlation_data, model_name):
    """Combine all lineage metrics following the 8-step logic."""
    log_message(f"Combining all lineage metrics for model: {model_name}")
    
    # Calculate basic metrics (consistency, certainty, presence, abundance)
    additional_metrics = calculate_additional_metrics(adata, ref_adata)
    if not additional_metrics:
        log_message("No additional metrics calculated")
        return pd.DataFrame()

    consistent_mask = get_consistent_cells_mask(adata)

    # Load and aggregate correlation metrics
    lineage_correlation_metrics = {}
    if model_name in expression_correlation_data and 'celltype' in expression_correlation_data[model_name]:
        celltype_corr_df = expression_correlation_data[model_name]['celltype']
        celltype_corr_dict = {
            str(row['Label']): {
                'Pearson_r': row.get('Pearson_r', 0),
                'Pearson_p': row.get('Pearson_p', 1),
                'Query_Cells': row.get('Query_Cells', 0)
            }
            for _, row in celltype_corr_df.iterrows()
        }
        lineage_correlation_metrics = aggregate_celltype_metrics_to_lineage_by_reanno_mapping(
            adata, consistent_mask, celltype_corr_dict, "expression correlation"
        )

    # Load and aggregate marker metrics
    lineage_marker_metrics = {}
    if model_name in marker_overlap_data and 'celltype' in marker_overlap_data[model_name]:
        celltype_marker_df = marker_overlap_data[model_name]['celltype']
        celltype_marker_dict = {
            str(row['Group']): {
                'Precision': row.get('Precision', 0),
                'Recall': row.get('Recall', 0),
                'F1_Score': row.get('F1_Score', 0),
                'Jaccard': row.get('Jaccard', 0),
                'Overlap': row.get('Overlap', 0),
                'Query_Markers': row.get('Query_Markers', 0),
                'Ref_Markers': row.get('Ref_Markers', 0)
            }
            for _, row in celltype_marker_df.iterrows()
        }
        lineage_marker_metrics = aggregate_celltype_metrics_to_lineage_by_reanno_mapping(
            adata, consistent_mask, celltype_marker_dict, "marker overlap"
        )

    combined_metrics = []

    for lineage in ALL_LINEAGES:
        log_message(f"Processing lineage for model {model_name}: {lineage}")

        # Start with a base row of zeros
        metrics_row = {
            'Model': model_name,
            'Lineage': lineage,
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

        # Only update values if metrics exist for this lineage
        if lineage in additional_metrics:
            metrics_row.update({
                'consistency': additional_metrics[lineage].get('consistency', 0.0),
                'mean_certainty': additional_metrics[lineage].get('mean_certainty', 0.0),
                'presence': additional_metrics[lineage].get('presence', 0.0),
                'abundance': additional_metrics[lineage].get('abundance', 0.0)
            })

        if lineage in lineage_correlation_metrics:
            corr_data = lineage_correlation_metrics[lineage]
            metrics_row['Pearson_r'] = corr_data.get('Pearson_r', 0.0)

        if lineage in lineage_marker_metrics:
            marker_data = lineage_marker_metrics[lineage]
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
def create_simplified_radar_plot(metrics_row, lineage, model_name, output_dir):
    """Create a radar plot for one model-lineage pair."""
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
    plt.title(f"{model_name} - {lineage}\nComposite Score: {metrics_row.get('Composite_Score', 0):.3f}",
              size=14, pad=20)
    plt.tight_layout()
    out_path = os.path.join(output_dir, 'radar_plots', f"{model_name}_{lineage}_radar.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_radar_plots_for_model(df, model_name, output_dir):
    """Create radar plots for all lineages in a model."""
    for _, row in df.iterrows():
        try:
            create_simplified_radar_plot(row, row['Lineage'], model_name, output_dir)
        except Exception as e:
            log_message(f"Failed to create radar plot for {model_name}-{row['Lineage']}: {e}")

def compare_models(all_metrics, output_dir):
    """Create comparison bar plots across models."""
    data = [df.assign(Model=name) for name, df in all_metrics.items() if isinstance(df, pd.DataFrame) and not df.empty]
    if not data:
        return
    combined = pd.concat(data, ignore_index=True)
    combined.to_csv(os.path.join(output_dir, 'all_models_lineage_metrics.csv'), index=False)

    for lineage in combined['Lineage'].unique():
        df = combined[combined['Lineage'] == lineage]
        if len(df) < 2:
            continue
        plt.figure(figsize=(10, 6))
        sns.barplot(data=df, x='Model', y='Composite_Score')
        plt.title(f"{lineage} - Composite Score Comparison")
        plt.ylim(0, 1)
        plt.xticks(rotation=45)
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'comparisons', f"{lineage}_composite.png"), dpi=300)
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

        df = combine_all_lineage_metrics(adata, ref_adata, marker_data, corr_data, model_name)
        if df.empty:
            return df

        out_file = os.path.join(output_dir, 'data', f"{model_name}_lineage_metrics.csv")
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
