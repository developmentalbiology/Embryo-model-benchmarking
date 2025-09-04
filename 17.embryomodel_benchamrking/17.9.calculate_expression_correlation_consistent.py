#!/usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse as sp_issparse
from scipy.stats import pearsonr
import scipy.sparse as sp

# ---------- FUNCTIONS ----------

def log_message(message, log_file='expression_correlation.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")


def calculate_matching_label_correlation(query_data, reference_data, 
                                       query_label_key, ref_label_key,
                                       consistency_key=None):
    """
    Calculate gene expression Pearson correlation between matching labels in query and reference
    
    Parameters
    ----------
    query_data : AnnData
        Query dataset with predictions and gene expression
    reference_data : AnnData
        Reference atlas with annotations and gene expression
    query_label_key : str
        Column in query_data.obs containing predicted labels (e.g., 'reanno_pred')
    ref_label_key : str
        Column in reference_data.obs containing reference labels (e.g., 'reanno')
    consistency_key : str, optional
        Column in query_data.obs indicating which cells are consistent
        
    Returns
    -------
    dict
        Dictionary with correlation metrics and detailed dataframe
    """
    log_message(f"Calculating Pearson correlation between matching labels using {query_label_key} and {ref_label_key}...")
    
    # Check if required columns exist
    if query_label_key not in query_data.obs.columns:
        log_message(f"Error: {query_label_key} not found in query data. Available columns: {query_data.obs.columns}")
        return {"error": f"{query_label_key} not found"}
    
    if ref_label_key not in reference_data.obs.columns:
        log_message(f"Error: {ref_label_key} not found in reference data. Available columns: {reference_data.obs.columns}")
        return {"error": f"{ref_label_key} not found"}
    
    # Get consistency mask - Only use consistent cells if available
    if consistency_key is not None and consistency_key in query_data.obs:
        consistent_mask = query_data.obs[consistency_key].values.astype(bool)
        log_message(f"Using {np.sum(consistent_mask)}/{len(consistent_mask)} consistent cells")
        if np.sum(consistent_mask) == 0:
            log_message(f"Error: No cells meet consistency criteria with key '{consistency_key}'")
            return {"error": "No consistent cells found"}
    else:
        # Check for consistency based on reanno_pred_lineage == lineage_pred
        if 'reanno_pred_lineage' in query_data.obs and 'lineage_pred' in query_data.obs:
            # Convert categorical columns to string to avoid category comparison issues
            reanno_pred_lineage_str = query_data.obs['reanno_pred_lineage'].astype(str)
            lineage_pred_str = query_data.obs['lineage_pred'].astype(str)
            consistent_mask = (reanno_pred_lineage_str == lineage_pred_str).values
            log_message(f"Using consistency based on reanno_pred_lineage==lineage_pred: {np.sum(consistent_mask)}/{len(consistent_mask)} cells")
        else:
            # Fallback to old certainty logic if consistency columns not available
            if 'is_celltype_certain' in query_data.obs and query_label_key == 'reanno_pred':
                consistent_mask = query_data.obs['is_celltype_certain'].values.astype(bool)
                log_message(f"Fallback: Using 'is_celltype_certain': {np.sum(consistent_mask)}/{len(consistent_mask)} cells")
            elif 'is_lineage_certain' in query_data.obs and query_label_key == 'lineage_pred':
                consistent_mask = query_data.obs['is_lineage_certain'].values.astype(bool)
                log_message(f"Fallback: Using 'is_lineage_certain': {np.sum(consistent_mask)}/{len(consistent_mask)} cells")
            else:
                consistent_mask = np.ones(query_data.n_obs, dtype=bool)
                log_message(f"No consistency filter found or provided, using all {np.sum(consistent_mask)} cells")
    
    # Get filtered query data - only consistent cells
    query_filtered = query_data[consistent_mask].copy()
    
    # Get unique labels from the query dataset (focus on labels in embryo models)
    query_labels = query_filtered.obs[query_label_key].unique().tolist()
    ref_labels = reference_data.obs[ref_label_key].unique().tolist()
    
    # Find common labels (only compare these)
    common_labels = list(set(query_labels) & set(ref_labels))
    log_message(f"Found {len(common_labels)} labels common to both query and reference")
    
    # Skip if no common labels
    if not common_labels:
        log_message("No common labels found between query and reference")
        return {"error": "No common labels found"}
    
    # Find common genes
    common_genes = np.intersect1d(query_data.var_names, reference_data.var_names)
    log_message(f"Found {len(common_genes)} common genes between query and reference")
    
    if len(common_genes) < 100:
        log_message("Warning: Too few common genes, correlation may be unreliable")
        return {"error": "Too few common genes"}
    
    # Calculate average gene expression for each label in both datasets
    correlations = {}
    
    # Create dataframe to hold detailed correlation values
    correlation_df_data = []
    
    for label in common_labels:
        # Extract cells with this label from query (already filtered for consistent cells)
        query_mask = query_filtered.obs[query_label_key] == label
    
        # Extract cells with this label from reference
        ref_mask = reference_data.obs[ref_label_key] == label
    
        # Only skip if there are no cells (not if there are fewer than 5)
        if np.sum(query_mask) == 0:
            log_message(f"Skipping {label}: No consistent cells in query")
            continue
    
        if np.sum(ref_mask) == 0:
            log_message(f"Skipping {label}: No cells in reference")
            continue
    
        # Log the count but proceed with calculation
        log_message(f"Processing label: {label} - {np.sum(query_mask)} consistent query cells, {np.sum(ref_mask)} reference cells")
        
        # Get average expression for common genes in query
        query_expr = query_filtered[query_mask, common_genes].X
        if sp.issparse(query_expr):
            query_expr = query_expr.toarray()
        query_avg = np.mean(query_expr, axis=0)
        
        # Get average expression for common genes in reference
        ref_expr = reference_data[ref_mask, common_genes].X
        if sp.issparse(ref_expr):
            ref_expr = ref_expr.toarray()
        ref_avg = np.mean(ref_expr, axis=0)
        
        # Calculate Pearson correlation
        try:
            # Pearson correlation
            pearson_r, pearson_p = pearsonr(query_avg, ref_avg)
            
            correlations[label] = {
                "pearson_r": float(pearson_r),
                "pearson_p": float(pearson_p),
                "query_cells": int(np.sum(query_mask)),
                "ref_cells": int(np.sum(ref_mask))
            }
            
            # Add to correlation dataframe
            correlation_df_data.append({
                "Label": label,
                "Pearson_r": float(pearson_r),
                "Pearson_p": float(pearson_p),
                "Query_Cells": int(np.sum(query_mask)),
                "Ref_Cells": int(np.sum(ref_mask)),
                "Significant": pearson_p < 0.05
            })
            
            log_message(f"  {label} - Pearson: {pearson_r:.4f}, p-value: {pearson_p:.4e}")
        except Exception as e:
            log_message(f"  Error calculating correlation for {label}: {str(e)}")
    
    # Create correlation dataframe
    correlation_df = pd.DataFrame(correlation_df_data)
    
    # Calculate mean correlation across all labels
    if correlations:
        mean_pearson = np.mean([c["pearson_r"] for c in correlations.values() if not np.isnan(c["pearson_r"])])
        log_message(f"Mean Pearson correlation: {mean_pearson:.4f}")
        
        # For significant correlations only
        sig_correlations = [c["pearson_r"] for k, c in correlations.items() if not np.isnan(c["pearson_r"]) and c["pearson_p"] < 0.05]
        if sig_correlations:
            mean_sig_pearson = np.mean(sig_correlations)
            log_message(f"Mean significant Pearson correlation (p<0.05): {mean_sig_pearson:.4f}")
        else:
            mean_sig_pearson = np.nan
    else:
        mean_pearson = np.nan
        mean_sig_pearson = np.nan
        log_message("No valid correlations calculated")
    
    # Compile results
    results = {
        "common_labels": len(common_labels),
        "common_genes": len(common_genes),
        "consistent_cells": int(np.sum(consistent_mask)),
        "mean_pearson": float(mean_pearson),
        "mean_sig_pearson": float(mean_sig_pearson),
        "label_correlations": correlations,
        "correlation_df": correlation_df
    }
    
    return results


def visualize_label_correlations(metrics, file_name, output_dir, label_type="Cell Type"):
    """
    Visualize label correlations
    
    Parameters
    ----------
    metrics : dict
        Dictionary with correlation metrics
    file_name : str
        Base name for output files
    output_dir : str
        Directory for output files
    label_type : str
        Type of label (e.g., "Cell Type" or "Lineage")
    """
    if "error" in metrics:
        log_message(f"Error in {label_type} correlation metrics: {metrics['error']}")
        return
    
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    
    # Get correlation dataframe
    if "correlation_df" in metrics:
        df = metrics["correlation_df"]
        
        # Skip if dataframe is empty
        if len(df) == 0:
            log_message(f"No valid correlations for {label_type}")
            return
        
        # Sort by Pearson correlation
        df = df.sort_values('Pearson_r', ascending=False)
        
        # 1. Bar plot of Pearson correlations by label
        plt.figure(figsize=(12, 10))
        ax = sns.barplot(data=df, y='Label', x='Pearson_r')
        
        # Color significant correlations
        for i, p in enumerate(df['Pearson_p']):
            if p < 0.05:
                ax.get_children()[i].set_facecolor('skyblue')
            if p < 0.01:
                ax.get_children()[i].set_facecolor('steelblue')
            if p < 0.001:
                ax.get_children()[i].set_facecolor('navy')
        
        # Add query cell counts
        for i, row in df.iterrows():
            ax.text(row['Pearson_r'] + 0.02, i, f"n={row['Query_Cells']} cells", 
                   va='center', ha='left', fontsize=9)
        
        plt.title(f"{label_type} Gene Expression Correlation with Reference\n{file_name}\nUsing {metrics['consistent_cells']} consistent cells and {metrics['common_genes']} genes")
        plt.xlabel('Pearson Correlation')
        plt.axvline(x=0, color='gray', linestyle='--')
        plt.xlim(-1.1, 1.1)
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_pearson_correlation.png"), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save correlation dataframe to CSV
        df_output_path = os.path.join(output_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_pearson_correlation.csv")
        df.to_csv(df_output_path, index=False)
        log_message(f"Saved {label_type} correlation dataframe to {df_output_path}")


def process_file(query_file, reference_file, output_dir):
    """Process a single query file against the reference"""
    query_name = os.path.basename(query_file).replace(".h5ad", "")
    log_message(f"Processing {query_name}")
    
    # Load query data
    query_adata = sc.read_h5ad(query_file)
    log_message(f"Loaded query data with {query_adata.n_obs} cells and {query_adata.n_vars} genes")
    
    # Load reference data
    reference_adata = sc.read_h5ad(reference_file)
    reference_adata.obsm['X_umap']=reference_adata.obsm['X_scANVI']  # manually define the umap which will be used
    log_message(f"Loaded reference data with {reference_adata.n_obs} cells and {reference_adata.n_vars} genes")
    
    metrics = {}
    
    # Calculate cell type correlation if available
    if 'reanno_pred' in query_adata.obs:
        log_message("Calculating cell type label correlation")
        # Check if we have consistency information
        consistency_key = None  # Let the function auto-detect consistency
        
        metrics['celltype_correlation'] = calculate_matching_label_correlation(
            query_adata, reference_adata,
            query_label_key='reanno_pred',
            ref_label_key='reanno',
            consistency_key=consistency_key
        )
        
        # Visualize cell type correlation
        visualize_label_correlations(
            metrics['celltype_correlation'], 
            query_name, 
            output_dir, 
            label_type="Cell Type"
        )
    else:
        log_message("Warning: 'reanno_pred' not found in query data")
    
    # Calculate lineage correlation if available
    if 'lineage_pred' in query_adata.obs:
        log_message("Calculating lineage label correlation")
        # Check if we have consistency information
        consistency_key = None  # Let the function auto-detect consistency
        
        metrics['lineage_correlation'] = calculate_matching_label_correlation(
            query_adata, reference_adata,
            query_label_key='lineage_pred',
            ref_label_key='lineage',
            consistency_key=consistency_key
        )
        
        # Visualize lineage correlation
        visualize_label_correlations(
            metrics['lineage_correlation'], 
            query_name, 
            output_dir, 
            label_type="Lineage"
        )
    else:
        log_message("Warning: 'lineage_pred' not found in query data")
    
    # Save metrics to JSON
    import json
    metrics_file = os.path.join(output_dir, f"{query_name}_expression_correlation.json")
    
    # Convert numpy values to Python types
    def convert_for_json(obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif pd.isna(obj):
            return None
        else:
            return obj
    
    # Process metrics dictionary for JSON (exclude DataFrames)
    json_metrics = {}
    for key, value in metrics.items():
        if key in ['celltype_correlation', 'lineage_correlation']:
            json_metrics[key] = {}
            for k, v in value.items():
                if k != 'correlation_df':  # Skip DataFrame
                    if isinstance(v, dict):
                        json_metrics[key][k] = {k2: convert_for_json(v2) for k2, v2 in v.items()}
                    else:
                        json_metrics[key][k] = convert_for_json(v)
    
    with open(metrics_file, 'w') as f:
        json.dump(json_metrics, f, indent=2)
    log_message(f"Saved correlation metrics to {metrics_file}")
    
    return metrics


def main():
    # Set paths
    data_dir = '/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output' # Change this to your data directory
    reference_file = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/data/human_clustering_20250724_v3.h5ad'  # Change this to your reference file
    output_dir = './expression_correlation_metrics'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find h5ad files with consistency information
    h5ad_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) 
                 if f.endswith('_with_certainty.h5ad') or f.endswith('_scPoli_annotated.h5ad')]
    log_message(f"Found {len(h5ad_files)} files to process")
    
    # Process each file
    all_metrics = {}
    
    for file_path in h5ad_files:
        file_name = os.path.basename(file_path).replace("_with_certainty.h5ad", "").replace("_scPoli_annotated.h5ad", "")
        try:
            metrics = process_file(file_path, reference_file, output_dir)
            all_metrics[file_name] = metrics
        except Exception as e:
            log_message(f"Error processing {file_name}: {str(e)}")
            import traceback
            log_message(traceback.format_exc())
    
    # Compare all models 
    if len(all_metrics) > 1:
        # Create aggregate correlation dataframes
        all_ct_correlations = []
        all_lin_correlations = []
        
        # Collect correlation data
        for model, metrics in all_metrics.items():
            # Cell type correlations
            if 'celltype_correlation' in metrics and 'error' not in metrics['celltype_correlation']:
                if 'correlation_df' in metrics['celltype_correlation']:
                    df = metrics['celltype_correlation']['correlation_df']
                    df['Model'] = model
                    all_ct_correlations.append(df)
            
            # Lineage correlations
            if 'lineage_correlation' in metrics and 'error' not in metrics['lineage_correlation']:
                if 'correlation_df' in metrics['lineage_correlation']:
                    df = metrics['lineage_correlation']['correlation_df']
                    df['Model'] = model
                    all_lin_correlations.append(df)
        
        # Combine and save all cell type correlations
        if all_ct_correlations:
            combined_ct_df = pd.concat(all_ct_correlations, ignore_index=True)
            combined_ct_path = os.path.join(output_dir, "all_models_celltype_correlations.csv")
            combined_ct_df.to_csv(combined_ct_path, index=False)
            log_message(f"Saved combined cell type correlations to {combined_ct_path}")
        
        # Combine and save all lineage correlations
        if all_lin_correlations:
            combined_lin_df = pd.concat(all_lin_correlations, ignore_index=True)
            combined_lin_path = os.path.join(output_dir, "all_models_lineage_correlations.csv")
            combined_lin_df.to_csv(combined_lin_path, index=False)
            log_message(f"Saved combined lineage correlations to {combined_lin_path}")
        
        # Cell type correlation comparison
        ct_data = []
        for model, metrics in all_metrics.items():
            if 'celltype_correlation' in metrics and 'error' not in metrics['celltype_correlation']:
                correlation = metrics['celltype_correlation']
                ct_data.append({
                    'Model': model,
                    'Pearson r': correlation.get('mean_pearson', np.nan),
                    'Common Labels': correlation.get('common_labels', 0),
                    'Common Genes': correlation.get('common_genes', 0),
                    'Consistent Cells': correlation.get('consistent_cells', 0)
                })
        
        if ct_data:
            ct_df = pd.DataFrame(ct_data).sort_values('Pearson r', ascending=False)
            
            # Plot comparison
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(x='Model', y='Pearson r', data=ct_df)
            
            # Add consistent cell counts
            for i, row in ct_df.iterrows():
                ax.text(i, row['Pearson r'] + 0.02, f"{row['Consistent Cells']} cells", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Cell Type Expression Correlation Comparison Across Models")
            plt.ylabel('Mean Pearson Correlation')
            plt.xticks(rotation=90)
            plt.ylim(-1, 1)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "celltype_correlation_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            ct_df.to_csv(os.path.join(output_dir, "celltype_correlation_comparison.csv"), index=False)
        
        # Lineage correlation comparison
        lin_data = []
        for model, metrics in all_metrics.items():
            if 'lineage_correlation' in metrics and 'error' not in metrics['lineage_correlation']:
                correlation = metrics['lineage_correlation']
                lin_data.append({
                    'Model': model,
                    'Pearson r': correlation.get('mean_pearson', np.nan),
                    'Common Labels': correlation.get('common_labels', 0),
                    'Common Genes': correlation.get('common_genes', 0),
                    'Consistent Cells': correlation.get('consistent_cells', 0)
                })
        
        if lin_data:
            lin_df = pd.DataFrame(lin_data).sort_values('Pearson r', ascending=False)
            
            # Plot comparison
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(x='Model', y='Pearson r', data=lin_df)
            
            # Add consistent cell counts
            for i, row in lin_df.iterrows():
                ax.text(i, row['Pearson r'] + 0.02, f"{row['Consistent Cells']} cells", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Lineage Expression Correlation Comparison Across Models")
            plt.ylabel('Mean Pearson Correlation')
            plt.xticks(rotation=90)
            plt.ylim(-1, 1)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "lineage_correlation_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            lin_df.to_csv(os.path.join(output_dir, "lineage_correlation_comparison.csv"), index=False)
    
    log_message("Completed expression correlation analysis")


if __name__ == "__main__":
    main()

