#!/usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ---------- FUNCTIONS ----------

def log_message(message, log_file='marker_overlap.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")


def calculate_markers(adata, groupby, min_cells=10, certain_key=None):
    """
    Calculate marker genes for each group using Scanpy's rank_genes_groups
    
    Parameters
    ----------
    adata : AnnData
        Dataset with labels
    groupby : str
        Column in adata.obs to use for grouping (e.g., 'final_anno_pred')
    min_cells : int
        Minimum number of cells required for a group to be considered
    certain_key : str, optional
        Column in adata.obs indicating which cells are certain
        
    Returns
    -------
    dict
        Dictionary with marker genes per group
    """
    log_message(f"Calculating marker genes using {groupby}")
    
    # Make a copy of the AnnData object to avoid modifying the original
    adata_copy = adata.copy()
    
    # Filter by certainty if available
    if certain_key is not None and certain_key in adata.obs:
        certain_mask = adata.obs[certain_key].values.astype(bool)
        log_message(f"Using {np.sum(certain_mask)}/{len(certain_mask)} certain cells for marker calculation")
        if np.sum(certain_mask) == 0:
            log_message(f"Error: No cells meet certainty criteria with key '{certain_key}'")
            return {"error": "No certain cells found"}
        adata_copy = adata_copy[certain_mask]
    
    # Get group counts and filter small groups
    group_counts = adata_copy.obs[groupby].value_counts()
    valid_groups = group_counts[group_counts >= min_cells].index.tolist()
    
    if not valid_groups:
        log_message(f"No valid groups with at least {min_cells} cells after filtering")
        return {"error": "No valid groups"}
    
    # Keep only cells from valid groups
    adata_copy = adata_copy[adata_copy.obs[groupby].isin(valid_groups)]

    
    # Run marker gene calculation
    try:    
        sc.tl.rank_genes_groups(
            adata_copy, 
            groupby=groupby, 
            method='wilcoxon', 
            key_added='rank_genes', 
            use_raw=True if 'raw' in adata_copy.layers else False
        )
        
        # Extract marker genes for each group
        results = {}
        
        for group in valid_groups:
            # Get marker genes with adjusted p-values
            markers_df = sc.get.rank_genes_groups_df(adata_copy, key='rank_genes', group=group)
            
            # Filter for significant markers
            sig_markers = markers_df[(markers_df['pvals_adj'] < 0.05) & (markers_df['logfoldchanges'] > 0.25)]
            
            # Store top markers
            results[group] = {
                'markers': sig_markers['names'].tolist(),
                'scores': sig_markers['scores'].tolist(),
                'logfoldchanges': sig_markers['logfoldchanges'].tolist(),
                'pvals_adj': sig_markers['pvals_adj'].tolist(),
                'n_cells': int(group_counts[group])
            }
            
            log_message(f"  {group}: {len(results[group]['markers'])} significant marker genes from {int(group_counts[group])} cells")
        
        return results
        
    except Exception as e:
        log_message(f"Error in marker gene calculation: {str(e)}")
        import traceback
        log_message(traceback.format_exc())
        return {"error": str(e)}


def calculate_marker_overlap(query_markers, reference_markers, all_genes=None):
    """
    Calculate overlap between marker genes in query and reference datasets
    
    Parameters
    ----------
    query_markers : dict
        Dictionary with marker genes per group from query dataset
    reference_markers : dict
        Dictionary with marker genes per group from reference dataset
    all_genes : list, optional
        List of all genes in the analysis (for hypergeometric test)
        
    Returns
    -------
    dict
        Dictionary with overlap metrics
    """
    log_message("Calculating marker gene overlap between query and reference")
    
    if "error" in query_markers or "error" in reference_markers:
        error_msg = query_markers.get("error", reference_markers.get("error", "Unknown error"))
        log_message(f"Error in marker calculation: {error_msg}")
        return {"error": error_msg}
    
    # Get common groups between query and reference
    query_groups = set(query_markers.keys())
    ref_groups = set(reference_markers.keys())
    common_groups = query_groups.intersection(ref_groups)
    
    log_message(f"Found {len(common_groups)} groups common to both query and reference")
    
    if not common_groups:
        log_message("No common groups found between query and reference")
        return {"error": "No common groups"}
    
    # Create a universe of all genes if not provided
    if all_genes is None:
        all_genes_set = set()
        for group in query_markers:
            all_genes_set.update(query_markers[group]['markers'])
        for group in reference_markers:
            all_genes_set.update(reference_markers[group]['markers'])
        all_genes = list(all_genes_set)
    
    # Calculate overlap for each common group
    overlap_metrics = {}
    overlap_df_data = []
    
    for group in common_groups:
        # Get marker sets
        query_set = set(query_markers[group]['markers'])
        ref_set = set(reference_markers[group]['markers'])
        
        # Calculate overlap
        overlap_set = query_set.intersection(ref_set)
        
        # Calculate different overlap metrics
        overlap_size = len(overlap_set)
        precision = overlap_size / len(query_set) if len(query_set) > 0 else 0
        recall = overlap_size / len(ref_set) if len(ref_set) > 0 else 0
        f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        jaccard = overlap_size / len(query_set.union(ref_set)) if len(query_set.union(ref_set)) > 0 else 0
        
        # Perform Fisher's exact test
        contingency_table = [
            [overlap_size, len(ref_set) - overlap_size], 
            [len(query_set) - overlap_size, len(all_genes) - len(query_set) - len(ref_set) + overlap_size]
        ]
        
        try:
            odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
        except:
            odds_ratio, p_value = np.nan, np.nan
            
        # Store metrics
        overlap_metrics[group] = {
            'query_markers': len(query_set),
            'ref_markers': len(ref_set),
            'overlap': overlap_size,
            'precision': precision,
            'recall': recall,
            'f1_score': f1_score,
            'jaccard': jaccard,
            'odds_ratio': float(odds_ratio),
            'p_value': float(p_value),
            'query_only': len(query_set - ref_set),
            'ref_only': len(ref_set - query_set),
            'overlap_genes': list(overlap_set)
        }
        
        # Add to dataframe data
        overlap_df_data.append({
            'Group': group,
            'Query_Markers': len(query_set),
            'Ref_Markers': len(ref_set),
            'Overlap': overlap_size,
            'Precision': precision,
            'Recall': recall,
            'F1_Score': f1_score,
            'Jaccard': jaccard,
            'Odds_Ratio': float(odds_ratio),
            'P_Value': float(p_value),
            'Overlap_Percentage': 100 * precision,
            'Query_Cells': query_markers[group]['n_cells']
        })
        
        log_message(f"  {group}: {overlap_size} overlapping markers out of {len(query_set)} query markers ({precision:.2%})")
    
    # Create overlap dataframe
    overlap_df = pd.DataFrame(overlap_df_data)
    
    # Apply multiple testing correction
    if len(overlap_df) > 1:
        _, adjusted_p_values, _, _ = multipletests(overlap_df['P_Value'].fillna(1).values, method='fdr_bh')
        overlap_df['Adjusted_P_Value'] = adjusted_p_values
        overlap_df['Significant'] = overlap_df['Adjusted_P_Value'] < 0.05
    else:
        overlap_df['Adjusted_P_Value'] = overlap_df['P_Value']
        overlap_df['Significant'] = overlap_df['P_Value'] < 0.05
    
    # Calculate overall metrics
    mean_precision = np.mean(overlap_df['Precision'])
    mean_recall = np.mean(overlap_df['Recall'])
    mean_f1 = np.mean(overlap_df['F1_Score'])
    mean_jaccard = np.mean(overlap_df['Jaccard'])
    total_overlap = np.sum(overlap_df['Overlap'])
    total_query_markers = np.sum(overlap_df['Query_Markers'])
    total_ref_markers = np.sum(overlap_df['Ref_Markers'])
    overall_precision = total_overlap / total_query_markers if total_query_markers > 0 else 0
    overall_recall = total_overlap / total_ref_markers if total_ref_markers > 0 else 0
    
    overall_metrics = {
        'mean_precision': float(mean_precision),
        'mean_recall': float(mean_recall),
        'mean_f1': float(mean_f1),
        'mean_jaccard': float(mean_jaccard),
        'weighted_precision': float(np.average(overlap_df['Precision'], weights=overlap_df['Query_Markers'])),
        'weighted_recall': float(np.average(overlap_df['Recall'], weights=overlap_df['Ref_Markers'])),
        'total_overlap': int(total_overlap),
        'total_query_markers': int(total_query_markers),
        'total_ref_markers': int(total_ref_markers),
        'overall_precision': float(overall_precision),
        'overall_recall': float(overall_recall),
        'common_groups': len(common_groups)
    }
    
    log_message(f"Overall marker overlap: {overall_precision:.2%} precision, {overall_recall:.2%} recall, {mean_f1:.2%} mean F1")
    
    # Return all metrics and DataFrame
    return {
        'overlap_df': overlap_df,
        'overlap_per_group': overlap_metrics,
        'overall_metrics': overall_metrics
    }


def visualize_marker_overlap(metrics, file_name, output_dir, label_type="Cell Type"):
    """
    Create visualizations of marker gene overlap
    
    Parameters
    ----------
    metrics : dict
        Dictionary with overlap metrics
    file_name : str
        Base name for output files
    output_dir : str
        Directory for output files
    label_type : str
        Type of label (Cell Type or Lineage)
    """
    if "error" in metrics:
        log_message(f"Error in {label_type} overlap metrics: {metrics['error']}")
        return
    
    if 'overlap_df' not in metrics:
        log_message(f"No overlap dataframe in {label_type} metrics")
        return
    
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    
    # Get overlap dataframe
    df = metrics['overlap_df']
    
    # 1. Bar plot of precision by group (sorted)
    plt.figure(figsize=(12, 10))
    plot_df = df.sort_values('Precision', ascending=False)
    
    # Limit to top 30 groups for readability
    if len(plot_df) > 30:
        plot_df = plot_df.head(30)
        
    ax = sns.barplot(data=plot_df, y='Group', x='Precision', palette='viridis')
    
    # Add count labels
    for i, row in plot_df.iterrows():
        label_text = f"{row['Overlap']}/{row['Query_Markers']} markers"
        ax.text(row['Precision'] + 0.02, i, label_text, va='center', ha='left', fontsize=9)
    
    plt.title(f"{label_type} Marker Overlap Precision with Reference\n{file_name}")
    plt.xlabel('Precision (Overlap / Query Markers)')
    plt.xlim(0, 1.1)
    plt.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_marker_precision.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Scatter plot of precision vs. number of markers
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=df, x='Query_Markers', y='Precision', hue='Query_Cells', 
                    size='Overlap', sizes=(20, 500), alpha=0.7, palette='viridis')
    
    # Add group labels for important points
    top_groups = df.sort_values('Precision', ascending=False).head(5)
    for _, row in top_groups.iterrows():
        plt.annotate(row['Group'], (row['Query_Markers'], row['Precision']), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    largest_groups = df.sort_values('Query_Markers', ascending=False).head(5)
    for _, row in largest_groups.iterrows():
        if row['Group'] not in top_groups['Group'].values:
            plt.annotate(row['Group'], (row['Query_Markers'], row['Precision']), 
                        xytext=(5, -10), textcoords='offset points', fontsize=9)
    
    plt.title(f"{label_type} Marker Precision vs. Number of Markers\n{file_name}")
    plt.xlabel('Number of Query Markers')
    plt.ylabel('Precision (Overlap / Query Markers)')
    plt.ylim(0, 1.05)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_marker_precision_vs_count.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Overall metrics plot
    overall = metrics['overall_metrics']
    metrics_to_plot = ['overall_precision', 'mean_precision', 'weighted_precision', 
                      'overall_recall', 'mean_recall', 'weighted_recall', 'mean_f1', 'mean_jaccard']
    
    plt.figure(figsize=(12, 6))
    plot_data = pd.DataFrame({
        'Metric': [m.replace('_', ' ').title() for m in metrics_to_plot],
        'Value': [overall[m] for m in metrics_to_plot]
    })
    
    sns.barplot(data=plot_data, x='Metric', y='Value')
    plt.title(f"Overall {label_type} Marker Overlap Metrics\n{file_name}")
    plt.ylabel('Score')
    plt.ylim(0, 1.05)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_overall_metrics.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()


def process_file(query_file, reference_file, output_dir):
    """Process a single query file against the reference"""
    query_name = os.path.basename(query_file).replace(".h5ad", "")
    log_message(f"Processing {query_name}")
    
    # Load query data
    query_adata = sc.read_h5ad(query_file)
    log_message(f"Loaded query data with {query_adata.n_obs} cells and {query_adata.n_vars} genes")
    
    # Load reference data
    reference_adata = sc.read_h5ad(reference_file)
    log_message(f"Loaded reference data with {reference_adata.n_obs} cells and {reference_adata.n_vars} genes")
    
    # Get common genes (universe for calculations)
    common_genes = np.intersect1d(query_adata.var_names, reference_adata.var_names)
    log_message(f"Found {len(common_genes)} common genes between query and reference")
    
    # Subset both datasets to common genes
    query_adata = query_adata[:, common_genes].copy()
    reference_adata = reference_adata[:, common_genes].copy()
    
    metrics = {}
    
    # Cell type marker overlap
    if 'final_anno_pred' in query_adata.obs and 'final_anno' in reference_adata.obs:
        log_message("Calculating cell type marker overlap")
        
        # Get certain cells mask for query dataset
        certainty_key = 'is_celltype_certain' if 'is_celltype_certain' in query_adata.obs else None
        
        # Calculate markers for query dataset
        query_markers = calculate_markers(
            query_adata, 
            groupby='final_anno_pred', 
            min_cells=10, 
            certain_key=certainty_key
        )
        
        # Calculate markers for reference dataset
        reference_markers = calculate_markers(
            reference_adata, 
            groupby='final_anno', 
            min_cells=10
        )
        
        # Calculate overlap metrics
        metrics['celltype_marker_overlap'] = calculate_marker_overlap(
            query_markers, 
            reference_markers, 
            all_genes=common_genes
        )
        
        # Visualize overlap metrics
        visualize_marker_overlap(
            metrics['celltype_marker_overlap'], 
            query_name, 
            output_dir, 
            label_type="Cell Type"
        )
        
        # Save cell type overlap dataframe
        if 'overlap_df' in metrics['celltype_marker_overlap']:
            ct_df_path = os.path.join(output_dir, f"{query_name}_celltype_marker_overlap.csv")
            metrics['celltype_marker_overlap']['overlap_df'].to_csv(ct_df_path, index=False)
            log_message(f"Saved cell type marker overlap dataframe to {ct_df_path}")
    else:
        if 'final_anno_pred' not in query_adata.obs:
            log_message("Warning: 'final_anno_pred' not found in query data")
        if 'final_anno' not in reference_adata.obs:
            log_message("Warning: 'final_anno' not found in reference data")
    
    # Lineage marker overlap
    if 'final_lineage_pred' in query_adata.obs and 'final_lineage' in reference_adata.obs:
        log_message("Calculating lineage marker overlap")
        
        # Get certain cells mask for query dataset
        certainty_key = 'is_lineage_certain' if 'is_lineage_certain' in query_adata.obs else None
        
        # Calculate markers for query dataset
        query_markers = calculate_markers(
            query_adata, 
            groupby='final_lineage_pred', 
            min_cells=10, 
            certain_key=certainty_key
        )
        
        # Calculate markers for reference dataset
        reference_markers = calculate_markers(
            reference_adata, 
            groupby='final_lineage', 
            min_cells=10
        )
        
        # Calculate overlap metrics
        metrics['lineage_marker_overlap'] = calculate_marker_overlap(
            query_markers, 
            reference_markers, 
            all_genes=common_genes
        )
        
        # Visualize overlap metrics
        visualize_marker_overlap(
            metrics['lineage_marker_overlap'], 
            query_name, 
            output_dir, 
            label_type="Lineage"
        )
        
        # Save lineage overlap dataframe
        if 'overlap_df' in metrics['lineage_marker_overlap']:
            lin_df_path = os.path.join(output_dir, f"{query_name}_lineage_marker_overlap.csv")
            metrics['lineage_marker_overlap']['overlap_df'].to_csv(lin_df_path, index=False)
            log_message(f"Saved lineage marker overlap dataframe to {lin_df_path}")
    else:
        if 'final_lineage_pred' not in query_adata.obs:
            log_message("Warning: 'final_lineage_pred' not found in query data")
        if 'final_lineage' not in reference_adata.obs:
            log_message("Warning: 'final_lineage' not found in reference data")
    
    # Save metrics to JSON
    import json
    metrics_file = os.path.join(output_dir, f"{query_name}_marker_overlap.json")
    
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
    
    for key in ['celltype_marker_overlap', 'lineage_marker_overlap']:
        if key in metrics:
            json_metrics[key] = {}
            for k, v in metrics[key].items():
                if k != 'overlap_df':  # Skip DataFrame
                    if isinstance(v, dict):
                        json_metrics[key][k] = {}
                        for k2, v2 in v.items():
                            if isinstance(v2, dict):
                                json_metrics[key][k][k2] = {k3: convert_for_json(v3) for k3, v3 in v2.items()}
                            else:
                                json_metrics[key][k][k2] = convert_for_json(v2)
                    else:
                        json_metrics[key][k] = convert_for_json(v)
    
    with open(metrics_file, 'w') as f:
        json.dump(json_metrics, f, indent=2)
    
    log_message(f"Saved marker overlap metrics to {metrics_file}")
    
    return metrics


def main():
    # Set paths
    data_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/output'  # Change this to your data directory
    reference_file = '/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/human_reanno_20250108.h5ad'  # Change this to your reference file
    output_dir = './marker_overlap_metrics'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find h5ad files with certainty information
    h5ad_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) 
                 if f.endswith('_scPoli_annotated.h5ad')]
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
        # Create aggregate marker overlap dataframes
        all_ct_overlaps = []
        all_lin_overlaps = []
        
        # Collect overlap data
        for model, metrics in all_metrics.items():
            # Cell type overlaps
            if 'celltype_marker_overlap' in metrics and 'overlap_df' in metrics['celltype_marker_overlap']:
                df = metrics['celltype_marker_overlap']['overlap_df'].copy()
                df['Model'] = model
                all_ct_overlaps.append(df)
            
            # Lineage overlaps
            if 'lineage_marker_overlap' in metrics and 'overlap_df' in metrics['lineage_marker_overlap']:
                df = metrics['lineage_marker_overlap']['overlap_df'].copy()
                df['Model'] = model
                all_lin_overlaps.append(df)
        
        # Combine and save all cell type overlaps
        if all_ct_overlaps:
            combined_ct_df = pd.concat(all_ct_overlaps, ignore_index=True)
            combined_ct_path = os.path.join(output_dir, "all_models_celltype_marker_overlap.csv")
            combined_ct_df.to_csv(combined_ct_path, index=False)
            log_message(f"Saved combined cell type marker overlaps to {combined_ct_path}")
        
        # Combine and save all lineage overlaps
        if all_lin_overlaps:
            combined_lin_df = pd.concat(all_lin_overlaps, ignore_index=True)
            combined_lin_path = os.path.join(output_dir, "all_models_lineage_marker_overlap.csv")
            combined_lin_df.to_csv(combined_lin_path, index=False)
            log_message(f"Saved combined lineage marker overlaps to {combined_lin_path}")
        
        # Cell type marker overlap comparison
        ct_data = []
        for model, metrics in all_metrics.items():
            if 'celltype_marker_overlap' in metrics and 'overall_metrics' in metrics['celltype_marker_overlap']:
                overlap = metrics['celltype_marker_overlap']['overall_metrics']
                ct_data.append({
                    'Model': model,
                    'Precision': overlap.get('mean_precision', np.nan) * 100,  # Convert to percentage
                    'Weighted Precision': overlap.get('weighted_precision', np.nan) * 100,
                    'F1 Score': overlap.get('mean_f1', np.nan) * 100,
                    'Total Overlap': overlap.get('total_overlap', 0),
                    'Total Markers': overlap.get('total_query_markers', 0),
                    'Common Groups': overlap.get('common_groups', 0)
                })
        
        if ct_data:
            ct_df = pd.DataFrame(ct_data).sort_values('Precision', ascending=False)
            
            # Plot comparison
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Precision', data=ct_df)
            
            # Add count labels
            for i, row in ct_df.iterrows():
                ax.text(i, row['Precision'] + 1, 
                       f"{row['Total Overlap']}/{row['Total Markers']}", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Cell Type Marker Overlap Comparison Across Models")
            plt.ylabel('Mean Precision (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "celltype_marker_overlap_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            ct_df.to_csv(os.path.join(output_dir, "celltype_marker_overlap_comparison.csv"), index=False)
        
        # Lineage marker overlap comparison
        lin_data = []
        for model, metrics in all_metrics.items():
            if 'lineage_marker_overlap' in metrics and 'overall_metrics' in metrics['lineage_marker_overlap']:
                overlap = metrics['lineage_marker_overlap']['overall_metrics']
                lin_data.append({
                    'Model': model,
                    'Precision': overlap.get('mean_precision', np.nan) * 100,  # Convert to percentage
                    'Weighted Precision': overlap.get('weighted_precision', np.nan) * 100,
                    'F1 Score': overlap.get('mean_f1', np.nan) * 100,
                    'Total Overlap': overlap.get('total_overlap', 0),
                    'Total Markers': overlap.get('total_query_markers', 0),
                    'Common Groups': overlap.get('common_groups', 0)
                })
        
        if lin_data:
            lin_df = pd.DataFrame(lin_data).sort_values('Precision', ascending=False)
            
            # Plot comparison
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Precision', data=lin_df)
            
            # Add count labels
            for i, row in lin_df.iterrows():
                ax.text(i, row['Precision'] + 1, 
                       f"{row['Total Overlap']}/{row['Total Markers']}", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Lineage Marker Overlap Comparison Across Models")
            plt.ylabel('Mean Precision (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "lineage_marker_overlap_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            lin_df.to_csv(os.path.join(output_dir, "lineage_marker_overlap_comparison.csv"), index=False)
    
    log_message("Completed marker overlap analysis")


if __name__ == "__main__":
    main()