#!/usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# ---------- FUNCTIONS ----------

def log_message(message, log_file='mean_certainty_metrics_certain_cells.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

def calculate_mean_certainty_all_cells(adata):
    """
    Calculate mean certainty (1 - uncertainty) for each label, 
    using ALL cells instead of only cells marked as "certain"
    
    Parameters
    ----------
    adata : AnnData
        Annotated dataset with uncertainty scores
        
    Returns
    -------
    dict
        Dictionary with mean certainty metrics for lineage and cell type (using all cells)
    """
    metrics = {}
    
    # Check if cell type predictions and uncertainty exist
    if 'final_anno_pred' in adata.obs and 'final_anno_uncert' in adata.obs:
        # Calculate overall mean certainty for ALL cells
        overall_ct_certainty = 1 - np.mean(adata.obs['final_anno_uncert'])
        
        # Calculate per cell type mean certainty (for ALL cells)
        ct_mean_certainty = {}
        for ct in adata.obs['final_anno_pred'].unique():
            ct_mask = (adata.obs['final_anno_pred'] == ct)
            if np.sum(ct_mask) > 0:
                mean_cert = 1 - np.mean(adata.obs.loc[ct_mask, 'final_anno_uncert'])
                ct_mean_certainty[ct] = {
                    'mean_certainty': float(mean_cert),
                    'cell_count': int(np.sum(ct_mask))
                }
        
        # Store metrics
        metrics['celltype_mean_certainty'] = {
            'overall': float(overall_ct_certainty),
            'per_label': ct_mean_certainty,
            'total_cells': adata.n_obs
        }
    else:
        metrics['celltype_mean_certainty'] = None
        print("Warning: Cell type predictions or uncertainty not found")
    
    # Check if lineage predictions and uncertainty exist
    if 'final_lineage_pred' in adata.obs and 'final_lineage_uncert' in adata.obs:
        # Calculate overall mean certainty for ALL cells
        overall_lin_certainty = 1 - np.mean(adata.obs['final_lineage_uncert'])
        
        # Calculate per lineage mean certainty (for ALL cells)
        lin_mean_certainty = {}
        for lin in adata.obs['final_lineage_pred'].unique():
            lin_mask = (adata.obs['final_lineage_pred'] == lin)
            if np.sum(lin_mask) > 0:
                mean_cert = 1 - np.mean(adata.obs.loc[lin_mask, 'final_lineage_uncert'])
                lin_mean_certainty[lin] = {
                    'mean_certainty': float(mean_cert),
                    'cell_count': int(np.sum(lin_mask))
                }
        
        # Store metrics
        metrics['lineage_mean_certainty'] = {
            'overall': float(overall_lin_certainty),
            'per_label': lin_mean_certainty,
            'total_cells': adata.n_obs
        }
    else:
        metrics['lineage_mean_certainty'] = None
        print("Warning: Lineage predictions or uncertainty not found")
    
    return metrics

def process_file(file_path, output_dir):
    """Process a single h5ad file to calculate mean certainty for all cells"""
    file_name = os.path.basename(file_path).replace(".h5ad", "")
    print(f"Processing {file_name}")
    
    # Load annotated data
    try:
        adata = sc.read_h5ad(file_path)
        print(f"Loaded {file_name} with {adata.n_obs} cells")
    except Exception as e:
        print(f"Error loading {file_path}: {str(e)}")
        return None
    
    # Calculate mean certainty metrics for ALL cells
    metrics = calculate_mean_certainty_all_cells(adata)
    
    # Create directory for individual model results
    model_dir = os.path.join(output_dir, file_name)
    os.makedirs(model_dir, exist_ok=True)
    
    # Save per cell type mean certainty to CSV
    if metrics['celltype_mean_certainty'] is not None:
        ct_data = []
        for ct, values in metrics['celltype_mean_certainty']['per_label'].items():
            ct_data.append({
                'Cell_Type': ct,
                'Mean_Certainty': values['mean_certainty'],
                'Cell_Count': values['cell_count']
            })
        
        if ct_data:
            ct_df = pd.DataFrame(ct_data)
            ct_df = ct_df.sort_values('Mean_Certainty', ascending=False)
            ct_csv_path = os.path.join(model_dir, f"{file_name}_celltype_mean_certainty_all_cells.csv")
            ct_df.to_csv(ct_csv_path, index=False)
            print(f"Saved cell type mean certainty to {ct_csv_path}")
            
            # Visualize cell type mean certainty (for all cells)
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(y='Cell_Type', x='Mean_Certainty', data=ct_df)
            
            # Add cell count
            for i, row in ct_df.iterrows():
                ax.text(row['Mean_Certainty'] + 0.01, i, f"n={row['Cell_Count']}", 
                       va='center', ha='left', fontsize=8)
            
            plt.title(f"Mean Certainty by Cell Type (All Cells)\n"
                      f"{file_name} - {adata.n_obs} cells")
            plt.xlabel('Mean Certainty (1 - Uncertainty)')
            plt.xlim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(model_dir, f"{file_name}_celltype_mean_certainty_all_cells.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    # Save per lineage mean certainty to CSV
    if metrics['lineage_mean_certainty'] is not None:
        lin_data = []
        for lin, values in metrics['lineage_mean_certainty']['per_label'].items():
            lin_data.append({
                'Lineage': lin,
                'Mean_Certainty': values['mean_certainty'],
                'Cell_Count': values['cell_count']
            })
        
        if lin_data:
            lin_df = pd.DataFrame(lin_data)
            lin_df = lin_df.sort_values('Mean_Certainty', ascending=False)
            lin_csv_path = os.path.join(model_dir, f"{file_name}_lineage_mean_certainty_all_cells.csv")
            lin_df.to_csv(lin_csv_path, index=False)
            print(f"Saved lineage mean certainty to {lin_csv_path}")
            
            # Visualize lineage mean certainty (for all cells)
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(y='Lineage', x='Mean_Certainty', data=lin_df)
            
            # Add cell count
            for i, row in lin_df.iterrows():
                ax.text(row['Mean_Certainty'] + 0.01, i, f"n={row['Cell_Count']}", 
                       va='center', ha='left', fontsize=8)
            
            plt.title(f"Mean Certainty by Lineage (All Cells)\n"
                      f"{file_name} - {adata.n_obs} cells")
            plt.xlabel('Mean Certainty (1 - Uncertainty)')
            plt.xlim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(model_dir, f"{file_name}_lineage_mean_certainty_all_cells.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    return metrics


def main():
    # Set paths
    data_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/output'
    output_dir = './mean_certainty_metrics_certain_cells'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all relevant h5ad files (either annotated or with_certainty)
    all_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) 
                if f.endswith('_scPoli_annotated.h5ad') or f.endswith('_with_certainty.h5ad')]
    
    # Prioritize _with_certainty files as they should have the certainty flags
    h5ad_files = []
    processed_models = set()
    
    # First add with_certainty files
    for file_path in all_files:
        if file_path.endswith('_with_certainty.h5ad'):
            model_name = os.path.basename(file_path).replace("_with_certainty.h5ad", "")
            if model_name not in processed_models:
                h5ad_files.append(file_path)
                processed_models.add(model_name)
    
    # Then add annotated files if we don't have them yet
    for file_path in all_files:
        if file_path.endswith('_scPoli_annotated.h5ad'):
            model_name = os.path.basename(file_path).replace("_scPoli_annotated.h5ad", "")
            if model_name not in processed_models:
                h5ad_files.append(file_path)
                processed_models.add(model_name)
    
    log_message(f"Found {len(h5ad_files)} files to process")
    
    # Process each file
    all_metrics = {}
    
    for file_path in h5ad_files:
        file_name = os.path.basename(file_path).replace("_with_certainty.h5ad", "").replace("_scPoli_annotated.h5ad", "")
        try:
            metrics = process_file(file_path, output_dir)
            if metrics:
                all_metrics[file_name] = metrics
        except Exception as e:
            log_message(f"Error processing {file_name}: {str(e)}")
            import traceback
            log_message(traceback.format_exc())
    
    # Create comparative analysis across all models
    if len(all_metrics) > 1:
        # Cell type overall mean certainty comparison
        ct_overall_data = []
        for model, metrics in all_metrics.items():
            if metrics and metrics['celltype_mean_certainty'] is not None:
                ct_overall_data.append({
                    'Model': model,
                    'Mean_Certainty': metrics['celltype_mean_certainty']['overall'],
                    'Certain_Cells': metrics['celltype_mean_certainty']['total_certain_cells'],
                    'Total_Cells': metrics['celltype_mean_certainty']['total_cells'],
                    'Certain_Percentage': metrics['celltype_mean_certainty']['certain_percentage']
                })
        
        if ct_overall_data:
            ct_overall_df = pd.DataFrame(ct_overall_data)
            ct_overall_df = ct_overall_df.sort_values('Mean_Certainty', ascending=False)
            ct_overall_path = os.path.join(output_dir, "celltype_overall_mean_certainty_certain_cells_comparison.csv")
            ct_overall_df.to_csv(ct_overall_path, index=False)
            log_message(f"Saved cell type overall mean certainty comparison to {ct_overall_path}")
            
            # Plot comparison
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(x='Model', y='Mean_Certainty', data=ct_overall_df)
            
            # Add cell count and percentage
            for i, row in ct_overall_df.iterrows():
                ax.text(i, row['Mean_Certainty'] + 0.01, 
                        f"{row['Certain_Cells']} cells\n({row['Certain_Percentage']:.1f}%)", 
                        ha='center', va='bottom', fontsize=8)
            
            plt.title("Overall Cell Type Mean Certainty Comparison (Certain Cells Only)")
            plt.ylabel('Mean Certainty (1 - Uncertainty)')
            plt.xticks(rotation=90)
            plt.ylim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "celltype_overall_mean_certainty_certain_cells_comparison.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
        
        # Lineage overall mean certainty comparison
        lin_overall_data = []
        for model, metrics in all_metrics.items():
            if metrics and metrics['lineage_mean_certainty'] is not None:
                lin_overall_data.append({
                    'Model': model,
                    'Mean_Certainty': metrics['lineage_mean_certainty']['overall'],
                    'Certain_Cells': metrics['lineage_mean_certainty']['total_certain_cells'],
                    'Total_Cells': metrics['lineage_mean_certainty']['total_cells'],
                    'Certain_Percentage': metrics['lineage_mean_certainty']['certain_percentage']
                })
        
        if lin_overall_data:
            lin_overall_df = pd.DataFrame(lin_overall_data)
            lin_overall_df = lin_overall_df.sort_values('Mean_Certainty', ascending=False)
            lin_overall_path = os.path.join(output_dir, "lineage_overall_mean_certainty_certain_cells_comparison.csv")
            lin_overall_df.to_csv(lin_overall_path, index=False)
            log_message(f"Saved lineage overall mean certainty comparison to {lin_overall_path}")
            
            # Plot comparison
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(x='Model', y='Mean_Certainty', data=lin_overall_df)
            
            # Add cell count and percentage
            for i, row in lin_overall_df.iterrows():
                ax.text(i, row['Mean_Certainty'] + 0.01, 
                        f"{row['Certain_Cells']} cells\n({row['Certain_Percentage']:.1f}%)", 
                        ha='center', va='bottom', fontsize=8)
            
            plt.title("Overall Lineage Mean Certainty Comparison (Certain Cells Only)")
            plt.ylabel('Mean Certainty (1 - Uncertainty)')
            plt.xticks(rotation=90)
            plt.ylim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "lineage_overall_mean_certainty_certain_cells_comparison.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
        
        # Compare per-label mean certainty across models
        # For cell types
        all_celltypes = set()
        for model, metrics in all_metrics.items():
            if metrics and metrics['celltype_mean_certainty'] is not None:
                all_celltypes.update(metrics['celltype_mean_certainty']['per_label'].keys())
        
        if all_celltypes:
            # Create a dataframe with models as columns and labels as rows
            ct_comparison = defaultdict(dict)
            for celltype in all_celltypes:
                for model, metrics in all_metrics.items():
                    if (metrics and metrics['celltype_mean_certainty'] is not None and 
                        celltype in metrics['celltype_mean_certainty']['per_label']):
                        ct_comparison[celltype][model] = metrics['celltype_mean_certainty']['per_label'][celltype]['mean_certainty']
                    else:
                        ct_comparison[celltype][model] = np.nan
            
            ct_comp_df = pd.DataFrame(ct_comparison).T
            ct_comp_df.index.name = 'Cell_Type'
            
            # Add a column with mean across models
            ct_comp_df['Mean_Across_Models'] = ct_comp_df.mean(axis=1)
            
            # Sort by mean certainty across models
            ct_comp_df = ct_comp_df.sort_values('Mean_Across_Models', ascending=False)
            
            # Save to CSV
            ct_comp_path = os.path.join(output_dir, "celltype_mean_certainty_by_label_certain_cells_comparison.csv")
            ct_comp_df.to_csv(ct_comp_path)
            log_message(f"Saved cell type mean certainty comparison by label to {ct_comp_path}")
        
        # For lineages
        all_lineages = set()
        for model, metrics in all_metrics.items():
            if metrics and metrics['lineage_mean_certainty'] is not None:
                all_lineages.update(metrics['lineage_mean_certainty']['per_label'].keys())
        
        if all_lineages:
            # Create a dataframe with models as columns and labels as rows
            lin_comparison = defaultdict(dict)
            for lineage in all_lineages:
                for model, metrics in all_metrics.items():
                    if (metrics and metrics['lineage_mean_certainty'] is not None and 
                        lineage in metrics['lineage_mean_certainty']['per_label']):
                        lin_comparison[lineage][model] = metrics['lineage_mean_certainty']['per_label'][lineage]['mean_certainty']
                    else:
                        lin_comparison[lineage][model] = np.nan
            
            lin_comp_df = pd.DataFrame(lin_comparison).T
            lin_comp_df.index.name = 'Lineage'
            
            # Add a column with mean across models
            lin_comp_df['Mean_Across_Models'] = lin_comp_df.mean(axis=1)
            
            # Sort by mean certainty across models
            lin_comp_df = lin_comp_df.sort_values('Mean_Across_Models', ascending=False)
            
            # Save to CSV
            lin_comp_path = os.path.join(output_dir, "lineage_mean_certainty_by_label_certain_cells_comparison.csv")
            lin_comp_df.to_csv(lin_comp_path)
            log_message(f"Saved lineage mean certainty comparison by label to {lin_comp_path}")
    
    log_message("Mean certainty analysis (for certain cells only) completed")


if __name__ == "__main__":
    main()