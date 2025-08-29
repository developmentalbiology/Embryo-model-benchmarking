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

def log_message(message, log_file='mean_certainty_metrics_consistent_cells.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")


def calculate_mean_certainty_consistent_cells(adata):
    """
    Calculate mean certainty (1 - uncertainty) for each label, 
    but only using cells that are marked as "consistent" (reanno_pred_lineage == lineage_pred)
    
    Parameters
    ----------
    adata : AnnData
        Annotated dataset with uncertainty scores and consistency information
        
    Returns
    -------
    dict
        Dictionary with mean certainty metrics for lineage and cell type (consistent cells only)
    """
    metrics = {}
    
    # Check for consistency based on reanno_pred_lineage == lineage_pred
    consistent_mask = None
    if 'reanno_pred_lineage' in adata.obs and 'lineage_pred' in adata.obs:
        # Convert categorical columns to string to avoid category comparison issues
        reanno_pred_lineage_str = adata.obs['reanno_pred_lineage'].astype(str)
        lineage_pred_str = adata.obs['lineage_pred'].astype(str)
        consistent_mask = (reanno_pred_lineage_str == lineage_pred_str)
        log_message(f"Using consistency based on reanno_pred_lineage==lineage_pred: {np.sum(consistent_mask)}/{len(consistent_mask)} cells")
    else:
        # Fallback to old certainty logic if consistency columns not available
        log_message("Consistency columns not found, falling back to certainty flags")
    
    # Check if cell type predictions and uncertainty exist
    if ('reanno_pred' in adata.obs and 'reanno_uncert' in adata.obs):
        
        # Determine which cells to use
        if consistent_mask is not None:
            # Use consistent cells
            mask_to_use = consistent_mask
            mask_name = "consistent"
        elif 'is_celltype_certain' in adata.obs:
            # Fallback to certain cells
            mask_to_use = adata.obs['is_celltype_certain'].astype(bool)
            mask_name = "certain"
            log_message("Fallback: Using 'is_celltype_certain' for cell type analysis")
        else:
            # Use all cells
            mask_to_use = np.ones(len(adata), dtype=bool)
            mask_name = "all"
            log_message("No consistency or certainty filter found, using all cells for cell type analysis")
        
        selected_cells = np.sum(mask_to_use)
        
        if selected_cells > 0:
            # Calculate overall mean certainty for selected cells only
            overall_ct_certainty = 1 - np.mean(adata.obs.loc[mask_to_use, 'reanno_uncert'])
            
            # Calculate per cell type mean certainty (for selected cells only)
            ct_mean_certainty = {}
            for ct in adata.obs.loc[mask_to_use, 'reanno_pred'].unique():
                ct_mask = (adata.obs['reanno_pred'] == ct) & mask_to_use
                if np.sum(ct_mask) > 0:
                    mean_cert = 1 - np.mean(adata.obs.loc[ct_mask, 'reanno_uncert'])
                    ct_mean_certainty[ct] = {
                        'mean_certainty': float(mean_cert),
                        'cell_count': int(np.sum(ct_mask))
                    }
            
            metrics['celltype_mean_certainty'] = {
                'overall': float(overall_ct_certainty),
                'per_label': ct_mean_certainty,
                'total_selected_cells': int(selected_cells),
                'total_cells': adata.n_obs,
                'selected_percentage': 100 * (selected_cells / adata.n_obs),
                'selection_criteria': mask_name
            }
        else:
            log_message(f"Warning: No {mask_name} cells found for cell type analysis")
            metrics['celltype_mean_certainty'] = None
    else:
        log_message("Warning: Cell type predictions or uncertainty not found")
        metrics['celltype_mean_certainty'] = None
    
    # Check if lineage predictions and uncertainty exist
    if ('lineage_pred' in adata.obs and 'lineage_uncert' in adata.obs):
        
        # Determine which cells to use
        if consistent_mask is not None:
            # Use consistent cells
            mask_to_use = consistent_mask
            mask_name = "consistent"
        elif 'is_lineage_certain' in adata.obs:
            # Fallback to certain cells
            mask_to_use = adata.obs['is_lineage_certain'].astype(bool)
            mask_name = "certain"
            log_message("Fallback: Using 'is_lineage_certain' for lineage analysis")
        else:
            # Use all cells
            mask_to_use = np.ones(len(adata), dtype=bool)
            mask_name = "all"
            log_message("No consistency or certainty filter found, using all cells for lineage analysis")
        
        selected_cells = np.sum(mask_to_use)
        
        if selected_cells > 0:
            # Calculate overall mean certainty for selected cells only
            overall_lin_certainty = 1 - np.mean(adata.obs.loc[mask_to_use, 'lineage_uncert'])
            
            # Calculate per lineage mean certainty (for selected cells only)
            lin_mean_certainty = {}
            for lin in adata.obs.loc[mask_to_use, 'lineage_pred'].unique():
                lin_mask = (adata.obs['lineage_pred'] == lin) & mask_to_use
                if np.sum(lin_mask) > 0:
                    mean_cert = 1 - np.mean(adata.obs.loc[lin_mask, 'lineage_uncert'])
                    lin_mean_certainty[lin] = {
                        'mean_certainty': float(mean_cert),
                        'cell_count': int(np.sum(lin_mask))
                    }
            
            metrics['lineage_mean_certainty'] = {
                'overall': float(overall_lin_certainty),
                'per_label': lin_mean_certainty,
                'total_selected_cells': int(selected_cells),
                'total_cells': adata.n_obs,
                'selected_percentage': 100 * (selected_cells / adata.n_obs),
                'selection_criteria': mask_name
            }
        else:
            log_message(f"Warning: No {mask_name} cells found for lineage analysis")
            metrics['lineage_mean_certainty'] = None
    else:
        log_message("Warning: Lineage predictions or uncertainty not found")
        metrics['lineage_mean_certainty'] = None
    
    return metrics


def process_file(file_path, output_dir):
    """Process a single h5ad file to calculate mean certainty for consistent cells only"""
    file_name = os.path.basename(file_path).replace(".h5ad", "")
    log_message(f"Processing {file_name}")
    
    # Load annotated data
    try:
        adata = sc.read_h5ad(file_path)
        log_message(f"Loaded {file_name} with {adata.n_obs} cells")
    except Exception as e:
        log_message(f"Error loading {file_path}: {str(e)}")
        return None
    
    # Check if certainty flags exist, if not we might need to load the with_certainty file
    if ('reanno_pred_lineage' not in adata.obs or 'lineage_pred' not in adata.obs) and \
       ('is_celltype_certain' not in adata.obs or 'is_lineage_certain' not in adata.obs):
        certainty_file = file_path.replace("_scPoli_annotated.h5ad", "_with_certainty.h5ad")
        if os.path.exists(certainty_file):
            log_message(f"Loading certainty information from {certainty_file}")
            try:
                certainty_adata = sc.read_h5ad(certainty_file)
                # Copy certainty flags and consistency columns
                if 'is_celltype_certain' in certainty_adata.obs:
                    adata.obs['is_celltype_certain'] = certainty_adata.obs['is_celltype_certain'].values
                if 'is_lineage_certain' in certainty_adata.obs:
                    adata.obs['is_lineage_certain'] = certainty_adata.obs['is_lineage_certain'].values
                if 'reanno_pred_lineage' in certainty_adata.obs:
                    adata.obs['reanno_pred_lineage'] = certainty_adata.obs['reanno_pred_lineage'].values
                if 'lineage_pred' in certainty_adata.obs:
                    adata.obs['lineage_pred'] = certainty_adata.obs['lineage_pred'].values
            except Exception as e:
                log_message(f"Error loading certainty file: {str(e)}")
        else:
            log_message(f"Warning: Consistency/certainty information not found and no certainty file available")
    
    # Calculate mean certainty metrics for consistent cells only
    metrics = calculate_mean_certainty_consistent_cells(adata)
    
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
            ct_csv_path = os.path.join(model_dir, f"{file_name}_celltype_mean_certainty_consistent_cells.csv")
            ct_df.to_csv(ct_csv_path, index=False)
            log_message(f"Saved cell type mean certainty to {ct_csv_path}")
            
            # Visualize cell type mean certainty (for consistent cells only)
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(y='Cell_Type', x='Mean_Certainty', data=ct_df)
            
            # Add cell count
            for i, row in ct_df.iterrows():
                ax.text(row['Mean_Certainty'] + 0.01, i, f"n={row['Cell_Count']}", 
                       va='center', ha='left', fontsize=8)
            
            selection_criteria = metrics['celltype_mean_certainty']['selection_criteria']
            plt.title(f"Mean Certainty by Cell Type ({selection_criteria.title()} Cells Only)\n"
                      f"{file_name} - {metrics['celltype_mean_certainty']['total_selected_cells']} {selection_criteria} cells "
                      f"({metrics['celltype_mean_certainty']['selected_percentage']:.1f}%)")
            plt.xlabel('Mean Certainty (1 - Uncertainty)')
            plt.xlim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(model_dir, f"{file_name}_celltype_mean_certainty_consistent_cells.png"), 
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
            lin_csv_path = os.path.join(model_dir, f"{file_name}_lineage_mean_certainty_consistent_cells.csv")
            lin_df.to_csv(lin_csv_path, index=False)
            log_message(f"Saved lineage mean certainty to {lin_csv_path}")
            
            # Visualize lineage mean certainty (for consistent cells only)
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(y='Lineage', x='Mean_Certainty', data=lin_df)
            
            # Add cell count
            for i, row in lin_df.iterrows():
                ax.text(row['Mean_Certainty'] + 0.01, i, f"n={row['Cell_Count']}", 
                       va='center', ha='left', fontsize=8)
            
            selection_criteria = metrics['lineage_mean_certainty']['selection_criteria']
            plt.title(f"Mean Certainty by Lineage ({selection_criteria.title()} Cells Only)\n"
                      f"{file_name} - {metrics['lineage_mean_certainty']['total_selected_cells']} {selection_criteria} cells "
                      f"({metrics['lineage_mean_certainty']['selected_percentage']:.1f}%)")
            plt.xlabel('Mean Certainty (1 - Uncertainty)')
            plt.xlim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(model_dir, f"{file_name}_lineage_mean_certainty_consistent_cells.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    return metrics


def main():
    # Set paths
    data_dir = '/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output'  # Change this to your data directory
    output_dir = './mean_certainty_metrics_consistent_cells'
    
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
                    'Selected_Cells': metrics['celltype_mean_certainty']['total_selected_cells'],
                    'Total_Cells': metrics['celltype_mean_certainty']['total_cells'],
                    'Selected_Percentage': metrics['celltype_mean_certainty']['selected_percentage'],
                    'Selection_Criteria': metrics['celltype_mean_certainty']['selection_criteria']
                })
        
        if ct_overall_data:
            ct_overall_df = pd.DataFrame(ct_overall_data)
            ct_overall_df = ct_overall_df.sort_values('Mean_Certainty', ascending=False)
            ct_overall_path = os.path.join(output_dir, "celltype_overall_mean_certainty_consistent_cells_comparison.csv")
            ct_overall_df.to_csv(ct_overall_path, index=False)
            log_message(f"Saved cell type overall mean certainty comparison to {ct_overall_path}")
            
            # Plot comparison
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(x='Model', y='Mean_Certainty', data=ct_overall_df)
            
            # Add cell count and percentage
            for i, row in ct_overall_df.iterrows():
                ax.text(i, row['Mean_Certainty'] + 0.01, 
                        f"{row['Selected_Cells']} cells\n({row['Selected_Percentage']:.1f}%)\n{row['Selection_Criteria']}", 
                        ha='center', va='bottom', fontsize=8)
            
            plt.title("Overall Cell Type Mean Certainty Comparison (Consistent Cells Only)")
            plt.ylabel('Mean Certainty (1 - Uncertainty)')
            plt.xticks(rotation=90)
            plt.ylim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "celltype_overall_mean_certainty_consistent_cells_comparison.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
        
        # Lineage overall mean certainty comparison
        lin_overall_data = []
        for model, metrics in all_metrics.items():
            if metrics and metrics['lineage_mean_certainty'] is not None:
                lin_overall_data.append({
                    'Model': model,
                    'Mean_Certainty': metrics['lineage_mean_certainty']['overall'],
                    'Selected_Cells': metrics['lineage_mean_certainty']['total_selected_cells'],
                    'Total_Cells': metrics['lineage_mean_certainty']['total_cells'],
                    'Selected_Percentage': metrics['lineage_mean_certainty']['selected_percentage'],
                    'Selection_Criteria': metrics['lineage_mean_certainty']['selection_criteria']
                })
        
        if lin_overall_data:
            lin_overall_df = pd.DataFrame(lin_overall_data)
            lin_overall_df = lin_overall_df.sort_values('Mean_Certainty', ascending=False)
            lin_overall_path = os.path.join(output_dir, "lineage_overall_mean_certainty_consistent_cells_comparison.csv")
            lin_overall_df.to_csv(lin_overall_path, index=False)
            log_message(f"Saved lineage overall mean certainty comparison to {lin_overall_path}")
            
            # Plot comparison
            plt.figure(figsize=(12, 6))
            ax = sns.barplot(x='Model', y='Mean_Certainty', data=lin_overall_df)
            
            # Add cell count and percentage
            for i, row in lin_overall_df.iterrows():
                ax.text(i, row['Mean_Certainty'] + 0.01, 
                        f"{row['Selected_Cells']} cells\n({row['Selected_Percentage']:.1f}%)\n{row['Selection_Criteria']}", 
                        ha='center', va='bottom', fontsize=8)
            
            plt.title("Overall Lineage Mean Certainty Comparison (Consistent Cells Only)")
            plt.ylabel('Mean Certainty (1 - Uncertainty)')
            plt.xticks(rotation=90)
            plt.ylim(0, 1.1)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "lineage_overall_mean_certainty_consistent_cells_comparison.png"), 
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
            ct_comp_path = os.path.join(output_dir, "celltype_mean_certainty_by_label_consistent_cells_comparison.csv")
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
            lin_comp_path = os.path.join(output_dir, "lineage_mean_certainty_by_label_consistent_cells_comparison.csv")
            lin_comp_df.to_csv(lin_comp_path)
            log_message(f"Saved lineage mean certainty comparison by label to {lin_comp_path}")
    
    log_message("Mean certainty analysis (for consistent cells only) completed")


if __name__ == "__main__":
    main()

