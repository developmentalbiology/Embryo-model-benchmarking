
#!/usr/bin/env python
# coding: utf-8

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# ---------- FUNCTIONS ----------

def log_message(message, log_file='certainty_metrics.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")


def calculate_certainty(adata, uncertainty_threshold=0.2):
    """
    Calculate percentage of certain cells (uncertainty < threshold)
    
    Parameters
    ----------
    adata : AnnData
        Annotated dataset with uncertainty scores
    uncertainty_threshold : float
        Threshold below which cells are considered "certain"
        
    Returns
    -------
    dict
        Dictionary with certainty metrics for lineage and cell type
    """
    metrics = {}
    
    # Check if cell type predictions exist
    if 'final_anno_pred' in adata.obs and 'final_anno_uncert' in adata.obs:
        celltype_certain = adata.obs['final_anno_uncert'] < uncertainty_threshold
        metrics['celltype_certainty'] = {
            'certain_percentage': 100 * np.mean(celltype_certain),
            'certain_count': np.sum(celltype_certain),
            'total_count': len(celltype_certain)
        }
        # Add a column to identify certain cells
        adata.obs['is_celltype_certain'] = celltype_certain
    else:
        metrics['celltype_certainty'] = {
            'certain_percentage': np.nan,
            'certain_count': 0,
            'total_count': adata.n_obs
        }
        adata.obs['is_celltype_certain'] = False
    
    # Check if lineage predictions exist
    if 'final_lineage_pred' in adata.obs and 'final_lineage_uncert' in adata.obs:
        lineage_certain = adata.obs['final_lineage_uncert'] < uncertainty_threshold
        metrics['lineage_certainty'] = {
            'certain_percentage': 100 * np.mean(lineage_certain),
            'certain_count': np.sum(lineage_certain),
            'total_count': len(lineage_certain)
        }
        # Add a column to identify certain cells
        adata.obs['is_lineage_certain'] = lineage_certain
    else:
        metrics['lineage_certainty'] = {
            'certain_percentage': np.nan,
            'certain_count': 0,
            'total_count': adata.n_obs
        }
        adata.obs['is_lineage_certain'] = False
    
    # Per cell type certainty
    if 'final_anno_pred' in adata.obs and 'final_anno_uncert' in adata.obs:
        per_celltype = {}
        for ct in adata.obs['final_anno_pred'].unique():
            mask = adata.obs['final_anno_pred'] == ct
            ct_certainty = np.mean(adata.obs.loc[mask, 'final_anno_uncert'] < uncertainty_threshold)
            count = np.sum(mask)
            
            per_celltype[ct] = {
                'certain_percentage': 100 * ct_certainty,
                'total_count': count
            }
        metrics['per_celltype_certainty'] = per_celltype
    
    # Per lineage certainty
    if 'final_lineage_pred' in adata.obs and 'final_lineage_uncert' in adata.obs:
        per_lineage = {}
        for lin in adata.obs['final_lineage_pred'].unique():
            mask = adata.obs['final_lineage_pred'] == lin
            lin_certainty = np.mean(adata.obs.loc[mask, 'final_lineage_uncert'] < uncertainty_threshold)
            count = np.sum(mask)
            
            per_lineage[lin] = {
                'certain_percentage': 100 * lin_certainty,
                'total_count': count
            }
        metrics['per_lineage_certainty'] = per_lineage
    
    return metrics

def compare_certainty_distributions(all_files, output_dir, certainty_type='both'):
    """
    Compare the distribution of certainty values across multiple models.
    
    Parameters
    ----------
    all_files : list
        List of paths to h5ad files
    output_dir : str
        Directory for output files
    certainty_type : str
        'both', 'celltype', or 'lineage'
    """
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scanpy as sc
    
    # Create output directory
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    
    # Initialize data collection
    celltype_data = []
    lineage_data = []
    
    # Process each file
    for file_path in all_files:
        model_name = os.path.basename(file_path).replace("_scPoli_annotated.h5ad", "")
        
        # Load the data
        try:
            adata = sc.read_h5ad(file_path)
            
            # Collect cell type uncertainty data
            if 'final_anno_uncert' in adata.obs and (certainty_type == 'both' or certainty_type == 'celltype'):
                for idx, uncert in enumerate(adata.obs['final_anno_uncert']):
                    celltype_data.append({
                        'Model': model_name,
                        'Uncertainty': uncert,
                        'Cell ID': adata.obs.index[idx]
                    })
            
            # Collect lineage uncertainty data
            if 'final_lineage_uncert' in adata.obs and (certainty_type == 'both' or certainty_type == 'lineage'):
                for idx, uncert in enumerate(adata.obs['final_lineage_uncert']):
                    lineage_data.append({
                        'Model': model_name,
                        'Uncertainty': uncert,
                        'Cell ID': adata.obs.index[idx]
                    })
                
        except Exception as e:
            print(f"Error processing {model_name}: {str(e)}")
    
    # Create dataframes and sort models by mean uncertainty
    if celltype_data:
        ct_df = pd.DataFrame(celltype_data)
        
        # Calculate mean uncertainty per model for sorting
        model_means = ct_df.groupby('Model')['Uncertainty'].mean().reset_index()
        model_means = model_means.sort_values('Uncertainty')  # Sort by mean uncertainty
        model_order = model_means['Model'].tolist()  # Get ordered list of models
        
        # 1. Violin plot of cell type uncertainty distributions, ordered by mean
        plt.figure(figsize=(14, 8))
        sns.violinplot(x='Model', y='Uncertainty', data=ct_df, 
                      order=model_order, cut=0)
        plt.title("Distribution of Cell Type Uncertainty Across Models (Ordered by Mean)")
        plt.ylabel("Uncertainty Score")
        plt.xticks(rotation=90)
        plt.axhline(y=0.2, color='r', linestyle='--', alpha=0.7)  # Add threshold line
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "celltype_uncertainty_distribution.pdf"), bbox_inches='tight')
        plt.close()

        
        # 2. ECDF plot for cell type uncertainty
        plt.figure(figsize=(10, 6))
        # Use the same model order for consistent colors
        for model in model_order:
            model_data = ct_df[ct_df['Model'] == model]
            sns.ecdfplot(data=model_data, x='Uncertainty', label=model)
        plt.title("ECDF of Cell Type Uncertainty Across Models")
        plt.xlabel("Uncertainty Score")
        plt.ylabel("Cumulative Probability")
        plt.axvline(x=0.2, color='r', linestyle='--', alpha=0.7)  # Add threshold line
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "celltype_uncertainty_ecdf.pdf"), bbox_inches='tight')
        plt.close()
        
        # Save to CSV with mean uncertainty included
        model_stats = ct_df.groupby('Model')['Uncertainty'].agg(['mean', 'std', 'median', 'min', 'max']).reset_index()
        model_stats.to_csv(os.path.join(output_dir, "celltype_uncertainty_stats.csv"), index=False)
        ct_df.to_csv(os.path.join(output_dir, "celltype_uncertainty_distribution.csv"), index=False)
    
    # Process lineage data similarly
    if lineage_data:
        lin_df = pd.DataFrame(lineage_data)
        
        # Calculate mean uncertainty per model for sorting
        model_means = lin_df.groupby('Model')['Uncertainty'].mean().reset_index()
        model_means = model_means.sort_values('Uncertainty')  # Sort by mean uncertainty
        model_order = model_means['Model'].tolist()  # Get ordered list of models
        
        # 1. Violin plot, ordered by mean
        plt.figure(figsize=(14, 8))
        sns.violinplot(x='Model', y='Uncertainty', data=lin_df, 
                      order=model_order, cut=0)
        plt.title("Distribution of Lineage Uncertainty Across Models (Ordered by Mean)")
        plt.ylabel("Uncertainty Score")
        plt.xticks(rotation=90)
        plt.axhline(y=0.2, color='r', linestyle='--', alpha=0.7)  # Add threshold line
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "lineage_uncertainty_distribution.pdf"), bbox_inches='tight')
        plt.close()

        
        # 2. ECDF plot
        plt.figure(figsize=(10, 6))
        # Use the same model order for consistent colors
        for model in model_order:
            model_data = lin_df[lin_df['Model'] == model]
            sns.ecdfplot(data=model_data, x='Uncertainty', label=model)
        plt.title("ECDF of Lineage Uncertainty Across Models")
        plt.xlabel("Uncertainty Score")
        plt.ylabel("Cumulative Probability")
        plt.axvline(x=0.2, color='r', linestyle='--', alpha=0.7)  # Add threshold line
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "lineage_uncertainty_ecdf.pdf"), bbox_inches='tight')
        plt.close()
        
        # Save to CSV with mean uncertainty included
        model_stats = lin_df.groupby('Model')['Uncertainty'].agg(['mean', 'std', 'median', 'min', 'max']).reset_index()
        model


def process_file(file_path, output_dir, uncertainty_threshold=0.2):
    """Process a single h5ad file and add certainty flags"""
    file_name = os.path.basename(file_path).replace(".h5ad", "")
    log_message(f"Processing {file_name}")
    
    # Load annotated data
    adata = sc.read_h5ad(file_path)
    log_message(f"Loaded data with {adata.n_obs} cells")
    
    # Calculate certainty metrics
    metrics = calculate_certainty(adata, uncertainty_threshold)
    
    # Save metrics to JSON
    import json
    metrics_file = os.path.join(output_dir, f"{file_name}_certainty.json")
    
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
    
    # Process metrics dictionary for JSON
    json_metrics = {}
    for key, value in metrics.items():
        if isinstance(value, dict):
            json_metrics[key] = {}
            for k, v in value.items():
                if isinstance(v, dict):
                    json_metrics[key][k] = {k2: convert_for_json(v2) for k2, v2 in v.items()}
                else:
                    json_metrics[key][k] = convert_for_json(v)
        else:
            json_metrics[key] = convert_for_json(value)
    
    with open(metrics_file, 'w') as f:
        json.dump(json_metrics, f, indent=2)
    
    log_message(f"Saved certainty metrics to {metrics_file}")
    
    # Save updated AnnData by overwriting the original file
    adata.write_h5ad(file_path)
    log_message(f"Updated original file with certainty flags: {file_path}")
    
    return adata, metrics


def main():
    # Set paths
    data_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/output'  # Change this to your data directory
    output_dir = './certainty_metrics'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find annotated h5ad files
    h5ad_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith('_scPoli_annotated.h5ad')]
    log_message(f"Found {len(h5ad_files)} annotated files")
    
    # Process each file
    all_metrics = {}
    
    for file_path in h5ad_files:
        file_name = os.path.basename(file_path).replace("_scPoli_annotated.h5ad", "")
        try:
            adata, metrics = process_file(file_path, output_dir, uncertainty_threshold=0.2)
            all_metrics[file_name] = metrics
        except Exception as e:
            log_message(f"Error processing {file_name}: {str(e)}")
            import traceback
            log_message(traceback.format_exc())

    # Compare certainty distributions across models
    compare_certainty_distributions(h5ad_files, output_dir, certainty_type='both')
    
    # Create comparison of all models
    if len(all_metrics) > 1:
        # Create DataFrame for cell type certainty
        ct_data = []
        for model, metrics in all_metrics.items():
            if 'celltype_certainty' in metrics:
                ct_data.append({
                    'Model': model,
                    'Certain Percentage': metrics['celltype_certainty']['certain_percentage'],
                    'Certain Count': metrics['celltype_certainty']['certain_count'],
                    'Total Count': metrics['celltype_certainty']['total_count']
                })
        
        if ct_data:
            ct_df = pd.DataFrame(ct_data).sort_values('Certain Percentage', ascending=False)
            
            # Plot comparison
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Certain Percentage', data=ct_df)
            
            # Add count labels
            for i, row in ct_df.iterrows():
                ax.text(i, row['Certain Percentage'] + 1, 
                        f"{row['Certain Count']}/{row['Total Count']}", 
                        ha='center', va='bottom', rotation=90)
            
            plt.title("Cell Type Certainty Comparison Across Models")
            plt.ylabel('Percentage of Certain Cells (Uncertainty < 0.2)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "celltype_certainty_comparison.pdf"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            ct_df.to_csv(os.path.join(output_dir, "celltype_certainty_comparison.csv"), index=False)
        
        # Create DataFrame for lineage certainty
        lin_data = []
        for model, metrics in all_metrics.items():
            if 'lineage_certainty' in metrics:
                lin_data.append({
                    'Model': model,
                    'Certain Percentage': metrics['lineage_certainty']['certain_percentage'],
                    'Certain Count': metrics['lineage_certainty']['certain_count'],
                    'Total Count': metrics['lineage_certainty']['total_count']
                })
        
        if lin_data:
            lin_df = pd.DataFrame(lin_data).sort_values('Certain Percentage', ascending=False)
            
            # Plot comparison
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Certain Percentage', data=lin_df)
            
            # Add count labels
            for i, row in lin_df.iterrows():
                ax.text(i, row['Certain Percentage'] + 1, 
                        f"{row['Certain Count']}/{row['Total Count']}", 
                        ha='center', va='bottom', rotation=90)
            
            plt.title("Lineage Certainty Comparison Across Models")
            plt.ylabel('Percentage of Certain Cells (Uncertainty < 0.2)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "lineage_certainty_comparison.pdf"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            lin_df.to_csv(os.path.join(output_dir, "lineage_certainty_comparison.csv"), index=False)
    
    log_message("Completed certainty analysis")


if __name__ == "__main__":
    main()
