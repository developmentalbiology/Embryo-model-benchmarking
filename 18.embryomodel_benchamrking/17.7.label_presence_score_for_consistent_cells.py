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

def log_message(message, log_file='cell_type_presence.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")


def calculate_presence_scores(query_data, reference_data, 
                            query_label_key='reanno_pred', 
                            ref_label_key='reanno',
                            consistency_key=None):
    """
    Calculate cell type presence scores comparing query to reference
    
    Parameters
    ----------
    query_data : AnnData
        Query dataset with predicted labels
    reference_data : AnnData
        Reference atlas with ground truth labels
    query_label_key : str
        Column in query_data.obs containing predicted labels (e.g., 'reanno_pred')
    ref_label_key : str
        Column in reference_data.obs containing reference labels (e.g., 'reanno')
    consistency_key : str, optional
        Column in query_data.obs indicating which cells are consistent
        
    Returns
    -------
    dict
        Dictionary with presence metrics
    """
    log_message(f"Calculating cell type presence scores using {query_label_key} and {ref_label_key}...")
    
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
    
    # Get unique labels from datasets
    query_labels = sorted(query_filtered.obs[query_label_key].unique().tolist())
    ref_labels = sorted(reference_data.obs[ref_label_key].unique().tolist())
    
    # Find common labels
    common_labels = sorted(list(set(query_labels) & set(ref_labels)))
    log_message(f"Found {len(common_labels)} labels common to both query and reference (out of {len(ref_labels)} reference labels)")
    
    if not common_labels:
        log_message("No common labels found between query and reference")
        return {"error": "No common labels"}
    
    # Calculate presence score
    presence_score = len(common_labels) / len(ref_labels) * 100
    
    # Get count of cells for each common label
    query_counts = {}
    ref_counts = {}
    
    # Calculate counts in query dataset (only consistent cells)
    total_query_cells = len(query_filtered)
    for label in common_labels:
        count = np.sum(query_filtered.obs[query_label_key] == label)
        query_counts[label] = int(count)
    
    # Calculate counts in reference dataset
    total_ref_cells = len(reference_data)
    for label in ref_labels:
        count = np.sum(reference_data.obs[ref_label_key] == label)
        ref_counts[label] = int(count)
    
    # Create DataFrames for labels analysis
    presence_df = pd.DataFrame({
        'Label': ref_labels,
        'Present_in_Query': [label in common_labels for label in ref_labels],
        'Query_Count': [query_counts.get(label, 0) for label in ref_labels],
        'Reference_Count': [ref_counts[label] for label in ref_labels],
        'Reference_Percentage': [ref_counts[label] / total_ref_cells * 100 for label in ref_labels]
    })
    
    # Add a column to identify major cell types (>1% of reference)
    major_threshold = 1.0  # 1% of cells
    presence_df['Is_Major_Type'] = presence_df['Reference_Percentage'] >= major_threshold
    
    # Calculate presence score for major cell types
    major_labels = presence_df[presence_df['Is_Major_Type']]['Label'].tolist()
    major_present = presence_df[(presence_df['Is_Major_Type']) & (presence_df['Present_in_Query'])]['Label'].tolist()
    
    if major_labels:
        major_presence_score = len(major_present) / len(major_labels) * 100
    else:
        major_presence_score = np.nan
    
    log_message(f"Overall presence score: {presence_score:.1f}% ({len(common_labels)}/{len(ref_labels)} reference labels)")
    log_message(f"Major types presence score: {major_presence_score:.1f}% ({len(major_present)}/{len(major_labels)} major reference labels)")
    
    # Calculate total percentage of reference cell types that are captured
    captured_ref_percentage = sum(presence_df[presence_df['Present_in_Query']]['Reference_Percentage'])
    log_message(f"Captured {captured_ref_percentage:.1f}% of reference cells by type")
    
    # Calculate total cells in query that map to reference
    query_common_cells = sum(query_counts.values())
    query_common_pct = query_common_cells / total_query_cells * 100 if total_query_cells > 0 else 0
    
    # Store overall metrics
    overall_metrics = {
        'presence_score': float(presence_score),
        'major_presence_score': float(major_presence_score),
        'common_labels': len(common_labels),
        'total_reference_labels': len(ref_labels),
        'query_cells_in_common_labels': int(query_common_cells),
        'query_cells_in_common_labels_pct': float(query_common_pct),
        'captured_reference_percentage': float(captured_ref_percentage),
        'major_label_threshold_pct': float(major_threshold)
    }
    
    # Return all metrics
    return {
        'presence_df': presence_df,
        'overall_metrics': overall_metrics,
        'consistent_cells': int(np.sum(consistent_mask)),
        'total_cells': int(len(query_data))
    }


def calculate_lineage_presence_scores(query_data, reference_data,
                                    query_lineage_key='lineage_pred', 
                                    ref_lineage_key='lineage',
                                    consistency_key=None):
    """Calculate presence scores for lineages"""
    return calculate_presence_scores(
        query_data, reference_data,
        query_label_key=query_lineage_key,
        ref_label_key=ref_lineage_key,
        consistency_key=consistency_key
    )


def visualize_presence(metrics, file_name, output_dir, label_type="Cell Type"):
    """
    Create visualizations of presence scores
    
    Parameters
    ----------
    metrics : dict
        Dictionary with presence metrics
    file_name : str
        Base name for output files
    output_dir : str
        Directory for output files
    label_type : str
        Type of label (Cell Type or Lineage)
    """
    if "error" in metrics:
        log_message(f"Error in {label_type} presence metrics: {metrics['error']}")
        return
    
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    
    # Get presence dataframe
    df = metrics['presence_df']
    
    # 1. Horizontal bar plot of reference cells by type, colored by presence in query
    plt.figure(figsize=(12, 10))
    
    # Sort by reference percentage
    plot_df = df.sort_values('Reference_Percentage', ascending=True)
    
    # Color bars based on presence
    colors = ['#d55e00' if not present else '#0072b2' for present in plot_df['Present_in_Query']]
    
    # Create barplot
    ax = plt.barh(plot_df['Label'], plot_df['Reference_Percentage'], color=colors)
    
    # Add count labels
    for i, row in plot_df.iterrows():
        if row['Present_in_Query']:
            ha = 'left'
            offset = 0.2
            plt.text(row['Reference_Percentage'] + offset, i, 
                    f"Present: {row['Query_Count']} cells", 
                    va='center', ha=ha, fontsize=8)
        else:
            ha = 'right'
            offset = -0.2
            plt.text(row['Reference_Percentage'] + offset, i, 
                    f"Missing", 
                    va='center', ha=ha, fontsize=8, color='darkred')
    
    # Add a vertical line at major type threshold
    plt.axvline(x=metrics['overall_metrics']['major_label_threshold_pct'], 
               color='black', linestyle='--', alpha=0.5,
               label=f"Major type threshold ({metrics['overall_metrics']['major_label_threshold_pct']}%)")
    
    plt.title(f"{label_type} Presence Analysis\n{file_name}\nUsing {metrics['consistent_cells']} consistent cells")
    plt.xlabel(f'Reference {label_type} Percentage')
    plt.grid(axis='x', alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_presence.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Present vs Missing pie chart
    plt.figure(figsize=(10, 10))
    
    # All labels
    all_labels_values = [
        df[df['Present_in_Query']].shape[0],  # Present 
        df[~df['Present_in_Query']].shape[0]  # Missing
    ]
    
    plt.subplot(1, 2, 1)
    plt.pie(
        all_labels_values, 
        labels=['Present', 'Missing'],
        colors=['#0072b2', '#d55e00'],
        autopct='%1.1f%%',
        startangle=90,
        wedgeprops={'edgecolor': 'w', 'linewidth': 1}
    )
    plt.title(f"All {label_type}s\n({metrics['overall_metrics']['common_labels']}/{metrics['overall_metrics']['total_reference_labels']} labels)")
    
    # Major labels only
    major_labels = df[df['Is_Major_Type']]
    major_labels_values = [
        major_labels[major_labels['Present_in_Query']].shape[0],  # Present 
        major_labels[~major_labels['Present_in_Query']].shape[0]  # Missing
    ]
    
    plt.subplot(1, 2, 2)
    plt.pie(
        major_labels_values, 
        labels=['Present', 'Missing'],
        colors=['#0072b2', '#d55e00'],
        autopct='%1.1f%%',
        startangle=90,
        wedgeprops={'edgecolor': 'w', 'linewidth': 1}
    )
    plt.title(f"Major {label_type}s (>{metrics['overall_metrics']['major_label_threshold_pct']}%)\n({major_labels_values[0]}/{sum(major_labels_values)} labels)")
    
    plt.suptitle(f"{label_type} Presence in {file_name}")
    plt.tight_layout()
    
    plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_presence_pie.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Overall metrics in textbox
    plt.figure(figsize=(8, 6))
    plt.text(0.5, 0.5, 
             f"Overall {label_type} Presence Metrics\n"
             f"-----------------------------------\n"
             f"Presence Score: {metrics['overall_metrics']['presence_score']:.1f}% ({metrics['overall_metrics']['common_labels']} of {metrics['overall_metrics']['total_reference_labels']} labels)\n"
             f"Major Types Presence: {metrics['overall_metrics']['major_presence_score']:.1f}%\n"
             f"Captured Reference Percentage: {metrics['overall_metrics']['captured_reference_percentage']:.1f}%\n\n"
             f"Number of Consistent Cells: {metrics['consistent_cells']} of {metrics['total_cells']} ({metrics['consistent_cells']/metrics['total_cells']*100:.1f}%)",
             horizontalalignment='center',
             verticalalignment='center',
             fontsize=12,
             transform=plt.gca().transAxes,
             bbox=dict(boxstyle='round,pad=1', facecolor='azure', alpha=0.8)
            )
    
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_presence_metrics.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Missing major cell types
    missing_major = df[(df['Is_Major_Type']) & (~df['Present_in_Query'])].sort_values('Reference_Percentage', ascending=False)
    
    if not missing_major.empty:
        plt.figure(figsize=(12, 6))
        ax = plt.bar(missing_major['Label'], missing_major['Reference_Percentage'], color='#d55e00')
        
        plt.title(f"Missing Major {label_type}s in {file_name}")
        plt.ylabel(f'Reference {label_type} Percentage')
        plt.xticks(rotation=90)
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(os.path.join(fig_dir, f"{file_name}_{label_type.lower().replace(' ', '_')}_missing_major.png"), 
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
    
    metrics = {}
    
    # Cell type presence
    if 'reanno_pred' in query_adata.obs and 'reanno' in reference_adata.obs:
        log_message("Calculating cell type presence")
        
        # Check if we have consistency information
        consistency_key = None  # Let the function auto-detect consistency
        
        metrics['celltype_presence'] = calculate_presence_scores(
            query_adata, reference_adata,
            query_label_key='reanno_pred',
            ref_label_key='reanno',
            consistency_key=consistency_key
        )
        
        # Visualize presence metrics
        visualize_presence(
            metrics['celltype_presence'], 
            query_name, 
            output_dir, 
            label_type="Cell Type"
        )
        
        # Save cell type presence dataframe
        if 'presence_df' in metrics['celltype_presence']:
            ct_df_path = os.path.join(output_dir, f"{query_name}_celltype_presence.csv")
            metrics['celltype_presence']['presence_df'].to_csv(ct_df_path, index=False)
            log_message(f"Saved cell type presence dataframe to {ct_df_path}")
    else:
        if 'reanno_pred' not in query_adata.obs:
            log_message("Warning: 'reanno_pred' not found in query data")
        if 'reanno' not in reference_adata.obs:
            log_message("Warning: 'reanno' not found in reference data")
    
    # Lineage presence
    if 'lineage_pred' in query_adata.obs and 'lineage' in reference_adata.obs:
        log_message("Calculating lineage presence")
        
        # Check if we have consistency information
        consistency_key = None  # Let the function auto-detect consistency
        
        metrics['lineage_presence'] = calculate_lineage_presence_scores(
            query_adata, reference_adata,
            query_lineage_key='lineage_pred',
            ref_lineage_key='lineage',
            consistency_key=consistency_key
        )
        
        # Visualize presence metrics
        visualize_presence(
            metrics['lineage_presence'], 
            query_name, 
            output_dir, 
            label_type="Lineage"
        )
        
        # Save lineage presence dataframe
        if 'presence_df' in metrics['lineage_presence']:
            lin_df_path = os.path.join(output_dir, f"{query_name}_lineage_presence.csv")
            metrics['lineage_presence']['presence_df'].to_csv(lin_df_path, index=False)
            log_message(f"Saved lineage presence dataframe to {lin_df_path}")
    else:
        if 'lineage_pred' not in query_adata.obs:
            log_message("Warning: 'lineage_pred' not found in query data")
        if 'lineage' not in reference_adata.obs:
            log_message("Warning: 'lineage' not found in reference data")
    
    # Save metrics to JSON
    import json
    metrics_file = os.path.join(output_dir, f"{query_name}_presence.json")
    
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
    
    for key in ['celltype_presence', 'lineage_presence']:
        if key in metrics:
            json_metrics[key] = {}
            for k, v in metrics[key].items():
                if k != 'presence_df':  # Skip DataFrame
                    if isinstance(v, dict):
                        json_metrics[key][k] = {k2: convert_for_json(v2) for k2, v2 in v.items()}
                    else:
                        json_metrics[key][k] = convert_for_json(v)
    
    with open(metrics_file, 'w') as f:
        json.dump(json_metrics, f, indent=2)
    
    log_message(f"Saved presence metrics to {metrics_file}")
    
    return metrics


def main():
    # Set paths
    data_dir = '/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output'  # Change this to your data directory
    reference_file = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/data/human_clustering_20250724_v3.h5ad'  # Change this to your reference file
    output_dir = './cell_type_presence_metrics'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find h5ad files with consistency information
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
        # Create aggregate presence dataframes
        all_ct_presence = []
        all_lin_presence = []
        
        # Collect presence data
        for model, metrics in all_metrics.items():
            # Cell type presence
            if 'celltype_presence' in metrics and 'presence_df' in metrics['celltype_presence']:
                df = metrics['celltype_presence']['presence_df'].copy()
                df['Model'] = model
                all_ct_presence.append(df)
            
            # Lineage presence
            if 'lineage_presence' in metrics and 'presence_df' in metrics['lineage_presence']:
                df = metrics['lineage_presence']['presence_df'].copy()
                df['Model'] = model
                all_lin_presence.append(df)
        
        # Combine and save all cell type presence
        if all_ct_presence:
            combined_ct_df = pd.concat(all_ct_presence, ignore_index=True)
            combined_ct_path = os.path.join(output_dir, "all_models_celltype_presence.csv")
            combined_ct_df.to_csv(combined_ct_path, index=False)
            log_message(f"Saved combined cell type presence to {combined_ct_path}")
        
        # Combine and save all lineage presence
        if all_lin_presence:
            combined_lin_df = pd.concat(all_lin_presence, ignore_index=True)
            combined_lin_path = os.path.join(output_dir, "all_models_lineage_presence.csv")
            combined_lin_df.to_csv(combined_lin_path, index=False)
            log_message(f"Saved combined lineage presence to {combined_lin_path}")
        
        # Cell type presence comparison
        ct_data = []
        for model, metrics in all_metrics.items():
            if 'celltype_presence' in metrics and 'overall_metrics' in metrics['celltype_presence']:
                presence = metrics['celltype_presence']['overall_metrics']
                ct_data.append({
                    'Model': model,
                    'Presence Score (%)': presence.get('presence_score', np.nan),
                    'Major Types Presence (%)': presence.get('major_presence_score', np.nan),
                    'Captured Reference (%)': presence.get('captured_reference_percentage', np.nan),
                    'Common Labels': presence.get('common_labels', 0),
                    'Total Reference Labels': presence.get('total_reference_labels', 0),
                    'Consistent Cells': metrics['celltype_presence'].get('consistent_cells', 0)
                })
        
        if ct_data:
            ct_df = pd.DataFrame(ct_data).sort_values('Presence Score (%)', ascending=False)
            
            # Plot comparison of presence score
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Presence Score (%)', data=ct_df)
            
            # Add label counts
            for i, row in ct_df.iterrows():
                ax.text(i, row['Presence Score (%)'] + 2, 
                       f"{row['Common Labels']}/{row['Total Reference Labels']} labels", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Cell Type Presence Score Comparison")
            plt.ylabel('Presence Score (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "celltype_presence_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Plot comparison of major types presence
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Major Types Presence (%)', data=ct_df)
            
            # Add consistent cell counts
            for i, row in ct_df.iterrows():
                ax.text(i, row['Major Types Presence (%)'] + 2, 
                       f"{row['Consistent Cells']} cells", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Major Cell Type Presence Comparison")
            plt.ylabel('Major Types Presence (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "major_celltype_presence_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            ct_df.to_csv(os.path.join(output_dir, "celltype_presence_comparison.csv"), index=False)
        
        # Lineage presence comparison
        lin_data = []
        for model, metrics in all_metrics.items():
            if 'lineage_presence' in metrics and 'overall_metrics' in metrics['lineage_presence']:
                presence = metrics['lineage_presence']['overall_metrics']
                lin_data.append({
                    'Model': model,
                    'Presence Score (%)': presence.get('presence_score', np.nan),
                    'Major Types Presence (%)': presence.get('major_presence_score', np.nan),
                    'Captured Reference (%)': presence.get('captured_reference_percentage', np.nan),
                    'Common Labels': presence.get('common_labels', 0),
                    'Total Reference Labels': presence.get('total_reference_labels', 0),
                    'Consistent Cells': metrics['lineage_presence'].get('consistent_cells', 0)
                })
        
        if lin_data:
            lin_df = pd.DataFrame(lin_data).sort_values('Presence Score (%)', ascending=False)
            
            # Plot comparison of presence score
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Presence Score (%)', data=lin_df)
            
            # Add label counts
            for i, row in lin_df.iterrows():
                ax.text(i, row['Presence Score (%)'] + 2, 
                       f"{row['Common Labels']}/{row['Total Reference Labels']} labels", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Lineage Presence Score Comparison")
            plt.ylabel('Presence Score (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "lineage_presence_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Plot comparison of major types presence
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Major Types Presence (%)', data=lin_df)
            
            # Add consistent cell counts
            for i, row in lin_df.iterrows():
                ax.text(i, row['Major Types Presence (%)'] + 2, 
                       f"{row['Consistent Cells']} cells", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Major Lineage Presence Comparison")
            plt.ylabel('Major Types Presence (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "major_lineage_presence_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            lin_df.to_csv(os.path.join(output_dir, "lineage_presence_comparison.csv"), index=False)
    
    log_message("Analysis complete!")


if __name__ == "__main__":
    main()

