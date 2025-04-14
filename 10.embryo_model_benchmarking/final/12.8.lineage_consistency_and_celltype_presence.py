#!/usr/bin/env python
# coding: utf-8
import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
# Removed matplotlib_venn import
# from matplotlib_venn import venn2, venn2_circles

# ---------- FUNCTIONS ----------
def log_message(message, log_file='lineage_consistency.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

# Define the cell type to lineage mapping
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

def add_predicted_lineage_from_celltype(adata, cell_type_col='final_anno_pred'):
    """
    Add a new column that maps cell types to their lineages
    
    Parameters
    ----------
    adata : AnnData
        Dataset with cell type predictions
    cell_type_col : str
        Column containing cell type predictions
    
    Returns
    -------
    bool
        Whether the operation was successful
    """
    if cell_type_col not in adata.obs.columns:
        log_message(f"Error: {cell_type_col} column not found")
        return False
    
    # Create new column for predicted lineage from cell type
    new_col = f"{cell_type_col}_lineage"
    
    # Map cell types to lineages
    adata.obs[new_col] = adata.obs[cell_type_col].map(CELL_TYPE_TO_LINEAGE)
    
    # Count how many were successfully mapped
    mapped_count = adata.obs[new_col].notna().sum()
    total_count = adata.n_obs
    unmapped = total_count - mapped_count
    
    if unmapped > 0:
        log_message(f"Warning: {unmapped} cells could not be mapped to a lineage")
        # Get list of unmapped cell types
        unmapped_types = adata.obs.loc[adata.obs[new_col].isna(), cell_type_col].unique()
        log_message(f"Unmapped cell types: {unmapped_types}")
    
    log_message(f"Added {new_col} column mapping cell types to lineages ({mapped_count} cells mapped)")
    
    return True

def calculate_lineage_consistency(adata, annotated_lineage_col='final_anno_pred_lineage', 
                               predicted_lineage_col='final_lineage_pred',
                               certain_cell_type_col='is_celltype_certain',
                               certain_lineage_col='is_lineage_certain'):
    """
    Calculate consistency between lineage derived from cell type and predicted lineage
    
    Parameters
    ----------
    adata : AnnData
        Dataset with lineage columns
    annotated_lineage_col : str
        Column containing lineage derived from cell type annotation
    predicted_lineage_col : str
        Column containing directly predicted lineage
    certain_cell_type_col : str
        Column indicating certain cell type predictions
    certain_lineage_col : str
        Column indicating certain lineage predictions
    
    Returns
    -------
    dict
        Dictionary with consistency metrics
    """
    # Check if required columns exist
    required_cols = [annotated_lineage_col, predicted_lineage_col]
    for col in required_cols:
        if col not in adata.obs.columns:
            log_message(f"Error: {col} column not found")
            return {"error": f"{col} not found"}
    
    # Define certainty filters
    use_cell_type_certainty = certain_cell_type_col in adata.obs.columns
    use_lineage_certainty = certain_lineage_col in adata.obs.columns
    
    # Get cells with both certain cell type and lineage predictions
    if use_cell_type_certainty and use_lineage_certainty:
        certain_mask = adata.obs[certain_cell_type_col] & adata.obs[certain_lineage_col]
        log_message(f"Using {certain_mask.sum()}/{len(certain_mask)} cells with both certain cell type and lineage predictions")
    elif use_cell_type_certainty:
        certain_mask = adata.obs[certain_cell_type_col]
        log_message(f"Using {certain_mask.sum()}/{len(certain_mask)} cells with certain cell type predictions")
    elif use_lineage_certainty:
        certain_mask = adata.obs[certain_lineage_col]
        log_message(f"Using {certain_mask.sum()}/{len(certain_mask)} cells with certain lineage predictions")
    else:
        certain_mask = pd.Series(True, index=adata.obs.index)
        log_message(f"No certainty information found, using all {len(certain_mask)} cells")
    
    # Get only cells with valid lineage mapping from cell type
    valid_mapping_mask = adata.obs[annotated_lineage_col].notna()
    analysis_mask = certain_mask & valid_mapping_mask
    
    log_message(f"Analyzing lineage consistency for {analysis_mask.sum()}/{len(analysis_mask)} cells")
    
    if analysis_mask.sum() == 0:
        log_message("Error: No cells meet criteria for analysis")
        return {"error": "No cells meet criteria for analysis"}
    
    # Get lineage values for filtered cells
    annotated_lineages = adata.obs.loc[analysis_mask, annotated_lineage_col]
    predicted_lineages = adata.obs.loc[analysis_mask, predicted_lineage_col]
    
    # Calculate overall consistency
    consistency_mask = annotated_lineages == predicted_lineages
    overall_consistency = consistency_mask.mean() * 100
    log_message(f"Overall lineage consistency: {overall_consistency:.2f}%")
    
    # Calculate consistency by cell type derived lineage
    annotated_lineage_consistency = {}
    for lineage in annotated_lineages.unique():
        lineage_mask = annotated_lineages == lineage
        if lineage_mask.sum() > 0:
            lineage_consistency = (predicted_lineages.loc[lineage_mask] == lineage).mean() * 100
            annotated_lineage_consistency[lineage] = {
                'total_cells': int(lineage_mask.sum()),
                'consistent_cells': int((predicted_lineages.loc[lineage_mask] == lineage).sum()),
                'consistency_pct': float(lineage_consistency)
            }
            log_message(f"  {lineage}: {lineage_consistency:.2f}% consistent ({lineage_mask.sum()} cells)")
    
    # Calculate consistency by predicted lineage
    predicted_lineage_consistency = {}
    for lineage in predicted_lineages.unique():
        lineage_mask = predicted_lineages == lineage
        if lineage_mask.sum() > 0:
            lineage_consistency = (annotated_lineages.loc[lineage_mask] == lineage).mean() * 100
            predicted_lineage_consistency[lineage] = {
                'total_cells': int(lineage_mask.sum()),
                'consistent_cells': int((annotated_lineages.loc[lineage_mask] == lineage).sum()),
                'consistency_pct': float(lineage_consistency)
            }
            log_message(f"  {lineage} (predicted): {lineage_consistency:.2f}% consistent ({lineage_mask.sum()} cells)")
    
    # Create confusion matrix
    confusion_data = pd.crosstab(annotated_lineages, predicted_lineages, 
                                 rownames=['Annotated Lineage'], 
                                 colnames=['Predicted Lineage'])
    
    # Create consistency DataFrame
    consistency_df = pd.DataFrame({
        'Annotated_Lineage': annotated_lineages,
        'Predicted_Lineage': predicted_lineages,
        'Is_Consistent': consistency_mask
    })
    
    # If we have cell types, add them to the DataFrame
    if 'final_anno_pred' in adata.obs:
        consistency_df['Cell_Type'] = adata.obs.loc[analysis_mask, 'final_anno_pred']
        
        # Calculate consistency by cell type
        cell_type_consistency = {}
        for cell_type in consistency_df['Cell_Type'].unique():
            ct_mask = consistency_df['Cell_Type'] == cell_type
            if ct_mask.sum() > 0:
                ct_consistency = consistency_df.loc[ct_mask, 'Is_Consistent'].mean() * 100
                ct_lineage = CELL_TYPE_TO_LINEAGE.get(cell_type, "Unknown")
                cell_type_consistency[cell_type] = {
                    'total_cells': int(ct_mask.sum()),
                    'consistent_cells': int(consistency_df.loc[ct_mask, 'Is_Consistent'].sum()),
                    'consistency_pct': float(ct_consistency),
                    'expected_lineage': ct_lineage
                }
    
    # Return all metrics
    return {
        'overall_consistency_pct': float(overall_consistency),
        'consistent_cells': int(consistency_mask.sum()),
        'total_cells_analyzed': int(analysis_mask.sum()),
        'annotated_lineage_consistency': annotated_lineage_consistency,
        'predicted_lineage_consistency': predicted_lineage_consistency,
        'cell_type_consistency': cell_type_consistency if 'final_anno_pred' in adata.obs else {},
        'confusion_matrix': confusion_data,
        'consistency_df': consistency_df
    }

def calculate_lineage_presence(adata, reference_adata, 
                             lineage_col='final_anno_pred_lineage',
                             ref_lineage_col='final_lineage',
                             cell_type_col='final_anno_pred',
                             ref_cell_type_col='final_anno',
                             certainty_col='is_celltype_certain'):
    """
    Calculate presence of cell types within each lineage compared to reference
    
    Parameters
    ----------
    adata : AnnData
        Query dataset
    reference_adata : AnnData
        Reference atlas
    lineage_col : str
        Column in query with lineage labels
    ref_lineage_col : str
        Column in reference with lineage labels
    cell_type_col : str
        Column in query with cell type labels
    ref_cell_type_col : str
        Column in reference with cell type labels
    certainty_col : str
        Column in query indicating certain predictions
    
    Returns
    -------
    dict
        Dictionary with presence metrics
    """
    # Check if required columns exist
    required_query_cols = [lineage_col, cell_type_col]
    required_ref_cols = [ref_lineage_col, ref_cell_type_col]
    
    for col in required_query_cols:
        if col not in adata.obs.columns:
            log_message(f"Error: {col} column not found in query")
            return {"error": f"{col} not found in query"}
    
    for col in required_ref_cols:
        if col not in reference_adata.obs.columns:
            log_message(f"Error: {col} column not found in reference")
            return {"error": f"{col} not found in reference"}
    
    # Get certainty mask
    if certainty_col in adata.obs:
        certain_mask = adata.obs[certainty_col]
        log_message(f"Using {certain_mask.sum()}/{len(certain_mask)} certain cells")
    else:
        certain_mask = pd.Series(True, index=adata.obs.index)
        log_message(f"No certainty column found, using all {len(certain_mask)} cells")
    
    # Only analyze cells with valid lineage mapping
    valid_mapping_mask = adata.obs[lineage_col].notna()
    analysis_mask = certain_mask & valid_mapping_mask
    
    log_message(f"Analyzing lineage presence for {analysis_mask.sum()}/{len(analysis_mask)} cells")
    
    # Get list of lineages in query and reference
    query_lineages = adata.obs.loc[analysis_mask, lineage_col].unique()
    ref_lineages = reference_adata.obs[ref_lineage_col].unique()
    
    # Find common lineages
    common_lineages = set(query_lineages) & set(ref_lineages)
    log_message(f"Found {len(common_lineages)} lineages common to both query and reference")
    
    # Overall lineage presence
    lineage_presence_score = len(common_lineages) / len(ref_lineages) * 100
    log_message(f"Overall lineage presence: {lineage_presence_score:.2f}% ({len(common_lineages)}/{len(ref_lineages)})")
    
    # Calculate presence of cell types within each lineage
    lineage_cell_type_presence = {}
    cell_type_presence_df_rows = []
    
    for lineage in sorted(common_lineages):
        # Get cell types in reference for this lineage
        ref_lineage_mask = reference_adata.obs[ref_lineage_col] == lineage
        ref_cell_types = reference_adata.obs.loc[ref_lineage_mask, ref_cell_type_col].unique()
        
        # Get cell types in query for this lineage
        query_lineage_mask = (adata.obs[lineage_col] == lineage) & analysis_mask
        query_cell_types = adata.obs.loc[query_lineage_mask, cell_type_col].unique()
        
        # Find common cell types
        common_cell_types = set(query_cell_types) & set(ref_cell_types)
        
        # Calculate presence score for this lineage
        if len(ref_cell_types) > 0:
            lineage_ct_presence = len(common_cell_types) / len(ref_cell_types) * 100
        else:
            lineage_ct_presence = np.nan
        
        # Store metrics for this lineage
        lineage_cell_type_presence[lineage] = {
            'query_cell_types': list(query_cell_types),
            'reference_cell_types': list(ref_cell_types),
            'common_cell_types': list(common_cell_types),
            'query_only_cell_types': list(set(query_cell_types) - set(ref_cell_types)),
            'reference_only_cell_types': list(set(ref_cell_types) - set(query_cell_types)),
            'query_cell_count': int(query_lineage_mask.sum()),
            'reference_cell_count': int(ref_lineage_mask.sum()),
            'presence_score': float(lineage_ct_presence)
        }
        
        log_message(f"  {lineage}: {lineage_ct_presence:.2f}% cell type presence ({len(common_cell_types)}/{len(ref_cell_types)})")
        
        # Add to dataframe rows for this lineage
        for ct in sorted(set(ref_cell_types) | set(query_cell_types)):
            cell_type_presence_df_rows.append({
                'Lineage': lineage,
                'Cell_Type': ct,
                'In_Query': ct in query_cell_types,
                'In_Reference': ct in ref_cell_types,
                'Query_Cell_Count': int(query_lineage_mask.sum()),
                'Reference_Cell_Count': int(ref_lineage_mask.sum()),
                'Query_Type_Count': len(query_cell_types),
                'Reference_Type_Count': len(ref_cell_types),
                'Common_Type_Count': len(common_cell_types),
                'Presence_Score': float(lineage_ct_presence)
            })
    
    # Create presence DataFrame
    cell_type_presence_df = pd.DataFrame(cell_type_presence_df_rows)
    
    # Calculate overall cell type presence
    all_ref_cell_types = set(reference_adata.obs[ref_cell_type_col].unique())
    all_query_cell_types = set(adata.obs.loc[analysis_mask, cell_type_col].unique())
    all_common_cell_types = all_ref_cell_types & all_query_cell_types
    
    overall_ct_presence = len(all_common_cell_types) / len(all_ref_cell_types) * 100
    log_message(f"Overall cell type presence: {overall_ct_presence:.2f}% ({len(all_common_cell_types)}/{len(all_ref_cell_types)})")
    
    # Calculate weighted lineage presence (weighted by reference cell count)
    ref_lineage_counts = reference_adata.obs[ref_lineage_col].value_counts()
    weighted_scores = []
    weighted_counts = []
    
    for lineage in common_lineages:
        if lineage in ref_lineage_counts:
            presence = lineage_cell_type_presence[lineage]['presence_score']
            weight = ref_lineage_counts[lineage]
            weighted_scores.append(presence * weight)
            weighted_counts.append(weight)
    
    if weighted_counts:
        weighted_lineage_presence = sum(weighted_scores) / sum(weighted_counts)
        log_message(f"Weighted lineage presence: {weighted_lineage_presence:.2f}%")
    else:
        weighted_lineage_presence = np.nan
        log_message("Warning: Could not calculate weighted lineage presence")
    
    # Return all metrics
    return {
        'overall_lineage_presence_pct': float(lineage_presence_score),
        'overall_cell_type_presence_pct': float(overall_ct_presence),
        'weighted_lineage_presence_pct': float(weighted_lineage_presence),
        'common_lineages': list(common_lineages),
        'query_only_lineages': list(set(query_lineages) - set(ref_lineages)),
        'reference_only_lineages': list(set(ref_lineages) - set(query_lineages)),
        'lineage_cell_type_presence': lineage_cell_type_presence,
        'total_reference_lineages': len(ref_lineages),
        'total_query_lineages': len(query_lineages),
        'cell_type_presence_df': cell_type_presence_df
    }

def visualize_lineage_consistency(metrics, file_name, output_dir):
    """
    Visualize lineage consistency results
    
    Parameters
    ----------
    metrics : dict
        Dictionary with consistency metrics
    file_name : str
        Base name for output files
    output_dir : str
        Directory for output files
    """
    if "error" in metrics:
        log_message(f"Error in lineage consistency metrics: {metrics['error']}")
        return
    
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    
    # 1. Overall consistency bar
    plt.figure(figsize=(6, 4))
    plt.bar(['Lineage Consistency'], [metrics['overall_consistency_pct']], color='#0072b2')
    plt.title(f"Overall Lineage Consistency\n{file_name}")
    plt.ylabel('Consistency (%)')
    plt.ylim(0, 105)
    
    # Add count label
    plt.text(0, metrics['overall_consistency_pct'] + 2, 
             f"{metrics['consistent_cells']}/{metrics['total_cells_analyzed']} cells", 
             ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{file_name}_overall_consistency.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Consistency by annotated lineage (derived from cell types)
    if metrics['annotated_lineage_consistency']:
        ann_data = []
        for lineage, data in metrics['annotated_lineage_consistency'].items():
            ann_data.append({
                'Lineage': lineage,
                'Consistency': data['consistency_pct'],
                'Cell Count': data['total_cells']
            })
        
        ann_df = pd.DataFrame(ann_data).sort_values('Consistency', ascending=False)
        
        plt.figure(figsize=(12, 8))
        ax = sns.barplot(data=ann_df, x='Lineage', y='Consistency')
        
        # Add count labels
        for i, row in ann_df.iterrows():
            plt.text(i, row['Consistency'] + 2, f"{row['Cell Count']} cells", 
                    ha='center', va='bottom')
        
        plt.title(f"Lineage Consistency by Annotated Lineage\n{file_name}")
        plt.ylabel('Consistency (%)')
        plt.ylim(0, 105)
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{file_name}_annotated_lineage_consistency.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Confusion matrix heatmap
    if 'confusion_matrix' in metrics:
        plt.figure(figsize=(14, 12))
        conf_matrix = metrics['confusion_matrix']
        
        # Normalize by row (annotated lineage)
        conf_matrix_norm = conf_matrix.div(conf_matrix.sum(axis=1), axis=0) * 100
        
        sns.heatmap(conf_matrix_norm, annot=conf_matrix, fmt='d', cmap='Blues',
                   cbar_kws={'label': 'Percentage of Annotated Lineage'})
        
        plt.title(f"Lineage Confusion Matrix\n{file_name}")
        plt.ylabel('Lineage from Cell Type Annotation')
        plt.xlabel('Predicted Lineage')
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{file_name}_lineage_confusion_matrix.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 4. Top mismatches by cell type
    if 'cell_type_consistency' in metrics and metrics['cell_type_consistency']:
        ct_data = []
        for cell_type, data in metrics['cell_type_consistency'].items():
            if data['total_cells'] >= 10:  # Only include cell types with enough cells
                ct_data.append({
                    'Cell Type': cell_type,
                    'Consistency': data['consistency_pct'],
                    'Cell Count': data['total_cells'],
                    'Consistent Count': data['consistent_cells'],
                    'Expected Lineage': data['expected_lineage']
                })
        
        if ct_data:
            ct_df = pd.DataFrame(ct_data).sort_values('Consistency')
            
            # Show top 20 most inconsistent cell types
            if len(ct_df) > 20:
                ct_df = ct_df.head(20)
            
            plt.figure(figsize=(12, 10))
            ax = sns.barplot(data=ct_df, y='Cell Type', x='Consistency')
            
            # Add count labels
            for i, row in ct_df.iterrows():
                plt.text(row['Consistency'] + 2, i, 
                        f"{row['Consistent Count']}/{row['Cell Count']}", 
                        va='center', ha='left', fontsize=8)
            
            plt.title(f"Cell Types with Lowest Lineage Consistency\n{file_name}")
            plt.xlabel('Lineage Consistency (%)')
            plt.xlim(0, 105)
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(fig_dir, f"{file_name}_cell_type_consistency.png"), dpi=300, bbox_inches='tight')
            plt.close()

def visualize_lineage_presence(metrics, file_name, output_dir):
    """
    Visualize lineage presence results
    
    Parameters
    ----------
    metrics : dict
        Dictionary with presence metrics
    file_name : str
        Base name for output files
    output_dir : str
        Directory for output files
    """
    if "error" in metrics:
        log_message(f"Error in lineage presence metrics: {metrics['error']}")
        return
    
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    
    # 1. Overall presence scores
    plt.figure(figsize=(10, 6))
    
    scores = [
        metrics['overall_lineage_presence_pct'],
        metrics['overall_cell_type_presence_pct'],
        metrics['weighted_lineage_presence_pct']
    ]
    
    labels = [
        f"Lineage Presence\n{len(metrics['common_lineages'])}/{metrics['total_reference_lineages']} lineages",
        f"Cell Type Presence",
        f"Weighted Lineage Presence"
    ]
    
    plt.bar(labels, scores, color=['#1f77b4', '#ff7f0e', '#2ca02c'])
    
    plt.title(f"Overall Presence Scores\n{file_name}")
    plt.ylabel('Presence Score (%)')
    plt.ylim(0, 105)
    plt.xticks(rotation=0)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{file_name}_overall_presence.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Replace Venn diagram with a simple bar chart
    try:
        plt.figure(figsize=(8, 8))
        
        # Calculate categories and counts
        common_count = len(metrics['common_lineages'])
        query_only_count = len(metrics['query_only_lineages'])
        ref_only_count = len(metrics['reference_only_lineages'])
        
        categories = ['Common', 'Query Only', 'Reference Only']
        counts = [common_count, query_only_count, ref_only_count]
        colors = ['#2ca02c', '#1f77b4', '#d62728']
        
        plt.bar(categories, counts, color=colors)
        
        # Add count labels
        for i, count in enumerate(counts):
            plt.text(i, count + 0.5, str(count), ha='center', va='bottom')
        
        plt.title(f"Lineage Distribution Between Query and Reference\n{file_name}")
        plt.ylabel('Number of Lineages')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{file_name}_lineage_overlap.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save overlap details to CSV
        overlap_details = pd.DataFrame({
            'Common_Lineages': pd.Series(list(metrics['common_lineages'])),
            'Query_Only_Lineages': pd.Series(list(metrics['query_only_lineages'])),
            'Reference_Only_Lineages': pd.Series(list(metrics['reference_only_lineages']))
        })
        
        overlap_details.to_csv(os.path.join(output_dir, f"{file_name}_lineage_overlap_details.csv"), index=False)
        
    except Exception as e:
        log_message(f"Error creating lineage overlap visualization: {str(e)}")
    
    # 3. Cell type presence by lineage
    if metrics['lineage_cell_type_presence']:
        lineage_data = []
        for lineage, data in metrics['lineage_cell_type_presence'].items():
            lineage_data.append({
                'Lineage': lineage,
                'Presence Score': data['presence_score'],
                'Query Count': data['query_cell_count'],
                'Reference Count': data['reference_cell_count'],
                'Query Types': len(data['query_cell_types']),
                'Reference Types': len(data['reference_cell_types']),
                'Common Types': len(data['common_cell_types'])
            })
        
        lineage_df = pd.DataFrame(lineage_data).sort_values('Presence Score', ascending=False)
        
        plt.figure(figsize=(12, 8))
        ax = sns.barplot(data=lineage_df, x='Lineage', y='Presence Score')
        
        # Add count labels
        for i, row in lineage_df.iterrows():
            plt.text(i, row['Presence Score'] + 2, 
                    f"{row['Common Types']}/{row['Reference Types']} types", 
                    ha='center', va='bottom', rotation=90, fontsize=9)
        
        plt.title(f"Cell Type Presence Score by Lineage\n{file_name}")
        plt.ylabel('Cell Type Presence Score (%)')
        plt.ylim(0, 105)
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{file_name}_cell_type_presence_by_lineage.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        # 4. Detailed cell type presence for each lineage
        for lineage, data in metrics['lineage_cell_type_presence'].items():
            if len(data['reference_cell_types']) > 0:
                plt.figure(figsize=(12, max(6, len(data['reference_cell_types']) * 0.4)))
                
                # Prepare data for plotting
                cell_types = sorted(set(data['reference_cell_types']) | set(data['query_cell_types']))
                in_ref = [ct in data['reference_cell_types'] for ct in cell_types]
                in_query = [ct in data['query_cell_types'] for ct in cell_types]
                
                # Create dataframe for plotting
                cell_type_df = pd.DataFrame({
                    'Cell Type': cell_types,
                    'In Reference': in_ref,
                    'In Query': in_query,
                    'Status': ['Both' if r and q else 'Reference Only' if r else 'Query Only' for r, q in zip(in_ref, in_query)]
                })
                
                # Sort by status
                cell_type_df['Status_Order'] = cell_type_df['Status'].map({'Both': 0, 'Reference Only': 1, 'Query Only': 2})
                cell_type_df = cell_type_df.sort_values(['Status_Order', 'Cell Type']).drop('Status_Order', axis=1)
                
                # Create color mapping for status
                colors = {'Both': '#2ca02c', 'Reference Only': '#d62728', 'Query Only': '#1f77b4'}
                
                # Create horizontal bar chart
                ax = sns.barplot(
                    y='Cell Type', 
                    x=[1] * len(cell_type_df),  # All bars same length
                    hue='Status',
                    palette=colors,
                    data=cell_type_df,
                    dodge=False
                )
                
                plt.title(f"Cell Type Presence in {lineage} Lineage\n{file_name}")
                plt.xlabel('')
                plt.xlim(0, 1.2)  # Make room for legend
                plt.legend(title='Status')
                ax.get_xaxis().set_visible(False)  # Hide x-axis
                plt.tight_layout()
                plt.savefig(os.path.join(fig_dir, f"{file_name}_{lineage}_cell_types.png"), dpi=300, bbox_inches='tight')
                plt.close()

def process_file(query_file, reference_file, output_dir):
    """Process a single query file against the reference"""
    query_name = os.path.basename(query_file).replace(".h5ad", "")
    log_message(f"Processing {query_name}")
    
    # Load query data
    query_adata = sc.read_h5ad(query_file)
    log_message(f"Loaded query data with {query_adata.n_obs} cells and {query_adata.n_vars} genes")
    
    # Add predicted lineage column based on cell type
    success = add_predicted_lineage_from_celltype(query_adata, cell_type_col='final_anno_pred')
    
    if not success:
        log_message(f"Error: Could not add lineage from cell type for {query_name}")
        return {"error": "Could not add lineage from cell type"}
    
    # Load reference data
    reference_adata = sc.read_h5ad(reference_file)
    log_message(f"Loaded reference data with {reference_adata.n_obs} cells and {reference_adata.n_vars} genes")
    
    metrics = {}
    
    # Calculate lineage consistency
    if 'final_anno_pred_lineage' in query_adata.obs and 'final_lineage_pred' in query_adata.obs:
        log_message("Calculating lineage consistency")
        
        metrics['lineage_consistency'] = calculate_lineage_consistency(
            query_adata,
            annotated_lineage_col='final_anno_pred_lineage',
            predicted_lineage_col='final_lineage_pred',
            certain_cell_type_col='is_celltype_certain' if 'is_celltype_certain' in query_adata.obs else None,
            certain_lineage_col='is_lineage_certain' if 'is_lineage_certain' in query_adata.obs else None
        )
        
        # Visualize consistency metrics
        visualize_lineage_consistency(
            metrics['lineage_consistency'], 
            query_name, 
            output_dir
        )
        
        # Save consistency dataframe
        if 'consistency_df' in metrics['lineage_consistency']:
            consistency_df_path = os.path.join(output_dir, f"{query_name}_lineage_consistency.csv")
            metrics['lineage_consistency']['consistency_df'].to_csv(consistency_df_path, index=False)
            log_message(f"Saved lineage consistency dataframe to {consistency_df_path}")
            
            # Save confusion matrix
            confusion_path = os.path.join(output_dir, f"{query_name}_lineage_confusion_matrix.csv")
            metrics['lineage_consistency']['confusion_matrix'].to_csv(confusion_path)
            log_message(f"Saved lineage confusion matrix to {confusion_path}")
    else:
        if 'final_anno_pred_lineage' not in query_adata.obs:
            log_message("Warning: 'final_anno_pred_lineage' not found in query data")
        if 'final_lineage_pred' not in query_adata.obs:
            log_message("Warning: 'final_lineage_pred' not found in query data")
    
    # Check available columns in reference data
    ref_cols = reference_adata.obs.columns.tolist()
    log_message(f"Available columns in reference data: {ref_cols}")
    
    # Determine reference column names
    ref_lineage_col = 'final_lineage' if 'final_lineage' in ref_cols else 'lineage'
    ref_celltype_col = 'final_anno' if 'final_anno' in ref_cols else 'cell_type'
    
    log_message(f"Using reference lineage column: {ref_lineage_col}")
    log_message(f"Using reference cell type column: {ref_celltype_col}")
    
    # Calculate lineage presence
    if 'final_anno_pred_lineage' in query_adata.obs and ref_lineage_col in reference_adata.obs:
        log_message("Calculating lineage presence")
        
        metrics['lineage_presence'] = calculate_lineage_presence(
            query_adata, 
            reference_adata,
            lineage_col='final_anno_pred_lineage',
            ref_lineage_col=ref_lineage_col,
            cell_type_col='final_anno_pred',
            ref_cell_type_col=ref_celltype_col,
            certainty_col='is_celltype_certain' if 'is_celltype_certain' in query_adata.obs else None
        )
        
        # Visualize presence metrics
        visualize_lineage_presence(
            metrics['lineage_presence'], 
            query_name, 
            output_dir
        )
        
        # Save presence dataframe
        if 'cell_type_presence_df' in metrics['lineage_presence']:
            presence_df_path = os.path.join(output_dir, f"{query_name}_lineage_cell_type_presence.csv")
            metrics['lineage_presence']['cell_type_presence_df'].to_csv(presence_df_path, index=False)
            log_message(f"Saved lineage cell type presence dataframe to {presence_df_path}")
    else:
        if 'final_anno_pred_lineage' not in query_adata.obs:
            log_message("Warning: 'final_anno_pred_lineage' not found in query data")
        if ref_lineage_col not in reference_adata.obs:
            log_message(f"Warning: '{ref_lineage_col}' not found in reference data")
    
    # Save metrics to JSON
    import json
    metrics_file = os.path.join(output_dir, f"{query_name}_lineage_analysis.json")
    
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
        elif isinstance(obj, set):
            return list(obj)
        else:
            return obj
    
    # Process metrics dictionary for JSON (exclude DataFrames)
    json_metrics = {}
    
    for key in ['lineage_consistency', 'lineage_presence']:
        if key in metrics:
            json_metrics[key] = {}
            for k, v in metrics[key].items():
                if k not in ['consistency_df', 'confusion_matrix', 'cell_type_presence_df']:  # Skip DataFrames
                    if isinstance(v, dict):
                        json_metrics[key][k] = {}
                        for k2, v2 in v.items():
                            if isinstance(v2, dict):
                                json_metrics[key][k][k2] = {k3: convert_for_json(v3) for k3, v3 in v2.items()}
                            else:
                                json_metrics[key][k][k2] = convert_for_json(v2)
                    else:
                        json_metrics[key][k] = convert_for_json(v)
    
    try:
        with open(metrics_file, 'w') as f:
            json.dump(json_metrics, f, indent=2)
        log_message(f"Saved lineage analysis metrics to {metrics_file}")
    except Exception as e:
        log_message(f"Error saving metrics to JSON: {str(e)}")
    
    # Save Updated AnnData with predicted lineage from cell type
    query_adata.write_h5ad(query_file)
    log_message(f"Saved updated AnnData with final_anno_pred_lineage column to {query_file}")
    
    return metrics

def main():
    # Set paths
    data_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/output'  # Change this to your data directory
    reference_file = '/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/human_reanno_20250108.h5ad'  # Change this to your reference file
    output_dir = './lineage_consistency_metrics'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find h5ad files with annotations
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
        # Lineage consistency comparison
        consistency_data = []
        
        for model, metrics in all_metrics.items():
            if 'lineage_consistency' in metrics and 'error' not in metrics['lineage_consistency']:
                consistency = metrics['lineage_consistency']
                consistency_data.append({
                    'Model': model,
                    'Consistency (%)': consistency.get('overall_consistency_pct', np.nan),
                    'Consistent Cells': consistency.get('consistent_cells', 0),
                    'Total Cells': consistency.get('total_cells_analyzed', 0)
                })
        
        if consistency_data:
            consistency_df = pd.DataFrame(consistency_data).sort_values('Consistency (%)', ascending=False)
            
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Consistency (%)', data=consistency_df)
            
            # Add count labels
            for i, row in consistency_df.iterrows():
                ax.text(i, row['Consistency (%)'] + 2, 
                       f"{row['Consistent Cells']}/{row['Total Cells']} cells", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Lineage Consistency Comparison Across Models")
            plt.ylabel('Consistency (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "lineage_consistency_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            consistency_df.to_csv(os.path.join(output_dir, "lineage_consistency_comparison.csv"), index=False)
        
        # Lineage presence comparison
        presence_data = []
        
        for model, metrics in all_metrics.items():
            if 'lineage_presence' in metrics and 'error' not in metrics['lineage_presence']:
                presence = metrics['lineage_presence']
                presence_data.append({
                    'Model': model,
                    'Lineage Presence (%)': presence.get('overall_lineage_presence_pct', np.nan),
                    'Cell Type Presence (%)': presence.get('overall_cell_type_presence_pct', np.nan),
                    'Weighted Lineage Presence (%)': presence.get('weighted_lineage_presence_pct', np.nan),
                    'Common Lineages': len(presence.get('common_lineages', [])),
                    'Total Reference Lineages': presence.get('total_reference_lineages', 0)
                })
        
        if presence_data:
            presence_df = pd.DataFrame(presence_data).sort_values('Weighted Lineage Presence (%)', ascending=False)
            
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='Model', y='Weighted Lineage Presence (%)', data=presence_df)
            
            # Add lineage counts
            for i, row in presence_df.iterrows():
                ax.text(i, row['Weighted Lineage Presence (%)'] + 2, 
                       f"{row['Common Lineages']}/{row['Total Reference Lineages']} lineages", 
                       ha='center', va='bottom', rotation=90, fontsize=8)
            
            plt.title("Weighted Lineage Presence Comparison Across Models")
            plt.ylabel('Weighted Lineage Presence (%)')
            plt.xticks(rotation=90)
            plt.ylim(0, 105)
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "weighted_lineage_presence_comparison.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save to CSV
            presence_df.to_csv(os.path.join(output_dir, "lineage_presence_comparison.csv"), index=False)
    
    log_message("Completed lineage consistency and presence analysis")

if __name__ == "__main__":
    main()