#!/usr/bin/env python
# coding: utf-8

"""
Cell Type Composition Similarity Analysis (Fixed Version)

This script analyzes cell type composition similarity between embryo model datasets
and reference atlas stages. It calculates correlations between dataset-stage combinations
based on their cell type composition vectors.

Fixed issues:
- Ensures all datasets (reference + embryo models) are included in correlation analysis
- Handles missing or incorrectly formatted orig.ident columns by using filename base names
- Properly merges data from different sources with different column structures
"""

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import warnings
import re

def log_message(message, log_file='composition_analysis.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

def is_valid_identifier(value):
    """Check if a value is a valid identifier (contains English letters)"""
    if pd.isna(value):
        return False
    str_value = str(value)
    # Check if contains at least one English letter
    return bool(re.search(r'[a-zA-Z]', str_value))

def calculate_cell_type_composition(adata, dataset_name, file_basename, cell_type_col='reanno', 
                                   stage_col='stage', orig_ident_col='orig.ident',
                                   use_consistent_cells=True):
    """
    Calculate cell type composition percentages for a dataset.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object
    dataset_name : str
        Name of the dataset
    file_basename : str
        Base name of the file (without extension)
    cell_type_col : str
        Column name for cell type annotations
    stage_col : str
        Column name for stage information (if available)
    orig_ident_col : str
        Column name for original identity/batch information (if available)
    use_consistent_cells : bool
        Whether to use only consistent cells (reanno_pred_lineage == lineage_pred)
    
    Returns:
    --------
    pd.DataFrame: DataFrame with cell type composition percentages
    """
    
    log_message(f"Calculating cell type composition for {dataset_name}")
    
    # Filter to consistent cells if requested
    if use_consistent_cells:
        if 'reanno_pred_lineage' in adata.obs and 'lineage_pred' in adata.obs:
            # Convert categorical to string to avoid comparison issues
            reanno_pred_lineage_str = adata.obs['reanno_pred_lineage'].astype(str)
            lineage_pred_str = adata.obs['lineage_pred'].astype(str)
            consistent_mask = (reanno_pred_lineage_str == lineage_pred_str)
            adata_filtered = adata[consistent_mask].copy()
            log_message(f"  Using {np.sum(consistent_mask)}/{len(adata)} consistent cells")
        else:
            log_message(f"  Consistency columns not found, using all {len(adata)} cells")
            adata_filtered = adata.copy()
    else:
        adata_filtered = adata.copy()
        log_message(f"  Using all {len(adata_filtered)} cells")
    
    # Check if required cell type column exists
    if cell_type_col not in adata_filtered.obs.columns:
        log_message(f"Warning: {cell_type_col} column not found in {dataset_name}")
        return pd.DataFrame()
    
    # Create metadata dataframe
    meta = adata_filtered.obs.copy()
    
    # Handle stage column
    if stage_col not in meta.columns:
        # If no stage column, create a default stage
        meta[stage_col] = 'embryo_model'
        log_message(f"  No {stage_col} column found, using default 'embryo_model'")
    
    # Handle orig.ident column with validation
    if orig_ident_col not in meta.columns:
        # If no orig.ident column, use file basename
        meta[orig_ident_col] = file_basename
        log_message(f"  No {orig_ident_col} column found, using file basename: {file_basename}")
    else:
        # Check if orig.ident values are valid (contain English letters)
        # Convert to string first to handle categorical data
        orig_ident_str = meta[orig_ident_col].astype(str)
        valid_orig_idents = orig_ident_str.apply(is_valid_identifier)
        
        # Convert to numpy array to avoid categorical issues with logical operations
        valid_orig_idents_array = valid_orig_idents.values
        invalid_count = (~valid_orig_idents_array).sum()
        
        if invalid_count > 0:
            log_message(f"  Found {invalid_count} invalid {orig_ident_col} values, replacing with file basename")
            # Use boolean indexing with the array
            meta.loc[~valid_orig_idents_array, orig_ident_col] = file_basename
        
        # If all values are invalid, replace all with file basename
        if invalid_count == len(meta):
            log_message(f"  All {orig_ident_col} values are invalid, using file basename: {file_basename}")
            meta[orig_ident_col] = file_basename
    
    log_message(f"  Final orig.ident values: {meta[orig_ident_col].unique()}")
    
    # Calculate frequency and percentage (equivalent to R code)
    freq_summary = (meta.groupby([stage_col, orig_ident_col, cell_type_col])
                   .size()
                   .reset_index(name='frequency'))
    
    # Calculate total per stage and orig.ident
    totals = (meta.groupby([stage_col, orig_ident_col])
             .size()
             .reset_index(name='total'))
    
    # Merge and calculate percentages
    freq_summary = freq_summary.merge(totals, on=[stage_col, orig_ident_col])
    freq_summary['percentage'] = (freq_summary['frequency'] / freq_summary['total']) * 100
    
    # Add dataset identifier
    freq_summary['dataset'] = dataset_name
    freq_summary['file_basename'] = file_basename
    
    # Sort as in R code
    freq_summary = freq_summary.sort_values([stage_col, orig_ident_col, cell_type_col])
    
    log_message(f"  Found {len(freq_summary)} cell type-stage-batch combinations")
    log_message(f"  Cell types: {sorted(freq_summary[cell_type_col].unique())}")
    log_message(f"  Stage-orig.ident combinations: {sorted(freq_summary[stage_col].astype(str) + '_' + freq_summary[orig_ident_col].astype(str))}")
    
    return freq_summary

def load_reference_composition(reference_file, cell_type_col='reanno', 
                              stage_col='stage', orig_ident_col='orig.ident'):
    """
    Load and calculate cell type composition for reference dataset.
    
    Parameters:
    -----------
    reference_file : str
        Path to reference h5ad file
    cell_type_col : str
        Column name for cell type annotations
    stage_col : str
        Column name for stage information
    orig_ident_col : str
        Column name for original identity/batch information
    
    Returns:
    --------
    pd.DataFrame: DataFrame with reference cell type composition percentages
    """
    
    log_message(f"Loading reference data from {reference_file}")
    
    try:
        reference_adata = sc.read_h5ad(reference_file)
        log_message(f"Loaded reference with {reference_adata.n_obs} cells and {reference_adata.n_vars} genes")
        
        # Get file basename for reference
        ref_basename = Path(reference_file).stem
        
        # Calculate composition for reference (don't filter for consistency)
        ref_composition = calculate_cell_type_composition(
            reference_adata, 
            'reference', 
            ref_basename,
            cell_type_col=cell_type_col,
            stage_col=stage_col,
            orig_ident_col=orig_ident_col,
            use_consistent_cells=False  # Don't filter reference
        )
        
        return ref_composition
        
    except Exception as e:
        log_message(f"Error loading reference data: {str(e)}")
        return pd.DataFrame()

def process_embryo_model_datasets(data_dir, output_dir, cell_type_col='reanno_pred', 
                                 use_consistent_cells=True):
    """
    Process all embryo model datasets and calculate compositions.
    
    Parameters:
    -----------
    data_dir : str
        Directory containing embryo model h5ad files
    output_dir : str
        Output directory for results
    cell_type_col : str
        Column name for cell type predictions in embryo models
    use_consistent_cells : bool
        Whether to use only consistent cells
    
    Returns:
    --------
    list: List of composition DataFrames for all datasets
    """
    
    data_dir = Path(data_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Find all h5ad files
    h5ad_files = list(data_dir.glob("*.h5ad"))
    log_message(f"Found {len(h5ad_files)} h5ad files in {data_dir}")
    
    all_compositions = []
    
    for file_path in h5ad_files:
        dataset_name = file_path.stem
        file_basename = file_path.stem  # Use stem as basename
        log_message(f"\nProcessing {dataset_name}")
        
        try:
            # Load dataset
            adata = sc.read_h5ad(file_path)
            
            # Calculate composition
            composition = calculate_cell_type_composition(
                adata, 
                dataset_name,
                file_basename,
                cell_type_col=cell_type_col,
                stage_col='stage',  # May not exist, will be handled
                orig_ident_col='orig.ident',  # May not exist, will be handled
                use_consistent_cells=use_consistent_cells
            )
            
            if not composition.empty:
                all_compositions.append(composition)
                
                # Save individual dataset composition
                csv_path = output_dir / f"{dataset_name}_composition.csv"
                composition.to_csv(csv_path, index=False)
                log_message(f"  Saved composition to {csv_path}")
            else:
                log_message(f"  No composition data generated for {dataset_name}")
                
        except Exception as e:
            log_message(f"Error processing {dataset_name}: {str(e)}")
            import traceback
            log_message(traceback.format_exc())
            continue
    
    return all_compositions

def standardize_cell_types(all_compositions, reference_composition):
    """
    Standardize cell type columns across all datasets to ensure compatibility.
    
    Parameters:
    -----------
    all_compositions : list
        List of composition DataFrames from embryo models
    reference_composition : pd.DataFrame
        Reference composition DataFrame
    
    Returns:
    --------
    tuple: (standardized_all_compositions, standardized_reference_composition)
    """
    
    log_message("Standardizing cell type columns across datasets")
    
    # Collect all cell type columns
    all_cell_type_cols = set()
    
    # From reference
    if not reference_composition.empty:
        ref_cell_type_cols = [col for col in reference_composition.columns if 'reanno' in col and col not in ['stage', 'orig.ident', 'dataset', 'file_basename']]
        all_cell_type_cols.update(ref_cell_type_cols)
        log_message(f"Reference cell type columns: {ref_cell_type_cols}")
    
    # From embryo models
    for i, comp in enumerate(all_compositions):
        model_cell_type_cols = [col for col in comp.columns if 'reanno' in col and col not in ['stage', 'orig.ident', 'dataset', 'file_basename']]
        all_cell_type_cols.update(model_cell_type_cols)
        log_message(f"Model {i} cell type columns: {model_cell_type_cols}")
    
    log_message(f"All cell type columns found: {all_cell_type_cols}")
    
    # Use a standard column name for cell types
    standard_cell_type_col = 'cell_type'
    
    # Standardize reference
    if not reference_composition.empty:
        ref_standardized = reference_composition.copy()
        # Find the cell type column in reference
        ref_cell_type_col = [col for col in reference_composition.columns if 'reanno' in col and col not in ['stage', 'orig.ident', 'dataset', 'file_basename']][0]
        ref_standardized = ref_standardized.rename(columns={ref_cell_type_col: standard_cell_type_col})
        log_message(f"Standardized reference: renamed {ref_cell_type_col} to {standard_cell_type_col}")
    else:
        ref_standardized = reference_composition
    
    # Standardize embryo models
    models_standardized = []
    for i, comp in enumerate(all_compositions):
        comp_standardized = comp.copy()
        # Find the cell type column in this composition
        model_cell_type_cols = [col for col in comp.columns if 'reanno' in col and col not in ['stage', 'orig.ident', 'dataset', 'file_basename']]
        if model_cell_type_cols:
            model_cell_type_col = model_cell_type_cols[0]
            comp_standardized = comp_standardized.rename(columns={model_cell_type_col: standard_cell_type_col})
            log_message(f"Standardized model {i}: renamed {model_cell_type_col} to {standard_cell_type_col}")
        models_standardized.append(comp_standardized)
    
    return models_standardized, ref_standardized

def create_wide_format_and_correlation(all_compositions, reference_composition, output_dir):
    """
    Create wide format data and calculate correlation matrix between dataset-stage combinations.
    
    Parameters:
    -----------
    all_compositions : list
        List of composition DataFrames
    reference_composition : pd.DataFrame
        Reference composition data
    output_dir : Path
        Output directory
    
    Returns:
    --------
    tuple: (wide_df, correlation_matrix)
    """
    
    log_message("Creating wide format data and calculating correlations between dataset-stage combinations")
    
    # Standardize cell type columns
    all_compositions_std, reference_composition_std = standardize_cell_types(all_compositions, reference_composition)
    
    # Combine all compositions
    all_data_list = []
    
    if not reference_composition_std.empty:
        all_data_list.append(reference_composition_std)
        log_message(f"Added reference data with {len(reference_composition_std)} rows")
    
    for i, comp in enumerate(all_compositions_std):
        if not comp.empty:
            all_data_list.append(comp)
            log_message(f"Added model {i} data with {len(comp)} rows")
    
    if not all_data_list:
        log_message("Error: No composition data available")
        return pd.DataFrame(), pd.DataFrame()
    
    combined_data = pd.concat(all_data_list, ignore_index=True)
    log_message(f"Combined data shape: {combined_data.shape}")
    log_message(f"Combined data columns: {list(combined_data.columns)}")
    
    # Create stage_orig.ident identifier (this will be our correlation matrix row/column names)
    combined_data['stage_orig_ident'] = (combined_data['stage'].astype(str) + '_' + 
                                        combined_data['orig.ident'].astype(str))
    
    log_message(f"Unique stage_orig_ident combinations: {sorted(combined_data['stage_orig_ident'].unique())}")
    
    # Use the standardized cell type column
    cell_type_col = 'cell_type'
    
    if cell_type_col not in combined_data.columns:
        log_message("Error: No standardized cell type column found")
        return pd.DataFrame(), pd.DataFrame()
    
    log_message(f"Using cell type column: {cell_type_col}")
    log_message(f"Unique cell types: {sorted(combined_data[cell_type_col].unique())}")
    
    # Reshape to wide format: rows = stage_orig_ident combinations, columns = cell types
    # This is the CORRECT orientation for calculating correlations between dataset-stage combinations
    wide_df = combined_data.pivot_table(
        index='stage_orig_ident',  # Each row is a dataset-stage combination
        columns=cell_type_col,     # Each column is a cell type
        values='percentage',
        fill_value=0
    )
    
    log_message(f"Wide format shape: {wide_df.shape}")
    log_message(f"Dataset-stage combinations: {len(wide_df.index)}")
    log_message(f"Cell types: {len(wide_df.columns)}")
    log_message(f"Dataset-stage combinations: {list(wide_df.index)}")
    
    # Save wide format data
    wide_csv_path = output_dir / "cell_type_composition_wide.csv"
    wide_df.to_csv(wide_csv_path)
    log_message(f"Saved wide format data to {wide_csv_path}")
    
    # Calculate correlation matrix between dataset-stage combinations
    # Each row in wide_df is a cell type composition vector for a dataset-stage combination
    # We want to correlate these vectors (i.e., correlate rows)
    correlation_matrix = wide_df.T.corr()  # Transpose and correlate to get correlations between rows (dataset-stage combinations)
    
    log_message(f"Correlation matrix shape: {correlation_matrix.shape}")
    log_message(f"Correlation matrix index: {list(correlation_matrix.index)}")
    
    # Save correlation matrix
    corr_csv_path = output_dir / "composition_correlation_matrix.csv"
    correlation_matrix.to_csv(corr_csv_path)
    log_message(f"Saved correlation matrix to {corr_csv_path}")
    
    return wide_df, correlation_matrix

def create_correlation_heatmap(correlation_matrix, output_dir, figsize=(15, 12)):
    """
    Create correlation heatmap with hierarchical clustering.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Correlation matrix
    output_dir : Path
        Output directory
    figsize : tuple
        Figure size
    """
    
    log_message("Creating correlation heatmap")
    
    if correlation_matrix.empty:
        log_message("Error: Empty correlation matrix")
        return
    
    # Create figure with clustering
    plt.figure(figsize=figsize)
    
    # Create clustermap (equivalent to ComplexHeatmap with clustering)
    g = sns.clustermap(
        correlation_matrix,
        cmap='RdBu_r',  # Red-Blue colormap (reversed)
        center=0.5,     # Center colormap at 0.5
        vmin=0,         # Minimum value
        vmax=1,         # Maximum value
        annot=False,    # Don't annotate with values (too crowded)
        fmt='.2f',      # Format for annotations
        cbar_kws={'label': 'Correlation'},
        dendrogram_ratio=0.1,
        figsize=figsize,
        linewidths=0.5
    )
    
    # Rotate labels for better visibility
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)
    
    # Add title
    g.fig.suptitle('Cell Type Composition Correlation Matrix\n(Dataset-Stage Combinations)', 
                   fontsize=14, y=0.98)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save heatmap
    heatmap_path = output_dir / "composition_correlation_heatmap.pdf"
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    log_message(f"Saved clustered heatmap to {heatmap_path}")
    
    # Also save as PNG
    heatmap_png_path = output_dir / "composition_correlation_heatmap.png"
    plt.savefig(heatmap_png_path, dpi=300, bbox_inches='tight')
    log_message(f"Saved clustered heatmap to {heatmap_png_path}")
    
    plt.close()
    
    # Create a simpler heatmap without clustering for easier interpretation
    plt.figure(figsize=figsize)
    
    # Separate reference and embryo model combinations
    ref_indices = [idx for idx in correlation_matrix.index if any(stage in idx for stage in ['E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E14', 'CS'])]
    model_indices = [idx for idx in correlation_matrix.index if idx not in ref_indices]
    
    log_message(f"Reference indices: {ref_indices}")
    log_message(f"Model indices: {model_indices}")
    
    # Reorder: reference first, then models
    ordered_indices = sorted(ref_indices) + sorted(model_indices)
    ordered_matrix = correlation_matrix.loc[ordered_indices, ordered_indices]
    
    # Create heatmap
    sns.heatmap(
        ordered_matrix,
        cmap='RdBu_r',
        center=0.5,
        vmin=0,
        vmax=1,
        annot=False,  # Too many values to annotate clearly
        fmt='.2f',
        cbar_kws={'label': 'Correlation'},
        square=True,
        linewidths=0.5
    )
    
    plt.title('Cell Type Composition Correlation Matrix\n(Dataset-Stage Combinations - Ordered)', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(rotation=0, fontsize=10)
    plt.tight_layout()
    
    # Save ordered heatmap
    ordered_heatmap_path = output_dir / "composition_correlation_heatmap_ordered.pdf"
    plt.savefig(ordered_heatmap_path, dpi=300, bbox_inches='tight')
    log_message(f"Saved ordered heatmap to {ordered_heatmap_path}")
    
    ordered_heatmap_png_path = output_dir / "composition_correlation_heatmap_ordered.png"
    plt.savefig(ordered_heatmap_png_path, dpi=300, bbox_inches='tight')
    log_message(f"Saved ordered heatmap to {ordered_heatmap_png_path}")
    
    plt.close()

def find_best_matches(correlation_matrix, output_dir):
    """
    Find best matching reference stages for each embryo model dataset-stage combination.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Correlation matrix
    output_dir : Path
        Output directory
    """
    
    log_message("Finding best matches between embryo models and reference dataset-stage combinations")
    
    # Separate reference and model combinations
    ref_indices = [idx for idx in correlation_matrix.index if any(stage in idx for stage in ['E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E14', 'CS'])]
    model_indices = [idx for idx in correlation_matrix.index if idx not in ref_indices]
    
    if not ref_indices or not model_indices:
        log_message("Warning: No reference or model combinations found for comparison")
        log_message(f"All indices: {list(correlation_matrix.index)}")
        return
    
    log_message(f"Reference combinations: {ref_indices}")
    log_message(f"Model combinations: {model_indices}")
    
    # Extract correlations between models and reference
    model_ref_corr = correlation_matrix.loc[model_indices, ref_indices]
    
    # Find best matches
    best_matches = []
    
    for model in model_indices:
        best_ref = model_ref_corr.loc[model].idxmax()
        best_corr = model_ref_corr.loc[model, best_ref]
        
        best_matches.append({
            'embryo_model_combination': model,
            'best_matching_reference': best_ref,
            'correlation': best_corr
        })
        
        log_message(f"  {model} -> {best_ref} (r = {best_corr:.3f})")
    
    # Create DataFrame and save
    best_matches_df = pd.DataFrame(best_matches)
    best_matches_df = best_matches_df.sort_values('correlation', ascending=False)
    
    matches_path = output_dir / "best_reference_matches.csv"
    best_matches_df.to_csv(matches_path, index=False)
    log_message(f"Saved best matches to {matches_path}")
    
    return best_matches_df

def main():
    """Main function to run the cell type composition similarity analysis."""
    
    # Configuration
    data_dir = "/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output"  # Change this to your embryo model data directory
    reference_file = "/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/data/human_clustering_20250724_v3.h5ad"  # Change this to your reference file
    output_dir = "./cell_type_composition_analysis"
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    log_message("="*60)
    log_message("CELL TYPE COMPOSITION SIMILARITY ANALYSIS (FIXED)")
    log_message("Dataset-Stage Combination Correlations")
    log_message("="*60)
    log_message(f"Data directory: {data_dir}")
    log_message(f"Reference file: {reference_file}")
    log_message(f"Output directory: {output_dir}")
    
    # Step 1: Load reference composition
    log_message("\nStep 1: Loading reference composition")
    reference_composition = load_reference_composition(
        reference_file,
        cell_type_col='reanno',  # Use 'reanno' for reference
        stage_col='stage',
        orig_ident_col='orig.ident'
    )
    
    # Step 2: Process embryo model datasets
    log_message("\nStep 2: Processing embryo model datasets")
    all_compositions = process_embryo_model_datasets(
        data_dir,
        output_dir,
        cell_type_col='reanno_pred',  # Use 'reanno_pred' for embryo models
        use_consistent_cells=True
    )
    
    if not all_compositions:
        log_message("Error: No embryo model compositions generated")
        return
    
    # Step 3: Create wide format and correlation matrix
    log_message("\nStep 3: Creating correlation matrix between dataset-stage combinations")
    wide_df, correlation_matrix = create_wide_format_and_correlation(
        all_compositions,
        reference_composition,
        output_dir
    )
    
    if correlation_matrix.empty:
        log_message("Error: Could not create correlation matrix")
        return
    
    # Step 4: Create heatmap
    log_message("\nStep 4: Creating correlation heatmap")
    create_correlation_heatmap(correlation_matrix, output_dir)
    
    # Step 5: Find best matches
    log_message("\nStep 5: Finding best reference matches")
    best_matches = find_best_matches(correlation_matrix, output_dir)
    
    # Summary
    log_message("\n" + "="*60)
    log_message("ANALYSIS COMPLETE")
    log_message("="*60)
    log_message(f"Processed {len(all_compositions)} embryo model datasets")
    log_message(f"Correlation matrix dimensions: {correlation_matrix.shape}")
    log_message(f"Generated files in: {output_dir}")
    log_message("Output files:")
    log_message("  - Individual dataset compositions: *_composition.csv")
    log_message("  - Wide format data: cell_type_composition_wide.csv")
    log_message("  - Correlation matrix: composition_correlation_matrix.csv")
    log_message("  - Heatmaps: composition_correlation_heatmap.pdf/png")
    log_message("  - Best matches: best_reference_matches.csv")
    log_message("\nCorrelation matrix rows/columns represent dataset-stage combinations (e.g., 'E14_IVC_Zhou')")

if __name__ == "__main__":
    main()

