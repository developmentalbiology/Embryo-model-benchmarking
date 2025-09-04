#!/usr/bin/env python
# coding: utf-8

"""
Independent script to calculate prediction consistency statistics across datasets.
Generates pred_consistency_comparison.csv for benchmarking purposes.

This script calculates the percentage of consistent cells per dataset, where
consistent cells are defined as cells where reanno_pred_lineage == lineage_pred.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
import time

def log_message(message, log_file='consistency_calculation.log'):
    """Log message to console and file"""
    print(message, flush=True)
    with open(log_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

def calculate_consistency_metrics(adata, dataset_name):
    """
    Calculate consistency metrics for a single dataset.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object
    dataset_name : str
        Name of the dataset for logging
    
    Returns:
    --------
    dict: Dictionary containing consistency metrics
    """
    
    log_message(f"Calculating consistency metrics for {dataset_name}")
    
    # Check if required columns exist
    required_cols = ['reanno_pred_lineage', 'lineage_pred']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        error_msg = f"Missing columns: {missing_cols}"
        log_message(f"Warning: {error_msg} in {dataset_name}")
        return {
            'dataset': dataset_name,
            'total_cells': len(adata),
            'consistent_cells': 0,
            'inconsistent_cells': 0,
            'consistency_percentage': 0.0,
            'inconsistency_percentage': 0.0,
            'error': error_msg
        }
    
    # Convert categorical columns to string to avoid category comparison issues
    try:
        reanno_pred_lineage_str = adata.obs['reanno_pred_lineage'].astype(str)
        lineage_pred_str = adata.obs['lineage_pred'].astype(str)
        
        # Calculate consistency
        consistent_cells = (reanno_pred_lineage_str == lineage_pred_str)
        consistent_count = consistent_cells.sum()
        inconsistent_count = (~consistent_cells).sum()
        total_count = len(adata)
        
        consistency_percentage = (consistent_count / total_count * 100) if total_count > 0 else 0.0
        inconsistency_percentage = (inconsistent_count / total_count * 100) if total_count > 0 else 0.0
        
        log_message(f"  {dataset_name}: {consistent_count}/{total_count} consistent cells ({consistency_percentage:.1f}%)")
        
        return {
            'dataset': dataset_name,
            'total_cells': total_count,
            'consistent_cells': consistent_count,
            'inconsistent_cells': inconsistent_count,
            'consistency_percentage': consistency_percentage,
            'inconsistency_percentage': inconsistency_percentage,
            'error': None
        }
        
    except Exception as e:
        error_msg = f"Error calculating consistency: {str(e)}"
        log_message(f"Error in {dataset_name}: {error_msg}")
        return {
            'dataset': dataset_name,
            'total_cells': len(adata),
            'consistent_cells': 0,
            'inconsistent_cells': 0,
            'consistency_percentage': 0.0,
            'inconsistency_percentage': 0.0,
            'error': error_msg
        }

def get_lineage_consistency_breakdown(adata, dataset_name):
    """
    Get detailed breakdown of consistency by lineage.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object
    dataset_name : str
        Name of the dataset
    
    Returns:
    --------
    pd.DataFrame: DataFrame with consistency breakdown by lineage
    """
    
    # Check if required columns exist
    required_cols = ['reanno_pred_lineage', 'lineage_pred']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        log_message(f"Cannot create lineage breakdown for {dataset_name}: missing columns {missing_cols}")
        return pd.DataFrame()
    
    try:
        # Convert to string to avoid categorical comparison issues
        reanno_pred_lineage_str = adata.obs['reanno_pred_lineage'].astype(str)
        lineage_pred_str = adata.obs['lineage_pred'].astype(str)
        
        # Calculate consistency
        consistent_cells = (reanno_pred_lineage_str == lineage_pred_str)
        
        # Create breakdown by lineage_pred
        breakdown_data = []
        
        for lineage in adata.obs['lineage_pred'].unique():
            if pd.isna(lineage):
                continue
                
            lineage_mask = (adata.obs['lineage_pred'] == lineage)
            lineage_total = lineage_mask.sum()
            lineage_consistent = (consistent_cells & lineage_mask).sum()
            lineage_consistency_pct = (lineage_consistent / lineage_total * 100) if lineage_total > 0 else 0.0
            
            breakdown_data.append({
                'dataset': dataset_name,
                'lineage': str(lineage),
                'total_cells': lineage_total,
                'consistent_cells': lineage_consistent,
                'consistency_percentage': lineage_consistency_pct
            })
        
        return pd.DataFrame(breakdown_data)
        
    except Exception as e:
        log_message(f"Error creating lineage breakdown for {dataset_name}: {str(e)}")
        return pd.DataFrame()

def process_h5ad_files_for_consistency(folder_path, output_folder=None):
    """
    Process all h5ad files in a folder and calculate consistency metrics.
    
    Parameters:
    -----------
    folder_path : str
        Path to folder containing h5ad files
    output_folder : str, optional
        Path to output folder (if None, uses the same folder as input)
    
    Returns:
    --------
    tuple: (consistency_df, lineage_breakdown_df)
    """
    
    folder_path = Path(folder_path)
    
    # Set output path
    if output_folder:
        output_path = Path(output_folder)
        output_path.mkdir(exist_ok=True)
    else:
        output_path = folder_path
    
    # Get all h5ad files in the folder
    h5ad_files = list(folder_path.glob("*.h5ad"))
    
    if not h5ad_files:
        log_message(f"No h5ad files found in {folder_path}")
        return pd.DataFrame(), pd.DataFrame()
    
    log_message(f"Found {len(h5ad_files)} h5ad files in {folder_path}")
    
    # List to store consistency metrics for all datasets
    consistency_metrics_list = []
    lineage_breakdown_list = []
    
    for file_path in h5ad_files:
        dataset_name = file_path.stem
        log_message(f"\nProcessing: {file_path.name}")
        
        try:
            # Load the h5ad file
            adata = sc.read_h5ad(file_path)
            log_message(f"Loaded {dataset_name} with {adata.n_obs} cells and {adata.n_vars} genes")
            
            # Calculate consistency metrics
            metrics = calculate_consistency_metrics(adata, dataset_name)
            metrics['filename'] = file_path.name
            consistency_metrics_list.append(metrics)
            
            # Get lineage breakdown
            lineage_breakdown = get_lineage_consistency_breakdown(adata, dataset_name)
            if not lineage_breakdown.empty:
                lineage_breakdown_list.append(lineage_breakdown)
            
        except Exception as e:
            log_message(f"Error processing {file_path.name}: {str(e)}")
            # Still add an entry for failed files
            consistency_metrics_list.append({
                'dataset': dataset_name,
                'filename': file_path.name,
                'total_cells': 0,
                'consistent_cells': 0,
                'inconsistent_cells': 0,
                'consistency_percentage': 0.0,
                'inconsistency_percentage': 0.0,
                'error': str(e)
            })
            continue
    
    # Create DataFrame with consistency metrics
    consistency_df = pd.DataFrame(consistency_metrics_list)
    
    # Sort by consistency percentage (descending)
    consistency_df = consistency_df.sort_values('consistency_percentage', ascending=False)
    
    # Combine lineage breakdowns
    if lineage_breakdown_list:
        lineage_breakdown_df = pd.concat(lineage_breakdown_list, ignore_index=True)
    else:
        lineage_breakdown_df = pd.DataFrame()
    
    # Save consistency comparison CSV
    csv_output_path = output_path / "pred_consistency_comparison.csv"
    consistency_df.to_csv(csv_output_path, index=False)
    log_message(f"Saved consistency comparison to: {csv_output_path}")
    
    # Save lineage breakdown CSV if available
    if not lineage_breakdown_df.empty:
        lineage_csv_path = output_path / "pred_consistency_by_lineage.csv"
        lineage_breakdown_df.to_csv(lineage_csv_path, index=False)
        log_message(f"Saved lineage breakdown to: {lineage_csv_path}")
    
    return consistency_df, lineage_breakdown_df

def print_consistency_summary(consistency_df):
    """
    Print a detailed summary of consistency statistics.
    
    Parameters:
    -----------
    consistency_df : pd.DataFrame
        DataFrame with consistency metrics
    """
    
    log_message("\n" + "="*60)
    log_message("PREDICTION CONSISTENCY SUMMARY")
    log_message("="*60)
    
    log_message(f"Total datasets processed: {len(consistency_df)}")
    
    # Filter out error cases for statistics
    valid_df = consistency_df[consistency_df['error'].isna()]
    error_df = consistency_df[consistency_df['error'].notna()]
    
    if len(valid_df) > 0:
        log_message(f"Valid datasets: {len(valid_df)}")
        log_message(f"Datasets with errors: {len(error_df)}")
        
        # Overall statistics
        log_message(f"\nOVERALL STATISTICS:")
        log_message(f"Mean consistency: {valid_df['consistency_percentage'].mean():.1f}%")
        log_message(f"Median consistency: {valid_df['consistency_percentage'].median():.1f}%")
        log_message(f"Standard deviation: {valid_df['consistency_percentage'].std():.1f}%")
        log_message(f"Min consistency: {valid_df['consistency_percentage'].min():.1f}%")
        log_message(f"Max consistency: {valid_df['consistency_percentage'].max():.1f}%")
        
        # Total cells statistics
        total_cells = valid_df['total_cells'].sum()
        total_consistent = valid_df['consistent_cells'].sum()
        overall_consistency = (total_consistent / total_cells * 100) if total_cells > 0 else 0.0
        log_message(f"\nAGGREGATE STATISTICS:")
        log_message(f"Total cells across all datasets: {total_cells:,}")
        log_message(f"Total consistent cells: {total_consistent:,}")
        log_message(f"Overall consistency rate: {overall_consistency:.1f}%")
        
        # Top performers
        log_message(f"\nTOP 5 MOST CONSISTENT DATASETS:")
        for i, (_, row) in enumerate(valid_df.head(5).iterrows(), 1):
            log_message(f"  {i}. {row['dataset']}: {row['consistency_percentage']:.1f}% ({row['consistent_cells']:,}/{row['total_cells']:,})")
        
        # Bottom performers
        if len(valid_df) > 5:
            log_message(f"\nBOTTOM 5 LEAST CONSISTENT DATASETS:")
            for i, (_, row) in enumerate(valid_df.tail(5).iterrows(), 1):
                log_message(f"  {i}. {row['dataset']}: {row['consistency_percentage']:.1f}% ({row['consistent_cells']:,}/{row['total_cells']:,})")
        
        # Consistency distribution
        log_message(f"\nCONSISTENCY DISTRIBUTION:")
        bins = [0, 20, 40, 60, 80, 100]
        for i in range(len(bins)-1):
            count = ((valid_df['consistency_percentage'] >= bins[i]) & 
                    (valid_df['consistency_percentage'] < bins[i+1])).sum()
            if i == len(bins)-2:  # Last bin includes 100%
                count = ((valid_df['consistency_percentage'] >= bins[i]) & 
                        (valid_df['consistency_percentage'] <= bins[i+1])).sum()
            log_message(f"  {bins[i]}-{bins[i+1]}%: {count} datasets")
    
    # Report any errors
    if len(error_df) > 0:
        log_message(f"\nDATASETS WITH ERRORS:")
        for _, row in error_df.iterrows():
            log_message(f"  {row['dataset']}: {row['error']}")
    
    log_message("="*60)

def main():
    """Main function to run the consistency analysis."""
    
    # Configuration
    folder_path = "/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output"  # Change this to your data directory
    output_folder = None  # Set to a specific path if you want output in a different location
    
    log_message("Starting prediction consistency analysis...")
    log_message(f"Input folder: {folder_path}")
    
    # Process all h5ad files and calculate consistency metrics
    consistency_df, lineage_breakdown_df = process_h5ad_files_for_consistency(
        folder_path=folder_path,
        output_folder=output_folder
    )
    
    if consistency_df.empty:
        log_message("No data processed. Exiting.")
        return
    
    # Print detailed summary
    print_consistency_summary(consistency_df)
    
    # Display the final table
    log_message("\nFINAL CONSISTENCY COMPARISON TABLE:")
    log_message("-" * 100)
    
    # Format the display
    display_df = consistency_df.copy()
    
    # Format columns for better display
    for col in ['consistency_percentage', 'inconsistency_percentage']:
        if col in display_df.columns:
            display_df[col] = display_df[col].round(1)
    
    # Select key columns for display
    display_cols = ['dataset', 'total_cells', 'consistent_cells', 'consistency_percentage']
    if 'error' in display_df.columns:
        # Only show error column if there are errors
        if display_df['error'].notna().any():
            display_cols.append('error')
    
    log_message(display_df[display_cols].to_string(index=False))
    
    log_message(f"\nAnalysis complete! Results saved to:")
    output_path = Path(output_folder) if output_folder else Path(folder_path)
    log_message(f"  - {output_path / 'pred_consistency_comparison.csv'}")
    if not lineage_breakdown_df.empty:
        log_message(f"  - {output_path / 'pred_consistency_by_lineage.csv'}")

if __name__ == "__main__":
    main()

