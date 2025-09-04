#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import scib

import warnings
warnings.filterwarnings('ignore')

def log_message(message):
    """Log message with timestamp"""
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def calculate_bio_conservation(adata, umap_config, batch_key="orig.ident", 
                              subsample=None):
    """
    Calculate biological conservation metrics for different UMAP embeddings
    
    Parameters:
    -----------
    adata: AnnData object
    umap_config: Dictionary mapping UMAP keys to list of (lineage_key, cluster_key) tuples
    batch_key: Key for batch information in adata.obs (default: "orig.ident")
    subsample: Number of cells to subsample for faster computation (None = use all)
    
    Returns:
    --------
    DataFrame with metrics for each UMAP version
    """
    
    # Initialize results dictionary - only keeping Silhouette, NMI, ARI
    results = {
        'UMAP': [],
        'Lineage_key': [],
        'Cluster_key': [],
        'Silhouette': [],
        'NMI': [],
        'ARI': []
    }
    
    # Check if batch_key exists in adata.obs
    if batch_key not in adata.obs:
        print(f"Warning: Batch key '{batch_key}' not found in adata.obs. Some metrics may fail.")
    
    # Process each UMAP configuration
    for umap_key in umap_config:
        for (label_key, cluster_key) in umap_config[umap_key]:
            start_time = time.time()
            
            # Check if UMAP key exists in adata.obsm
            if umap_key not in adata.obsm:
                print(f"Warning: UMAP key '{umap_key}' not found in adata.obsm. Skipping.")
                continue

            # Check if label_key exists in adata.obs
            if label_key not in adata.obs:
                print(f"Warning: Label key '{label_key}' not found in adata.obs. Skipping.")
                continue

            # Check if cluster_key exists in adata.obs
            if cluster_key not in adata.obs:
                print(f"Warning: Cluster key '{cluster_key}' not found in adata.obs. Skipping.")
                continue

            print(f"Processing {umap_key} with lineage: {label_key}, clusters: {cluster_key}")
            results['UMAP'].append(umap_key)
            results['Lineage_key'].append(label_key)
            results['Cluster_key'].append(cluster_key)

            # Filter out cells with missing or NA labels or clusters
            mask = ~(adata.obs[label_key].isna() | adata.obs[cluster_key].isna())
            if batch_key in adata.obs:
                mask = mask & ~(adata.obs[batch_key].isna())
            adata_filtered = adata[mask].copy()

            # Subsample if needed (for speed)
            if subsample is not None and subsample < adata_filtered.n_obs:
                print(f"  Subsampling from {adata_filtered.n_obs} to {subsample} cells")
                sc.pp.subsample(adata_filtered, n_obs=subsample, random_state=42)

            print(f"  Working with {adata_filtered.n_obs} cells")

            # Calculate Silhouette score (using existing UMAP coordinates)
            try:
                # Create a temporary AnnData with the correct UMAP key for scib
                temp_adata = adata_filtered.copy()
                temp_adata.obsm['X_emb'] = temp_adata.obsm[umap_key].copy()
                
                # Use scib silhouette method
                sil_score = scib.metrics.silhouette(temp_adata, label_key=label_key, embed='X_emb')
                
                results['Silhouette'].append(sil_score)
                print(f"  Silhouette score: {sil_score}")
            except Exception as e:
                print(f"  Error calculating silhouette for {umap_key}: {e}")
                results['Silhouette'].append(np.nan)

            # Calculate NMI and ARI
            try:
                nmi_score = scib.metrics.nmi(adata_filtered, cluster_key=cluster_key, label_key=label_key)
                results['NMI'].append(nmi_score)
                print(f"  NMI score: {nmi_score}")
            except Exception as e:
                print(f"  Error calculating NMI for {umap_key}: {e}")
                results['NMI'].append(np.nan)

            try:
                ari_score = scib.metrics.ari(adata_filtered, cluster_key=cluster_key, label_key=label_key)
                results['ARI'].append(ari_score)
                print(f"  ARI score: {ari_score}")
            except Exception as e:
                print(f"  Error calculating ARI for {umap_key}: {e}")
                results['ARI'].append(np.nan)

            elapsed_time = time.time() - start_time
            print(f"  Completed in {elapsed_time:.2f} seconds")
    
    return pd.DataFrame(results)

import gc # Import the garbage collector module

def process_single_file(file_path, batch_key="orig.ident", subsample=None):
    """
    Process a single h5ad file and calculate metrics
    """
    log_message(f"Processing file: {os.path.basename(file_path)}")
    try:
        # Load the data
        adata = sc.read_h5ad(file_path)
        log_message(f"Loaded data with shape: {adata.shape}")
        
        # Print available keys for debugging
        print(f"Available obs keys: {list(adata.obs.columns)}")
        print(f"Available obsm keys: {list(adata.obsm.keys())}")
        
        # Define UMAP configuration based on your requirements
        umap_config = {
            'X_scANVI': [
                ('human_ref_lineage_pred', 'scANVI_res_0.5'),
                ('human_ref_reanno_pred', 'scANVI_res_0.5')
            ]
        }
        
        # Calculate metrics
        metrics_df = calculate_bio_conservation(
            adata, 
            umap_config, 
            batch_key=batch_key,
            subsample=subsample
        )
        
        # Add file information
        metrics_df['File'] = os.path.basename(file_path)
        metrics_df['File_path'] = file_path
        
        # --- MEMORY RELEASE ---
        # 1. Delete the large AnnData object
        del adata
        # 2. Force garbage collection to reclaim the memory
        gc.collect()
        # 3. (Optional) Log the current memory usage to confirm release
        # import psutil
        # process = psutil.Process(os.getpid())
        # log_message(f"Memory usage after processing: {process.memory_info().rss / 1024 / 1024:.2f} MB")
        
        log_message(f"Completed processing {os.path.basename(file_path)}")
        return metrics_df
        
    except Exception as e:
        log_message(f"Error processing {file_path}: {str(e)}")
        # Ensure memory is released even if an error occurs
        if 'adata' in locals():
            del adata
        gc.collect()
        return None

def process_folder(folder_path, output_csv="bio_conservation_metrics.csv", 
                  batch_key="orig.ident", subsample=None, file_pattern="*.h5ad"):
    """
    Process all h5ad files in a folder and combine results
    """
    log_message(f"Starting batch processing in folder: {folder_path}")
    
    # Find all h5ad files
    h5ad_files = glob.glob(os.path.join(folder_path, file_pattern))
    
    if not h5ad_files:
        log_message(f"No h5ad files found in {folder_path}")
        return None
    
    log_message(f"Found {len(h5ad_files)} h5ad files to process")
    
    all_results = []
    
    # Process each file
    for i, file_path in enumerate(h5ad_files, 1):
        log_message(f"Processing file {i}/{len(h5ad_files)}: {os.path.basename(file_path)}")
        
        result_df = process_single_file(
            file_path, 
            batch_key=batch_key,
            subsample=subsample
        )
        
        if result_df is not None and not result_df.empty:
            all_results.append(result_df)
        
        log_message(f"Completed file {i}/{len(h5ad_files)}")
    
    # Combine all results
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Reorder columns for better readability
        column_order = ['File', 'UMAP', 'Lineage_key', 'Cluster_key', 'Silhouette', 'NMI', 'ARI', 'File_path']
        
        combined_df = combined_df[column_order]
        
        # Save to CSV
        output_path = os.path.join(folder_path, output_csv)
        combined_df.to_csv(output_path, index=False)
        log_message(f"Results saved to: {output_path}")
        
        # Print summary
        print("\n" + "="*80)
        print("PROCESSING SUMMARY")
        print("="*80)
        print(f"Total files processed: {len(all_results)}")
        print(f"Total metrics calculated: {len(combined_df)}")
        print(f"Output saved to: {output_path}")
        print("\nMetrics summary:")
        print(combined_df.groupby(['UMAP', 'Lineage_key']).agg({
            'Silhouette': ['mean', 'std'],
            'NMI': ['mean', 'std'], 
            'ARI': ['mean', 'std']
        }).round(4))
        print("="*80)
        
        return combined_df
    else:
        log_message("No valid results to combine")
        return None

def main():
    """
    Main function - modify the folder path here
    """
    # TODO: Modify this path to your folder containing h5ad files
    folder_path = "/path/to/your/folder"  # Change this to your actual folder path
    
    # Configuration
    output_csv = "bio_conservation_metrics.csv"
    batch_key = "orig.ident"
    subsample = None  # Set to a number (e.g., 5000) for faster computation on large datasets
    
    # Process all files
    results = process_folder(
        folder_path=folder_path,
        output_csv=output_csv,
        batch_key=batch_key,
        subsample=subsample
    )
    
    return results

if __name__ == "__main__":
    # Example usage with command line arguments
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate biological conservation metrics for h5ad files')
    parser.add_argument('--folder', '-f', required=True, help='Path to folder containing h5ad files')
    parser.add_argument('--output', '-o', default='bio_conservation_metrics.csv', help='Output CSV filename')
    parser.add_argument('--batch_key', '-b', default='orig.ident', help='Batch key in adata.obs')
    parser.add_argument('--subsample', '-s', type=int, default=None, help='Subsample cells for faster computation')
    
    args = parser.parse_args()
    
    # Process files
    results = process_folder(
        folder_path=args.folder,
        output_csv=args.output,
        batch_key=args.batch_key,
        subsample=args.subsample
    )
    
    if results is not None:
        print(f"\nFinal results shape: {results.shape}")
        print("\nFirst few rows:")
        print(results.head())
    else:
        print("No results generated.")