#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
from scipy.sparse import issparse
import os
import matplotlib.pyplot as plt
import numpy as np
import glob
from pathlib import Path
import pandas as pd
from datetime import datetime

def preprocess_data(file_path: str, resolution: float = 1.0) -> dict:
    """
    Preprocesses single-cell data for downstream analysis.
    
    Parameters:
    -----------
    file_path : str
        Path to the input .h5ad file.
    resolution : float, optional
        Resolution parameter for clustering, controls the number of clusters.
        Default is 1.0.
    
    Returns:
    --------
    dict
        Summary information about the preprocessing results.
    """
    try:
        # Set output directory
        output_dir = os.path.dirname(file_path)
        os.makedirs(output_dir, exist_ok=True)

        # Extract the base filename without extension
        base_filename = os.path.splitext(os.path.basename(file_path))[0]
        output_filename = f"{base_filename}_preprocessed.h5ad"
        
        # Load query data
        try:
            query_adata = sc.read_h5ad(file_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {file_path}")
        
        print(f"Query dataset loaded successfully: {base_filename}")
        
        # Sanitize column names in .obs and .var to avoid reserved names
        if "_index" in query_adata.obs.columns:
            query_adata.obs.rename(columns={"_index": "cell_index"}, inplace=True)
        if "_index" in query_adata.var.columns:
            query_adata.var.rename(columns={"_index": "gene_index"}, inplace=True)
        
        # Remove "empty" genes
        sc.pp.filter_genes(query_adata, min_cells=1)

        # Convert adata.X to a numpy array if it's sparse or has 'todense'/'toarray' methods
        if hasattr(query_adata.X, 'todense'):
            X_data = np.asarray(query_adata.X.todense())  # Convert to dense matrix and then to NumPy array
        elif hasattr(query_adata.X, 'toarray'):
            X_data = np.asarray(query_adata.X.toarray())  # Convert to dense matrix and then to NumPy array
        else:
            X_data = np.asarray(query_adata.X)  # Already dense, convert to NumPy array

        # Check for NaN values
        if np.isnan(X_data).any():
            nan_count = np.isnan(X_data).sum()
            print(f"Found {nan_count} NaN values in data matrix")
            replace_nans = True  # Set to False if you want to skip instead

            if replace_nans:
                print(f"Replacing NaN values with zeros...")
                X_data = np.nan_to_num(X_data, nan=0.0)
                query_adata.X = X_data
            else:
                print(f"NaN values found in data, skipping...")
                return {
                    "status": "error",
                    "message": f"Found {nan_count} NaN values in data matrix",
                    "file_path": file_path
                }

        # Preprocess query dataset
        sc.settings.seed = 42  # Set random seed for reproducibility
        
        # Save raw counts
        query_adata.layers["counts"] = query_adata.X.copy()
        
        # Normalize total counts
        sc.pp.normalize_total(query_adata, target_sum=1e4)
        
        # Log-transform the data
        sc.pp.log1p(query_adata)
        query_adata.layers["logcounts"] = query_adata.X.copy()
        
        # Identify highly variable genes
        sc.pp.highly_variable_genes(query_adata, n_top_genes=2000, flavor="cell_ranger")
        
        # Perform PCA
        sc.tl.pca(query_adata, n_comps=30, use_highly_variable=True)
        
        # Compute neighborhood graph
        sc.pp.neighbors(query_adata)

        # Generate UMAP embeddings
        sc.tl.umap(query_adata)
        
        # Cluster cells using Leiden algorithm
        print(f"Clustering with resolution {resolution}...")
        sc.tl.leiden(query_adata, resolution=resolution, key_added=f"leiden_r{resolution}", random_state=42)
        
        # Export UMAP plot with dataset-specific name
        umap_plot_filename = f"{base_filename}_umap_plot_resolution_{resolution}.png"
        umap_plot_path = os.path.join(output_dir, umap_plot_filename)
        plt.figure(figsize=(6, 4))
        sc.pl.umap(query_adata, color=[f"leiden_r{resolution}"], show=False, frameon=False)  
        plt.savefig(umap_plot_path, bbox_inches='tight')
        plt.close()
        print(f"UMAP plot saved at: {umap_plot_path}")
        
        # Save preprocessed AnnData object with dataset-specific name
        preprocessed_h5ad_path = os.path.join(output_dir, output_filename)
        
        # Fix any potential _index issues before saving
        if query_adata.raw is not None:
            if '_index' in query_adata.raw.var.columns:
                query_adata.raw.var.rename(columns={'_index': 'index'}, inplace=True)
        else:
            if '_index' in query_adata.var.columns:
                query_adata.var.rename(columns={'_index': 'index'}, inplace=True)
        
        query_adata.write(preprocessed_h5ad_path)
        print(f"Preprocessed .h5ad file saved at: {preprocessed_h5ad_path}")
        
        # Remove sparsity
        def remove_sparsity(adata):
            if issparse(adata.X):
                adata.X = adata.X.toarray()
            return adata
        
        query_adata = remove_sparsity(query_adata)
        
        # Return summary information
        clusters = query_adata.obs[f'leiden_r{resolution}'].value_counts().to_dict()
        
        return {
            "Number of observations (cells)": query_adata.n_obs,
            "Number of variables (genes)": query_adata.n_vars,
            "Original shape": query_adata.shape,
            "Processed shape": query_adata.shape,
            "Highly Variable Genes": query_adata.var['highly_variable'].sum(),
            "PCA Components": query_adata.obsm['X_pca'].shape[1],
            "Resolution used": resolution,
            "Leiden Clusters": len(clusters),
            "Cluster sizes": clusters,
            "UMAP Embeddings": query_adata.obsm['X_umap'].shape[1],
            "UMAP Plot Path": umap_plot_path,
            "Preprocessed h5ad Path": preprocessed_h5ad_path,
            "Base Filename": base_filename,
            "Summary": str(query_adata),
            "status": "success",
            "error_message": None
        }
    
    except Exception as e:
        return {
            "status": "error",
            "output": None,
            "plots": [],
            "error_message": str(e),
            "file_path": file_path
        }


def process_folder(folder_path: str, resolution: float = 1.0, output_summary: bool = True) -> dict:
    """
    Process all .h5ad files in a given folder.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing .h5ad files.
    resolution : float, optional
        Resolution parameter for clustering. Default is 1.0.
    output_summary : bool, optional
        Whether to create a summary CSV file. Default is True.
    
    Returns:
    --------
    dict
        Summary of batch processing results.
    """
    folder_path = Path(folder_path)
    
    if not folder_path.exists() or not folder_path.is_dir():
        raise ValueError(f"Invalid folder path: {folder_path}")
    
    # Find all .h5ad files in the folder
    h5ad_files = list(folder_path.glob("*.h5ad"))
    
    if not h5ad_files:
        print(f"No .h5ad files found in {folder_path}")
        return {
            "status": "warning",
            "message": "No .h5ad files found",
            "processed_files": 0,
            "successful_files": 0,
            "failed_files": 0
        }
    
    print(f"Found {len(h5ad_files)} .h5ad files to process...")
    
    results = []
    successful_count = 0
    failed_count = 0
    
    for i, file_path in enumerate(h5ad_files, 1):
        print(f"\n{'='*50}")
        print(f"Processing file {i}/{len(h5ad_files)}: {file_path.name}")
        print(f"{'='*50}")
        
        try:
            result = preprocess_data(str(file_path), resolution)
            
            if result["status"] == "success":
                successful_count += 1
                print(f"✓ Successfully processed: {file_path.name}")
            else:
                failed_count += 1
                print(f"✗ Failed to process: {file_path.name}")
                print(f"  Error: {result.get('error_message', 'Unknown error')}")
            
            # Add file information to result
            result["file_name"] = file_path.name
            result["file_path"] = str(file_path)
            results.append(result)
            
        except Exception as e:
            failed_count += 1
            error_result = {
                "status": "error",
                "file_name": file_path.name,
                "file_path": str(file_path),
                "error_message": str(e)
            }
            results.append(error_result)
            print(f"✗ Failed to process: {file_path.name}")
            print(f"  Error: {str(e)}")
    
    # Create summary report
    if output_summary and results:
        summary_data = []
        for result in results:
            if result["status"] == "success":
                summary_row = {
                    "File Name": result["file_name"],
                    "Status": result["status"],
                    "Cells": result.get("Number of observations (cells)", "N/A"),
                    "Genes": result.get("Number of variables (genes)", "N/A"),
                    "Highly Variable Genes": result.get("Highly Variable Genes", "N/A"),
                    "Clusters": result.get("Leiden Clusters", "N/A"),
                    "Resolution": result.get("Resolution used", "N/A"),
                    "UMAP Plot": result.get("UMAP Plot Path", "N/A"),
                    "Preprocessed File": result.get("Preprocessed h5ad Path", "N/A"),
                    "Error Message": ""
                }
            else:
                summary_row = {
                    "File Name": result["file_name"],
                    "Status": result["status"],
                    "Cells": "N/A",
                    "Genes": "N/A",
                    "Highly Variable Genes": "N/A",
                    "Clusters": "N/A",
                    "Resolution": "N/A",
                    "UMAP Plot": "N/A",
                    "Preprocessed File": "N/A",
                    "Error Message": result.get("error_message", "Unknown error")
                }
            summary_data.append(summary_row)
        
        # Save summary to CSV
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        summary_filename = f"preprocessing_summary_{timestamp}.csv"
        summary_path = folder_path / summary_filename
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(summary_path, index=False)
        print(f"\n📊 Summary report saved to: {summary_path}")
    
    # Print final summary
    print(f"\n{'='*60}")
    print(f"BATCH PROCESSING COMPLETE")
    print(f"{'='*60}")
    print(f"Total files processed: {len(h5ad_files)}")
    print(f"Successful: {successful_count}")
    print(f"Failed: {failed_count}")
    print(f"Success rate: {successful_count/len(h5ad_files)*100:.1f}%")
    
    return {
        "status": "complete",
        "total_files": len(h5ad_files),
        "successful_files": successful_count,
        "failed_files": failed_count,
        "success_rate": successful_count/len(h5ad_files)*100,
        "results": results,
        "summary_file": summary_path if output_summary else None
    }


if __name__ == "__main__":
    import argparse

    # Argument parsing for standalone execution
    parser = argparse.ArgumentParser(description="Preprocess single-cell RNA-seq dataset(s).")
    parser.add_argument("--file_path", type=str, help="Path to a single .h5ad file.")
    parser.add_argument("--folder_path", type=str, help="Path to folder containing multiple .h5ad files.")
    parser.add_argument("--resolution", type=float, default=1.0, help="Clustering resolution parameter.")
    parser.add_argument("--no_summary", action="store_true", help="Skip creating summary CSV file for batch processing.")
    
    args = parser.parse_args()

    # Validate arguments
    if not args.file_path and not args.folder_path:
        parser.error("Either --file_path or --folder_path must be specified.")
    
    if args.file_path and args.folder_path:
        parser.error("Cannot specify both --file_path and --folder_path. Choose one.")

    try:
        if args.file_path:
            # Single file processing
            print(f"Processing single file: {args.file_path}")
            result = preprocess_data(args.file_path, args.resolution)

            # Check if preprocessing was successful
            if result.get("status") == "success":
                # Prepare a concise summary for LLM
                summary = (
                    f"Preprocessing completed successfully for file: {args.file_path}\n"
                    f"- Number of observations (cells): {result['Number of observations (cells)']}\n"
                    f"- Number of variables (genes): {result['Number of variables (genes)']}\n"
                    f"- Highly Variable Genes: {result['Highly Variable Genes']}\n"
                    f"- PCA Components: {result['PCA Components']}\n"
                    f"- Clustering Resolution: {result['Resolution used']}\n"
                    f"- Number of Leiden Clusters: {result['Leiden Clusters']}\n"
                    f"- UMAP Embeddings Dimensions: {result['UMAP Embeddings']}\n"
                    f"- UMAP Plot saved at: {result['UMAP Plot Path']}\n"
                    f"- Preprocessed h5ad file saved at: {result['Preprocessed h5ad Path']}\n"
                )
            else:
                # Handle errors and provide an error message
                summary = f"Error during preprocessing: {result.get('error_message')}"

            # Print the summary as plain text
            print(summary)
            
        elif args.folder_path:
            # Batch processing
            print(f"Processing all .h5ad files in folder: {args.folder_path}")
            batch_result = process_folder(
                args.folder_path, 
                args.resolution, 
                output_summary=not args.no_summary
            )
            
            # Print batch summary
            if batch_result["status"] == "complete":
                print(f"\nBatch processing summary:")
                print(f"- Total files: {batch_result['total_files']}")
                print(f"- Successful: {batch_result['successful_files']}")
                print(f"- Failed: {batch_result['failed_files']}")
                print(f"- Success rate: {batch_result['success_rate']:.1f}%")
                if batch_result.get('summary_file'):
                    print(f"- Summary report: {batch_result['summary_file']}")
            else:
                print(f"Batch processing warning: {batch_result.get('message', 'Unknown issue')}")

    except Exception as e:
        # Handle unexpected errors
        error_message = f"Unexpected error during preprocessing: {str(e)}"
        print(error_message)

# example ussage  
# Single file processing (original functionality):      
#python preprocess_data.py --file_path /path/to/single_file.h5ad --resolution 1.0

# Batch processing all files in a folder:
#python preprocess_data.py --folder_path /path/to/folder/ --resolution 1.0

#run
#python preprocess_data.py --folder_path /storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250801_final_model_test_jingJiang/data

