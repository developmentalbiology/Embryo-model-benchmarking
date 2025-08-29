#!/usr/bin/env python
# coding: utf-8

import os
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
from scarches.dataset.trvae.data_handling import remove_sparsity
import argparse
from pathlib import Path

import warnings
# Suppress warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

def get_h5ad_files(input_path):
    """
    Get list of h5ad files from input path.
    
    Args:
        input_path (str): Path to either a single h5ad file or a directory containing h5ad files
        
    Returns:
        list: List of h5ad file paths
    """
    input_path = Path(input_path)
    
    if input_path.is_file() and input_path.suffix == '.h5ad':
        # Single h5ad file
        return [str(input_path)]
    elif input_path.is_dir():
        # Directory containing h5ad files
        h5ad_files = list(input_path.glob('*.h5ad'))
        if not h5ad_files:
            raise ValueError(f"No .h5ad files found in directory: {input_path}")
        return [str(f) for f in h5ad_files]
    else:
        raise ValueError(f"Input path must be either a .h5ad file or a directory containing .h5ad files: {input_path}")

def process_query_dataset(h5ad_file, scpoli_model, source_adata, figures_folder, condition_key='orig.ident', cell_type_key='reanno'):
    """
    Process a single query dataset with scPoli model.
    
    Args:
        h5ad_file (str): Path to the h5ad file
        scpoli_model: Trained scPoli model
        source_adata: Reference dataset
        figures_folder (str): Output directory for figures and results
        condition_key (str): Column name for batch/condition information
        cell_type_key (str): Column name for cell type annotation
    """
    file_name = os.path.basename(h5ad_file)  # Extract just the file name
    print(f"Processing {file_name}...")

    try:
        # Load the query dataset
        query_adata = sc.read_h5ad(h5ad_file)
        
        # Check if condition_key column exists, if not, create it
        if condition_key not in query_adata.obs.columns:
            # Extract dataset name from filename (remove .h5ad extension)
            default_dataset_name = Path(file_name).stem
            query_adata.obs[condition_key] = default_dataset_name
            print(f"Added missing '{condition_key}' column with value: {default_dataset_name}")

        # 1. Save original query information
        original_query_cell_ids = query_adata.obs_names.copy()
        original_cell_count = query_adata.n_obs
        print(f"Original query cell count: {original_cell_count}")
        print(f"Original query cell ID examples: {original_query_cell_ids[:5].tolist()}")

        # Create a set of unique values
        unique_values = set(query_adata.obs[condition_key])
        if len(unique_values) == 1:
            query_name = unique_values.pop()
        else:
            raise ValueError(f"Expected only one unique value in '{condition_key}', but found multiple or none.")
        
        print(f"Query name: {query_name}")

        # 2. Check overlap in source_adata using the dynamic query_name
        source_overlap_mask = source_adata.obs[condition_key] == query_name
        source_overlap_count = source_overlap_mask.sum()
        source_overlap_cells = source_adata.obs_names[source_overlap_mask]
        print(f"{query_name} cells in source data: {source_overlap_count}")
        if source_overlap_count > 0:
            print(f"{query_name} cell ID examples in source data: {source_overlap_cells[:5].tolist()}")

        # 3. Check for identical cell IDs (unlikely but need to confirm)
        overlapping_ids = set(original_query_cell_ids) & set(source_overlap_cells)
        print(f"Number of identical cell IDs: {len(overlapping_ids)}")

        # Preprocess query dataset
        sc.settings.seed = 42
        query_adata.layers["counts"] = query_adata.X.copy()
        sc.pp.normalize_total(query_adata, target_sum=1e4)
        sc.pp.log1p(query_adata)
        query_adata.layers["logcounts"] = query_adata.X.copy()
        sc.pp.highly_variable_genes(query_adata, n_top_genes=2000, flavor="cell_ranger", batch_key=condition_key)
        sc.tl.pca(query_adata, n_comps=30, use_highly_variable=True)
        counts_matrix = query_adata.layers["counts"].toarray()

        query_adata = sc.AnnData(
            X=counts_matrix,
            obs=query_adata.obs.copy(),
            var=query_adata.var.copy(),
            layers={'counts': counts_matrix}
        )
        query_adata = remove_sparsity(query_adata)

        # Important: Add unique identifiers to query_adata
        query_adata.obs['is_original_query'] = True
        query_adata.obs['original_query_id'] = range(len(query_adata))

        # Reorganize query dataset to match genes in the reference dataset
        all_genes = source_adata.var_names
        missing_genes = all_genes.difference(query_adata.var_names)
        missing_data = np.zeros((query_adata.shape[0], len(missing_genes)))

        query_adata_df = pd.DataFrame(query_adata.X, columns=query_adata.var_names, index=query_adata.obs_names)
        missing_df = pd.DataFrame(missing_data, columns=missing_genes, index=query_adata.obs_names)
        query_adata_combined_df = pd.concat([query_adata_df, missing_df], axis=1)[all_genes]

        query_adata_extended = sc.AnnData(
            X=query_adata_combined_df.values,
            obs=query_adata.obs.copy(),  # This will preserve is_original_query and original_query_id
            var=pd.DataFrame(index=all_genes),
            layers={'counts': query_adata_combined_df.values}
        )

        # Check if 'features' column exists in the original query_adata.var before copying
        if 'features' in query_adata.var.columns:
            # If features column exists, reindex and copy it
            query_adata_extended.var['features'] = query_adata.var.reindex(all_genes)['features']
            print("Copied existing 'features' column from original query data")
        else:
            # If features column doesn't exist, check if source_adata has it and copy from there
            if 'features' in source_adata.var.columns:
                query_adata_extended.var['features'] = source_adata.var['features'].copy()
                print("Copied 'features' column from source_adata")
            else:
                # If neither has features column, create a default one using gene names
                query_adata_extended.var['features'] = query_adata_extended.var_names
                print("Created 'features' column using gene names as default")

        print(f"Extended query data shape: {query_adata_extended.shape}")
        print(f"is_original_query column exists: {'is_original_query' in query_adata_extended.obs.columns}")
        print(f"features column exists: {'features' in query_adata_extended.var.columns}")
        print(f"Features column sample: {query_adata_extended.var['features'].head().tolist()}")

        # Label transfer to query dataset
        query_adata_extended.obs[cell_type_key] = 'Unknown'
        scpoli_query = scPoli.load_query_data(
            adata=query_adata_extended,
            reference_model=scpoli_model,
            labeled_indices=[],
        )
        scpoli_query.train(
            n_epochs=50,
            pretraining_epochs=40,
            eta=5
        )
        query_adata_extended.X = query_adata_extended.X.astype(np.float32)

        # Label transfer from reference to query
        results_dict = scpoli_query.classify(query_adata_extended, scale_uncertainties=True)

        # Check the label transfer performance achieved
        for i in range(len(cell_type_key)):
            preds = results_dict[cell_type_key]["preds"]
            results_dict[cell_type_key]["uncert"]
            classification_df = pd.DataFrame(
                classification_report(
                    y_true=query_adata_extended.obs[cell_type_key],
                    y_pred=preds,
                    output_dict=True,
                )
            ).transpose()
        print(classification_df)

        # Get latent representation of reference data
        scpoli_query.model.eval()
        data_latent_source = scpoli_query.get_latent(source_adata, mean=True)
        adata_latent_source = sc.AnnData(data_latent_source)
        adata_latent_source.obs = source_adata.obs.copy()
        # Add identifiers for source data
        adata_latent_source.obs['is_original_query'] = False
        adata_latent_source.obs['original_query_id'] = -1  # Use -1 to mark non-query data

        # Get latent representation of query data
        data_latent = scpoli_query.get_latent(query_adata_extended, mean=True)
        adata_latent = sc.AnnData(data_latent)
        adata_latent.obs = query_adata_extended.obs.copy()

        adata_latent.obs[f'{cell_type_key}_pred'] = results_dict[cell_type_key]['preds'].tolist()
        adata_latent.obs[f'{cell_type_key}_uncert'] = results_dict[cell_type_key]['uncert'].tolist()
        adata_latent.obs['classifier_outcome'] = (
            adata_latent.obs[f'{cell_type_key}_pred'] == adata_latent.obs[cell_type_key]
        )

        # Get prototypes
        labeled_prototypes = scpoli_query.get_prototypes_info()
        labeled_prototypes.obs['study'] = 'labeled prototype'
        labeled_prototypes.obs['is_original_query'] = False
        labeled_prototypes.obs['original_query_id'] = -2  # Use -2 to mark labeled prototypes

        unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
        unlabeled_prototypes.obs['study'] = 'unlabeled prototype'
        unlabeled_prototypes.obs['is_original_query'] = False
        unlabeled_prototypes.obs['original_query_id'] = -3  # Use -3 to mark unlabeled prototypes

        # Join AnnDatas
        adata_latent_full = adata_latent_source.concatenate(
            [adata_latent, labeled_prototypes, unlabeled_prototypes],
            batch_key='query'
        )

        print(f"\n=== Post-concatenation data analysis ===")
        print(f"Total data size: {adata_latent_full.n_obs}")
        print(f"Query identifier distribution: {adata_latent_full.obs['query'].value_counts()}")
        print(f"is_original_query distribution: {adata_latent_full.obs['is_original_query'].value_counts()}")

        # Check distribution of query_name cells
        query_name_mask = adata_latent_full.obs[condition_key] == query_name
        print(f"Total {query_name} cells: {query_name_mask.sum()}")

        query_name_by_query = adata_latent_full.obs[query_name_mask]['query'].value_counts()
        print(f"{query_name} cells by query group: {query_name_by_query}")

        query_name_original_query = adata_latent_full.obs[query_name_mask]['is_original_query'].value_counts()
        print(f"True original query among {query_name} cells: {query_name_original_query}")

        adata_latent_full.obs[f'{cell_type_key}_pred'][adata_latent_full.obs['query'].isin(['0'])] = np.nan
        sc.pp.neighbors(adata_latent_full, n_neighbors=15)
        sc.tl.umap(adata_latent_full)

        # Get AnnData without prototypes
        adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]

        # Plot UMAP
        sc.pl.umap(
            adata_no_prototypes,
            color=f'{cell_type_key}_pred',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_{cell_type_key}.png'
        )
        sc.pl.umap(
            adata_no_prototypes,
            color=condition_key,
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_dataset.png'
        )
        sc.pl.umap(
            adata_no_prototypes,
            color=f'{cell_type_key}_uncert',
            show=False,
            frameon=False,
            cmap='magma',
            vmax=1,
            save=f'_{file_name}_scPoli_uncert.png'
        )

        sc.pp.neighbors(adata_latent)
        sc.tl.leiden(adata_latent)
        sc.tl.umap(adata_latent)
        sc.pl.umap(
            adata_latent,
            color=f'{cell_type_key}_pred',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_query.png'
        )

        # Fixed final filtering - get only true original query cells
        print(f"\n=== Precise filtering of original query cells ===")

        # Method 1: Use is_original_query identifier (most accurate)
        original_query_mask = adata_latent_full.obs['is_original_query'] == True
        full_latent_correct = adata_latent_full.obs[original_query_mask]
        print(f"Filtered by is_original_query: {len(full_latent_correct)} rows")

        # Method 2: Use original cell IDs (for validation)
        query_cells_mask = adata_latent_full.obs.index.isin(original_query_cell_ids)
        full_latent_by_ids = adata_latent_full.obs[query_cells_mask]
        print(f"Filtered by original cell IDs: {len(full_latent_by_ids)} rows")

        # Method 3: Use query identifier but exclude overlapping cells from source (old incorrect method)
        old_method_mask = adata_latent_full.obs[condition_key] == query_name
        full_latent_old_method = adata_latent_full.obs[old_method_mask]
        print(f"Old method (includes overlapping data): {len(full_latent_old_method)} rows")

        # Validate results
        if len(full_latent_correct) == original_cell_count:
            full_latent = full_latent_correct
            print(f"✅ Success! Got correct cell count: {len(full_latent)}")
            print(f"Number of NaN values: {full_latent.isna().sum().sum()}")
        elif len(full_latent_by_ids) == original_cell_count:
            full_latent = full_latent_by_ids
            print(f"✅ Got correct result through cell IDs: {len(full_latent)}")
        else:
            print(f"❌ Still have issues, need further debugging")
            print(f"Expected: {original_cell_count}, Actual: {len(full_latent_correct)}")
            # Use the best available result
            full_latent = full_latent_correct

        # Additional validation: check if full_latent still contains overlapping cells from source
        if 'full_latent' in locals():
            source_contamination = full_latent['original_query_id'] == -1
            if source_contamination.any():
                print(f"⚠️  Still have {source_contamination.sum()} cells from source")
            else:
                print(f"✅ Successfully excluded all overlapping cells from source")

        print(f"\nFinal results summary:")
        print(f"- Original query cells: {original_cell_count}")
        print(f"- Final obtained cells: {len(full_latent) if 'full_latent' in locals() else 'N/A'}")
        print(f"- Overlapping cells in source: {source_overlap_count}")
        print(f"- Problem solved: {'✅' if 'full_latent' in locals() and len(full_latent) == original_cell_count else '❌'}")

        # Transfer predictions to original adata
        print("Transferring annotation predictions to original adata")
        matching_indices = adata_latent.obs.index.intersection(query_adata.obs.index)
        if f'{cell_type_key}_pred' not in query_adata.obs.columns:
            query_adata.obs[f'{cell_type_key}_pred']  = np.nan
        if f'{cell_type_key}_uncert' not in query_adata.obs.columns:
            query_adata.obs[f'{cell_type_key}_uncert'] = np.nan
            
        query_adata.obs.loc[matching_indices, f'{cell_type_key}_pred'] = adata_latent.obs.loc[matching_indices, f'{cell_type_key}_pred']
        query_adata.obs.loc[matching_indices, f'{cell_type_key}_uncert'] = adata_latent.obs.loc[matching_indices, f'{cell_type_key}_uncert']
        
        # Store UMAP coordinates
        print("Storing annotation UMAP coordinates")
        query_adata.obsm[f'X_umap_{cell_type_key}'] = adata_latent.obsm['X_umap']
        
        # Save UMAP visualizations
        print("Saving UMAP visualizations")
        sc.pl.umap(adata_latent, color=[f'{cell_type_key}_pred', f'{cell_type_key}_uncert'], 
                  save=f'_{file_name}_{cell_type_key}.pdf', show=False)
        
        # Save the annotated dataset
        #output_file = os.path.join(figures_folder, f"{file_name}_scPoli_annotated.h5ad")
        #print(f"Saving annotated dataset to {output_file}")
        print(f"Saving annotated dataset")
        query_adata.write_h5ad(filename=h5ad_file)

        # Save results - only original query cells, excluding source_adata information
        if 'full_latent' in locals():
            output_csv = os.path.join(figures_folder, f'{file_name}_{cell_type_key}_scPoli_query.csv')
            # Save the .obs DataFrame (metadata) to CSV
            full_latent.to_csv(output_csv, index=True)
            print(f"Saved clean query results to {output_csv}")
            print(f"Saved {len(full_latent)} rows (original query cells only, no source contamination)")
            
            
            # Optional: Save a summary of what was excluded
            summary_info = {
                'file_name': file_name,
                'query_name': query_name,
                'original_query_cells': original_cell_count,
                'source_overlapping_cells': source_overlap_count,
                'final_saved_cells': len(full_latent),
                'cells_excluded': source_overlap_count,
                'data_clean': len(full_latent) == original_cell_count
            }
            
            summary_csv = os.path.join(figures_folder, f'{file_name}_processing_summary.csv')
            pd.DataFrame([summary_info]).to_csv(summary_csv, index=False)
            print(f"Saved processing summary to {summary_csv}")
        else:
            print("❌ Cannot save results - full_latent not properly generated")

        print(f"Completed processing {file_name}\n" + "="*50)

    except Exception as e:
        print(f"Error processing {file_name}: {e}")
        import traceback
        print(traceback.format_exc())

def main():
    """
    Main function to run scPoli query processing.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process query datasets with scPoli model')
    parser.add_argument('--query_path', '-q', type=str, required=True,
                       help='Path to query data (single .h5ad file or directory containing .h5ad files)')
    parser.add_argument('--reference_model', '-r', type=str, required=True,
                       help='Path to trained scPoli reference model directory')
    parser.add_argument('--reference_data', '-d', type=str, required=True,
                       help='Path to reference dataset (.h5ad file)')
    parser.add_argument('--output_dir', '-o', type=str, default='./figures',
                       help='Output directory for figures and results (default: ./figures)')
    parser.add_argument('--condition_key', '-c', type=str, default='orig.ident',
                       help='Column name for batch/condition information (default: orig.ident)')
    parser.add_argument('--cell_type_key', '-t', type=str, default='reanno',
                       help='Column name for cell type annotation (default: reanno)')
    parser.add_argument('--working_dir', '-w', type=str, default=None,
                       help='Working directory to change to (optional)')
    
    args = parser.parse_args()
    
    # Change working directory if specified
    if args.working_dir:
        os.chdir(args.working_dir)
        print(f"Changed working directory to: {os.getcwd()}")
    else:
        print(f"Current working directory: {os.getcwd()}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory: {args.output_dir}")
    
    # Load reference model and data
    print("Loading reference model and data...")
    try:
        source_adata = sc.read(args.reference_data)
        print(f"Reference data loaded: {source_adata.shape}")
        
        scpoli_model = scPoli.load(args.reference_model)
        print("Reference model loaded successfully")
    except Exception as e:
        print(f"Error loading reference model or data: {e}")
        return
    
    # Get list of h5ad files to process
    try:
        h5ad_files = get_h5ad_files(args.query_path)
        print(f"Found {len(h5ad_files)} .h5ad file(s) to process:")
        for f in h5ad_files:
            print(f"  - {f}")
    except Exception as e:
        print(f"Error getting h5ad files: {e}")
        return
    
    # Process each query dataset
    for h5ad_file in h5ad_files:
        process_query_dataset(
            h5ad_file=h5ad_file,
            scpoli_model=scpoli_model,
            source_adata=source_adata,
            figures_folder=args.output_dir,
            condition_key=args.condition_key,
            cell_type_key=args.cell_type_key
        )
    
    print(f"\nProcessing completed! Results saved in: {args.output_dir}")

def run_with_default_params():
    """
    Run with default hardcoded parameters (original behavior).
    """
    print("Running with default parameters (hardcoded paths)...")
    
    # Change the working directory to the Garfield folder (if needed)
    os.chdir('/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250729_scpoli_optimization_v3/reanno_batch_enhanced_pred_2ndround')
    print(f"Working directory: {os.getcwd()}")
    
    # Load trained ref model and data
    source_adata = sc.read('/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250729_scpoli_optimization_v3/enhanced_reference_model_reanno_2ndround/adata.h5ad')
    print(f"Reference data loaded: {source_adata.shape}")
    
    # Train reference model for "lineage" (only once)
    condition_key = 'orig.ident'
    cell_type_key = "reanno"
    early_stopping_kwargs = {
        "early_stopping_metric": "val_prototype_loss",
        "mode": "min",
        "threshold": 0,
        "patience": 20,
        "reduce_lr": True,
        "lr_patience": 13,
        "lr_factor": 0.1,
    }
    
    reference_model_dir = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250729_scpoli_optimization_v3/enhanced_reference_model_reanno_2ndround/'
    scpoli_model = scPoli.load(reference_model_dir)
    print("Reference model loaded successfully")
    
    # Process query datasets from folder
    query_folder = '/storage2/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/garfield/in_vitro_embryo_models/processed/data'
    figures_folder = './figures'
    os.makedirs(figures_folder, exist_ok=True)
    
    # Get h5ad files (modified to be more flexible)
    try:
        h5ad_files = get_h5ad_files(query_folder)
        print(f"Found {len(h5ad_files)} .h5ad file(s) to process")
    except Exception as e:
        print(f"Error: {e}")
        return
    
    # Process each query dataset
    for h5ad_file in h5ad_files:
        process_query_dataset(
            h5ad_file=h5ad_file,
            scpoli_model=scpoli_model,
            source_adata=source_adata,
            figures_folder=figures_folder,
            condition_key=condition_key,
            cell_type_key=cell_type_key
        )
    
    print(f"\nProcessing completed! Results saved in: {figures_folder}")

if __name__ == "__main__":
    # If running without command line arguments, use the original hardcoded paths
    if len(os.sys.argv) == 1:
        run_with_default_params()
    else:
        # Run with command line arguments
        main()
        
# Process a single h5ad file
#python script.py -q /path/to/single_file.h5ad -r /path/to/model -d /path/to/reference.h5ad -o ./output

# Process all h5ad files in a directory
#python script.py -q /path/to/query_folder -r /path/to/model -d /path/to/reference.h5ad -o ./output

# With custom parameters
#python script.py -q /path/to/data -r /path/to/model -d /path/to/reference.h5ad -o ./output -c "batch" -t "cell_type"