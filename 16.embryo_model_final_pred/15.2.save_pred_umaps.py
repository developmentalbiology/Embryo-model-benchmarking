#!/usr/bin/env python
# coding: utf-8

import os
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Set up scanpy settings
sc.settings.set_figure_params(dpi=300, frameon=False)

# Define lineage colors
lineage_color = {
    "Amniotic_ecto": "#1f77b4",    # Blue
    "Notochord": "#aa40fc",        # Purple
    "Endoderm": "#ff7f0e",         # Orange
    "PGC": "#8c564b",              # Brown
    "ExE_endo": "#279e68",         # Green
    "Primitive.streak": "#e377c2", # Pink
    "NMP": "#d62728",              # Red
    "TE_TrB": "#b5bd61",           # Yellow-Green
    "epi": "#17becf",              # Cyan
    "hemogenic": "#aec7e8",        # Light Blue
    "meso_Exe.meso": "#ffbb78",    # Orange-Yellow
    "neural_ecto": "#98df8a"       # Light Green
}

# No need for custom reanno colors - we'll use scanpy defaults with consistent ordering

# Create an ordered list of lineage labels based on the color dictionary keys
lineage_ordered_labels = list(lineage_color.keys())

# Your ordered labels list for annotation visualization
ordered_labels = [
    'TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
    'Epiblast_1','Epiblast_2','Epiblast_3','Ectoderm',
    'Amniontic.epi','Amniontic.ectoderm',
    'PGC',
    'Primitive.streak',
    'Neuromesodermal.progenitor',
    'Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 'Neural.ectoderm.midbrain','Spinal.cord',
    'Paraxial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 'Rostral.mesoderm', 'Lateral.plate.mesoderm_1',
    'Lateral.plate.mesoderm_2','Lateral.plate.mesoderm_3','Cardiac.mesoderm','Amniotic.mesoderm','Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2',
    'Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm',
    'DE','Gut',
    'Notochord',
    'Hemogenic.endothelial.progenitor','Endothelium','Erythroid','Primitive.megakaryocyte','Myeloid.progenitor'
]

# Define input and output directories
input_folder = '/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output'
figure_output_folder = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250801_embryomodel_export_plots/figures'

# Create figure output directory if it doesn't exist
os.makedirs(figure_output_folder, exist_ok=True)

# List all h5ad files in the processed directory
h5ad_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('_scPoli_annotated.h5ad')]

print(f"Found {len(h5ad_files)} files to process:")
for file in h5ad_files:
    print(f"  - {os.path.basename(file)}")

# Process each file
for file_path in h5ad_files:
    try:
        # Load the data
        print(f"\nLoading: {os.path.basename(file_path)}")
        adata = sc.read_h5ad(file_path)
        
        # Extract dataset name
        dataset_name = os.path.basename(file_path).replace('_with_scpoli.h5ad', '').replace('corrected_processed_', '')
        
        print(f"Processing {dataset_name}...")
    
        
        # Create figure directory for this dataset
        dataset_figure_folder = os.path.join(figure_output_folder, dataset_name)
        os.makedirs(dataset_figure_folder, exist_ok=True)
        
        # Plot lineage prediction
        if 'lineage_pred' in adata.obs.columns:
            print(f"  Creating lineage UMAP plot...")
            
            adata.obsm['X_umap'] = adata.obsm['X_umap_anno']
            
            # Set up categorical data with proper ordering and colors
            adata.obs['lineage_pred'] = adata.obs['lineage_pred'].astype('category')
            present_lineages = adata.obs['lineage_pred'].cat.categories.tolist()
            ordered_present_lineages = [l for l in lineage_ordered_labels if l in present_lineages]
            ordered_present_lineages.extend([l for l in present_lineages if l not in lineage_ordered_labels])
            adata.obs['lineage_pred'] = adata.obs['lineage_pred'].cat.reorder_categories(ordered_present_lineages)
            
            # Create color palette for present lineages
            lineage_colors_present = [lineage_color.get(lineage, '#cccccc') for lineage in ordered_present_lineages]
            
            sc.pl.umap(
                adata,
                color='lineage_pred',
                palette=lineage_colors_present,
                show=False,
                frameon=False,
                size=25,    #size=25, for zheng_2019 use 40
                save=f'_{dataset_name}_lineage.pdf'
            )
            # Move the saved file to our organized folder
            source_file = f'umap_{dataset_name}_lineage.pdf'
            if os.path.exists(source_file):
                os.rename(source_file, os.path.join(dataset_figure_folder, f'{dataset_name}_umap_lineage.pdf'))
                
            # Plot lineage uncert
            sc.pl.umap(
                adata,
                color='lineage_uncert',
                show=False,
                frameon=False,
                cmap='magma',
                vmax=1,
                size=25,
                save=f'_{dataset_name}_lineage_uncert.pdf'
            )
            # Move the saved file to our organized folder
            source_file = f'umap_{dataset_name}_lineage_uncert.pdf'
            if os.path.exists(source_file):
                os.rename(source_file, os.path.join(dataset_figure_folder, f'{dataset_name}_umap_lineage_uncert.pdf'))
        
        # Plot reannotation prediction
        if 'reanno_pred' in adata.obs.columns:
            print(f"  Creating reannotation UMAP plot...")
            
            adata.obsm['X_umap'] = adata.obsm['X_umap_anno']
            
            # Set up categorical data with consistent ordering across all datasets
            adata.obs['reanno_pred'] = adata.obs['reanno_pred'].astype('category')
            present_reanno = adata.obs['reanno_pred'].cat.categories.tolist()
            
            # Always use the same order: first the ordered_labels that are present, then any extras
            ordered_present_reanno = [ra for ra in ordered_labels if ra in present_reanno]
            ordered_present_reanno.extend([ra for ra in present_reanno if ra not in ordered_labels])
            
            # Add all missing categories from ordered_labels to maintain consistent color mapping
            # This ensures that each category always gets the same default scanpy color
            all_categories_ordered = ordered_labels.copy()
            for extra in present_reanno:
                if extra not in all_categories_ordered:
                    all_categories_ordered.append(extra)
            
            # Set categories to the full ordered list (even if some are missing from data)
            adata.obs['reanno_pred'] = adata.obs['reanno_pred'].cat.set_categories(all_categories_ordered)
            
            sc.pl.umap(
                adata,
                color='reanno_pred',
                show=False,
                frameon=False,
                size=25,
                save=f'_{dataset_name}_reanno.pdf'
            )
            # Move the saved file to our organized folder
            source_file = f'umap_{dataset_name}_reanno.pdf'
            if os.path.exists(source_file):
                os.rename(source_file, os.path.join(dataset_figure_folder, f'{dataset_name}_umap_reanno.pdf'))
                
        
            # Plot reanno uncert
            sc.pl.umap(
                adata,
                color='reanno_uncert',
                show=False,
                frameon=False,
                cmap='magma',
                vmax=1,
                size=25,
                save=f'_{dataset_name}_reanno_uncert.pdf'
            )
            # Move the saved file to our organized folder
            source_file = f'umap_{dataset_name}_reanno_uncert.pdf'
            if os.path.exists(source_file):
                os.rename(source_file, os.path.join(dataset_figure_folder, f'{dataset_name}_umap_reanno_uncert.pdf'))
                
            # Export expression of marker genes
            print(f"  Creating gene expression UMAP plots...")
            
            # Check if data is log-transformed, if not, normalize and log-transform
            # Check if data is log-transformed, if not, normalize and log-transform
            # Check if X contains negative values (indicating already log-transformed) or if max is very high (raw counts)
            import numpy as np
            from scipy import sparse
            
            # Handle both sparse and dense matrices
            if sparse.issparse(adata.X):
                max_val = adata.X.max()
                has_negative = (adata.X < 0).nnz > 0  # Check if any values are negative in sparse matrix
            else:
                max_val = np.max(adata.X)
                has_negative = np.any(adata.X < 0)
            
            print(f"    Data check - Max value: {max_val:.2f}, Has negative values: {has_negative}")
            
            # If max value is very high (>50) and no negative values, likely raw counts
            if max_val > 50 and not has_negative:
                print(f"    Data appears to be raw counts, applying normalization and log transformation...")
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
                print(f"    Normalization and log transformation completed")
            else:
                print(f"    Data appears to already be normalized/log-transformed")
            
            # Define list of genes
            import seaborn as sns
            
            genes_to_plot = [
                "POU5F1", "NANOG",  # epiblast
                "SOX2", "TTYH1",    # neural ectoderm
                "GATA3", "TFAP2A",
                "TBXT", "CDX1", "PDGFRA", "APOA2", "FOXA2", "NANOS3",
                "PECAM1", "HBZ", "PTPRC", 
                "GABRP",  "COL6A1", "COL6A2"
            ]

            # Create gene expression directory for this dataset
            gene_output_dir = os.path.join(dataset_figure_folder, "gene_expression")
            os.makedirs(gene_output_dir, exist_ok=True)

            # Loop over each gene and plot/save individually
            for gene in genes_to_plot:
                # 1. Create a shuffled index array
                shuffled_indices = np.random.permutation(adata.n_obs) # or adata.shape[0]

                # 2. Create a temporary AnnData object with shuffled cells
                #    This avoids modifying the original adata object
                adata_shuffled = adata[shuffled_indices, :].copy()

                # 3. Plot using the shuffled AnnData object
                fig, ax = plt.subplots(figsize=(3, 3))  # Create a new figure for each gene
                sc.pl.umap(
                    adata_shuffled, # Use the shuffled AnnData
                    color=gene,
                    use_raw=False,
                    cmap=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True),
                    ax=ax,
                    show=False,  # Prevents Scanpy from calling plt.show()
                    size=15
                )
                plt.tight_layout()
                plt.savefig(os.path.join(gene_output_dir, f"{dataset_name}_{gene}_umap.pdf"), dpi=300)
                plt.close()
                # Optional: Delete the temporary object to free memory if processing many genes
                # del adata_shuffled
        
        # Overwrite the h5ad file with updated data (not needed since data is already processed)
        # print(f"  Saving updated h5ad file...")
        # adata.write_h5ad(file_path)
        
        print(f"  Successfully processed {dataset_name}")
        print(f"    Generated plots: lineage_pred, lineage_uncert, reanno_pred, reanno_uncert")
        
    except Exception as e:
        print(f"  Error processing {file_path}: {str(e)}")
        continue

print(f"\nProcessing complete!")
print(f"UMAP plots saved to: {figure_output_folder}")
print(f"Generated 4 plots per dataset: lineage_pred, lineage_uncert, reanno_pred, reanno_uncert")