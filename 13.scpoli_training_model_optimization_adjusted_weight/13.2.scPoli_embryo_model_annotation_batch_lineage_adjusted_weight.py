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

import warnings
# Suppress warnings
warnings.filterwarnings('ignore')


sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)

# Change the working directory to the Garfield folder (if needed)
os.chdir('/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250731_scpoli_optimization_comparasion_v3/ref_training_balance_weight')
os.getcwd()

# Load trained ref model and data
source_adata = sc.read('/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250731_scpoli_optimization_comparasion_v3/ref_training_balance_weight/method3_balanced_weighted_lineage/adata.h5ad')
source_adata


# Train reference model for "lineage" (only once)
condition_key = 'orig.ident'
cell_type_key = "lineage"
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

reference_model_dir = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250731_scpoli_optimization_comparasion_v3/ref_training_balance_weight/method3_balanced_weighted_lineage/'
scpoli_model = scPoli.load(reference_model_dir)

# Process query datasets from folder
query_folder = '/storage2/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/garfield/in_vitro_embryo_models/processed/data'
figures_folder = './method3_balanced_weighted_lineage/figures'
os.makedirs(figures_folder, exist_ok=True)  # Create the figures directory if it doesn't exist

# List all h5ad files in the directory
h5ad_files = [os.path.join(query_folder, f) for f in os.listdir(query_folder) if f.startswith('corrected_processed_')]

for h5ad_file in h5ad_files:
    file_name = os.path.basename(h5ad_file)  # Extract just the file name
    print(f"Processing {file_name}...")  # Print only the file name

    try:
        # Load the query dataset
        query_adata = sc.read_h5ad(h5ad_file)

        # Create a set of unique values
        unique_values = set(query_adata.obs["orig.ident"])
        # Ensure there is only one unique value
        if len(unique_values) == 1:
            query_name = unique_values.pop()  # Extract the single value from the set
        else:
            raise ValueError("Expected only one unique value in 'orig.ident', but found multiple or none.")
            
        print(f"Query name: {query_name}")

        # Preprocess query dataset
        sc.settings.seed = 42
        query_adata.layers["counts"] = query_adata.X.copy()
        sc.pp.normalize_total(query_adata, target_sum=1e4)
        sc.pp.log1p(query_adata)
        query_adata.layers["logcounts"] = query_adata.X.copy()
        sc.pp.highly_variable_genes(query_adata, n_top_genes=2000, flavor="cell_ranger", batch_key="orig.ident")
        sc.tl.pca(query_adata, n_comps=30, use_highly_variable=True)
        counts_matrix = query_adata.layers["counts"].toarray()

        query_adata = sc.AnnData(
            X=counts_matrix,
            obs=query_adata.obs.copy(),
            var=query_adata.var.copy(),
            layers={'counts': counts_matrix}
        )
        query_adata = remove_sparsity(query_adata)

        # Reorganize query dataset to match genes in the reference dataset
        all_genes = source_adata.var_names
        missing_genes = all_genes.difference(query_adata.var_names)
        missing_data = np.zeros((query_adata.shape[0], len(missing_genes)))
        query_adata_df = pd.DataFrame(query_adata.X, columns=query_adata.var_names, index=query_adata.obs_names)
        missing_df = pd.DataFrame(missing_data, columns=missing_genes, index=query_adata.obs_names)
        query_adata_combined_df = pd.concat([query_adata_df, missing_df], axis=1)[all_genes]
        query_adata_extended = sc.AnnData(
            X=query_adata_combined_df.values,
            obs=query_adata.obs,
            var=pd.DataFrame(index=all_genes),
            layers={'counts': query_adata_combined_df.values}
        )
        query_adata_extended.var['features'] = query_adata.var.reindex(all_genes)['features']

        # Label transfer to query dataset
        query_adata_extended.obs['lineage'] = 'Unknown'
        scpoli_query = scPoli.load_query_data(
            adata=query_adata_extended,
            reference_model=scpoli_model,
            labeled_indices=[],
        )
        scpoli_query.train(
            n_epochs=100,
            pretraining_epochs=40,
            eta=10
        )
        query_adata_extended.X = query_adata_extended.X.astype(np.float32)

        #Label transfer from reference to query
        results_dict = scpoli_query.classify(query_adata_extended, scale_uncertainties=True)

        # check the label transfer performance achieved
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

        #get latent representation of reference data
        scpoli_query.model.eval()
        data_latent_source = scpoli_query.get_latent(
            source_adata,
            mean=True
        )
        adata_latent_source = sc.AnnData(data_latent_source)
        adata_latent_source.obs = source_adata.obs.copy()

        # Get latent representation of query data
        data_latent = scpoli_query.get_latent(query_adata_extended, mean=True)
        adata_latent = sc.AnnData(data_latent)
        adata_latent.obs = query_adata_extended.obs.copy()

        adata_latent.obs['lineage_pred'] = results_dict['lineage']['preds'].tolist()
        adata_latent.obs['lineage_uncert'] = results_dict['lineage']['uncert'].tolist()
        adata_latent.obs['classifier_outcome'] = (
            adata_latent.obs['lineage_pred'] == adata_latent.obs['lineage']
        )

        # Get prototypes
        labeled_prototypes = scpoli_query.get_prototypes_info()
        labeled_prototypes.obs['study'] = 'labeled prototype'
        unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
        unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

        # Join AnnDatas
        adata_latent_full = adata_latent_source.concatenate(
            [adata_latent, labeled_prototypes, unlabeled_prototypes],
            batch_key='query'
        )
        adata_latent_full.obs['lineage_pred'][adata_latent_full.obs['query'].isin(['0'])] = np.nan
        sc.pp.neighbors(adata_latent_full, n_neighbors=15)
        sc.tl.umap(adata_latent_full)

        # Get AnnData without prototypes
        adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]

        # Plot UMAP
        sc.pl.umap(
            adata_no_prototypes,
            color='lineage_pred',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_lineage.png'
        )
        sc.pl.umap(
            adata_no_prototypes,
            color='orig.ident',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_dataset.png'
        )
        sc.pl.umap(
            adata_no_prototypes,
            color='lineage_uncert',
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
            color='lineage_pred',
            show=False,
            frameon=False,
            save=f'_{file_name}_scPoli_query.png'
        )

        adata_latent_full
        full_latent = adata_latent_full.obs[adata_latent_full.obs["orig.ident"] == query_name]
        full_latent

        # Save results
        output_csv = os.path.join(figures_folder, f'{file_name}_scPoli_query.csv')
        # Save the .obs DataFrame (metadata) to CSV
        full_latent.to_csv(output_csv, index=True)
        print(f"Saved results to {output_csv}")

    except Exception as e:
        print(f"Error processing {file_name}: {e}")