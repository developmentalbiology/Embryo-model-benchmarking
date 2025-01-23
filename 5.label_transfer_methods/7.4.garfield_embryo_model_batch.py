#!/usr/bin/env python
# coding: utf-8

import sys
import os

# Add the Garfield directory to the Python path
sys.path.append('/storage/liuxiaodongLab/fanxueying/mayanalysis/Garfield_run')

# Change the working directory to the Garfield folder (if needed)
os.chdir('/storage/liuxiaodongLab/fanxueying/mayanalysis/Garfield_run')
os.getcwd()

# Import Garfield
import Garfield as gf

# load packages
import os
import warnings
import pandas as pd
import scanpy as sc
from mudata import MuData
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

gf.__version__


# Load dataset
os.chdir("/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug")
adata = sc.read('human_reanno_20250108.h5ad')

print(adata)
print(adata.layers["counts"])

adata.X = adata.layers["counts"].copy()


# Set working directory for Garfield
os.chdir('/storage/liuxiaodongLab/fanxueying/mayanalysis/Garfield_run/20250118_Garfield_embryo_model_batch')


# Copy data for reference
ref_adata = adata.copy()
ref_adata.X = ref_adata.layers['counts'].copy()

# Inspect the batches contained in the dataset
print(ref_adata.obs['orig.ident'].value_counts())
print(ref_adata.X.max())

# Configure Garfield settings
workdir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/Garfield_run/20250118_Garfield_embryo_model_batch'
gf.settings.set_workdir(workdir)

### modify parameter
user_config = dict(
    ## Input options
    data_dir=workdir,  # STR     Location of the dataset to be used.         Default is `data`.
    project_name='human_embryo_model',  # STR     Name of the dataset to be used.             Default is `name`.
    adata_list=ref_adata,  # STR     adata object of single-cell dataset.         Default is `adata`.
    profile='RNA',   # STR     Type of single-cell dataset.                  Default is `RNA`.
    sample_col='orig.ident',  # STR     Column name of sample in adata.obs.       Default is `batch`.

    ## Preprocessing options
    # used_feat=False,
    rna_n_top_features=3000,
    metric='euclidean',  # STR     Metric for clustering.                   Default is `euclidean`.

    ## Model options
    conv_type='GATv2Conv', # GAT or GATv2Conv or GCN
    gnn_layer=2,
    hidden_dims=[128, 128],
    svd_q=5,  # default=5, type=int, help='rank'
    cluster_num=20, # default=20, type=int, help='number of clusters for contrastive loss'
    test_split=0.2,
    val_split=0.2,
    used_edge_weight=False,
    used_recon_exp=True,
    used_DSBN=False,
    used_mmd=True,
    loader_type='graphsaint',
    batch_size=128,  # INT   batch size of model training
    num_neighbors=[5, 5],
    epochs=100,  # INT       Number of epochs.                        Default is 100.
    mmd_temperature=0.8,  ## mmd regu
    instance_temperature=1.0,
    cluster_temperature=0.5,
    monitor_only_val_losses=False,
    learning_rate=0.0001
)

dict_config = gf.settings.set_gf_params(user_config)

# Start training
from Garfield.model import GarfieldTrainer
trainer = GarfieldTrainer(dict_config)
trainer.fit()


# Save the model
ref_path = '20250118_Garfield_embryo_model_batch/human_reference_model/'
trainer.save(ref_path, overwrite=True)

# Get latent representation
adata_final = trainer.get_latent_representation()
# Plot and save reference UMAP
sc.pp.neighbors(adata_final, use_rep='X_gf')
sc.tl.umap(adata_final)
sc.pl.umap(adata_final, color=['orig.ident', 'final_anno'], wspace=0.20, save='_reference_umap.png')

# Process query datasets from folder
query_folder = '/storage/liuxiaodongLab/fanxueying/mayanalysis/2024Aug/garfield/in_vitro_embryo_models/processed/data'
#figures_folder = os.path.join(query_folder, 'figures')

# Create figures directory if it doesn't exist
#os.makedirs(figures_folder, exist_ok=True)

# List all h5ad files in the directory
h5ad_files = [os.path.join(query_folder, f) for f in os.listdir(query_folder) if f.startswith('corrected_processed_')]

for h5ad_file in h5ad_files:
    file_name = os.path.basename(h5ad_file)  # Extract just the file name
    print(f"Processing {file_name}...")  # Print only the file name

   
    try:
        # Load the data
        query_adata = sc.read_h5ad(h5ad_file)
        query_adata.layers["counts"] = query_adata.X.copy()

        query_adata
        query_adata.X = query_adata.layers['counts'].copy()
        
        
        # Load query data and perform label transfer
        new_model = trainer.load_query_data(adata=query_adata, reference_model=ref_path, project_name=f"{file_name}_query")
        new_model.fit()
        
        # Get latent representation for the query
        adata_concat = new_model.get_latent_representation()
        adata_concat


        # Plot UMAP
        sc.pp.neighbors(adata_concat, use_rep='X_gf')
        sc.tl.umap(adata_concat)
        sc.pl.umap(adata_concat, color=['projection', 'orig.ident'], wspace=0.20, edges=False, save=f'{file_name}_garfield_query.png')

        # Separate reference and query data
        adata_ref = adata_concat[adata_concat.obs['projection'] == 'reference', :]
        adata_query = adata_concat[adata_concat.obs['projection'] == 'query', :]

        # Perform label transfer
        adata_query = new_model.label_transfer(ref_adata=adata_ref,
                                                ref_adata_emb='X_gf',
                                                query_adata=adata_query,
                                                query_adata_emb='X_gf',
                                                n_neighbors=5,
                                                ref_adata_obs=adata_ref.obs,
                                                label_keys='final_anno')
        
        adata_query = new_model.label_transfer(ref_adata=adata_ref,
                                        ref_adata_emb='X_gf',
                                        query_adata=adata_query,
                                        query_adata_emb='X_gf',
                                        n_neighbors=5,
                                        ref_adata_obs=adata_ref.obs,
                                        label_keys='final_lineage')
 
        
        # Export the results
        obs_df = pd.DataFrame(adata_query.obs)
        output_file = f'{file_name}_garfield_query.obs.csv'
        obs_df.to_csv(output_file, index=True)

        print(f"Finished saving csv {file_name}")

        ## predicted label
        sc.pl.umap(adata_query,color=['transferred_lineage_unfiltered', 'transferred_reanno_unfiltered'], save=f'{file_name}__transfer.png')

        # Save the query AnnData object
        #adata_query.raw.var.rename(columns={'_index': 'index'}, inplace=True)

        #output_filename = f"{file_name}_garfield_transfer"
        #adata_query.write_h5ad(filename=output_filename)


        print(f"Finished processing {file_name}")

    except Exception as e:
        print(f"An error occurred while processing {file_name}: {e}")
        print(f"Skipping {file_name}.")

print("All query datasets have been processed.")
