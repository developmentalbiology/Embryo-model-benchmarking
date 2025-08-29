import scanpy as sc
import pandas as pd
import numpy as np
from scib import metrics
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import scPoli
try:
    from scarches.models.scpoli import scPoli
except ImportError:
    print("scPoli not found. Please install scPoli first.")
    print("pip install scpoli")

def extract_latent_and_compute_umap(model_dir, adata_path, use_rep_key="X_scpoli"):
    """
    Extract latent representation using scPoli model and compute UMAP
    
    Parameters:
    -----------
    model_dir : str
        Directory containing the scPoli model
    adata_path : str
        Path to the adata.h5ad file
    use_rep_key : str
        Key to store latent representation
    
    Returns:
    --------
    adata : AnnData
        AnnData object with latent representation and UMAP
    """
    
    print(f"Loading scPoli model from: {model_dir}")
    scpoli_model = scPoli.load(model_dir)
    
    print(f"Loading data from: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    
    print("Extracting latent representation...")
    # Set model to evaluation mode
    scpoli_model.model.eval()
    
    # Get latent representation
    data_latent = scpoli_model.get_latent(adata, mean=True)
    adata.obsm[use_rep_key] = data_latent
    
    print("Computing neighbors and UMAP...")
    # Compute neighbors using latent representation
    sc.pp.neighbors(adata, use_rep=use_rep_key)
    sc.tl.umap(adata)
    
    print(f"Latent representation shape: {data_latent.shape}")
    print(f"Available obsm keys: {list(adata.obsm.keys())}")
    
    return adata

def evaluate_scpoli_clustering(adata, 
                             clustering_key,
                             reference_clustering_key=None,
                             batch_key='orig.ident',
                             embedding_key='X_umap',
                             latent_key='X_scpoli'):
    """
    Evaluate scPoli clustering results using both latent space and UMAP
    """
    
    results = {}
    
    print(f"Evaluating clustering: {clustering_key}")
    
    # Check if clustering results exist
    if clustering_key not in adata.obs.columns:
        raise ValueError(f"Clustering key '{clustering_key}' not found in adata.obs")
    
    try:
        # 1. Silhouette Score on latent space (more meaningful for scPoli)
        print("Computing silhouette score on latent space...")
        sil_latent = metrics.silhouette(
            adata,
            group_key=clustering_key,
            embed=latent_key
        )
        results['silhouette_latent'] = sil_latent
        
        # 2. Silhouette Score on UMAP (for visualization comparison)
        print("Computing silhouette score on UMAP...")
        sil_umap = metrics.silhouette(
            adata,
            group_key=clustering_key,
            embed=embedding_key
        )
        results['silhouette_umap'] = sil_umap
        
    except Exception as e:
        print(f"Error in silhouette calculation: {e}")
    
    # 3. NMI and ARI if reference clustering is provided
    if reference_clustering_key and reference_clustering_key in adata.obs.columns:
        try:
            nmi_score = metrics.nmi(adata, group1=clustering_key, group2=reference_clustering_key)
            ari_score = metrics.ari(adata, group1=clustering_key, group2=reference_clustering_key)
            
            results['NMI_vs_reference'] = nmi_score
            results['ARI_vs_reference'] = ari_score
            
            print(f"NMI vs {reference_clustering_key}: {nmi_score:.4f}")
            print(f"ARI vs {reference_clustering_key}: {ari_score:.4f}")
            
        except Exception as e:
            print(f"Error in NMI/ARI calculation: {e}")
    
    # 4. Clustering statistics
    try:
        n_clusters = len(adata.obs[clustering_key].unique())
        cluster_sizes = adata.obs[clustering_key].value_counts()
        
        results['n_clusters'] = n_clusters
        results['min_cluster_size'] = cluster_sizes.min()
        results['max_cluster_size'] = cluster_sizes.max()
        results['mean_cluster_size'] = cluster_sizes.mean()
        results['cluster_size_std'] = cluster_sizes.std()
        
        print(f"Number of clusters: {n_clusters}")
        print(f"Cluster size range: {cluster_sizes.min()} - {cluster_sizes.max()}")
        
    except Exception as e:
        print(f"Error in clustering statistics: {e}")
    
    # 5. Batch effect evaluation
    if batch_key in adata.obs.columns:
        try:
            # Graph connectivity
            graph_conn = metrics.graph_connectivity(adata, label_key=batch_key)
            results['graph_connectivity'] = graph_conn
            
            # kBET
            kbet_score = metrics.kbet(adata, batch_key=batch_key, label_key=clustering_key)
            results['kBET'] = kbet_score
            
            print(f"Graph connectivity: {graph_conn:.4f}")
            print(f"kBET score: {kbet_score:.4f}")
            
        except Exception as e:
            print(f"Error in batch effect evaluation: {e}")
    
    return results

def batch_evaluate_scpoli_models(model_configs, output_dir='scpoli_evaluation_results'):
    """
    Batch evaluate multiple scPoli models
    
    Parameters:
    -----------
    model_configs : list of dict
        Configuration for each model
    output_dir : str
        Output directory for results
    
    Returns:
    --------
    pd.DataFrame : Comparison results
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    all_results = []
    
    for config in model_configs:
        print(f"\n{'='*60}")
        print(f"Processing: {config['name']}")
        print(f"{'='*60}")
        
        try:
            # Extract latent representation and compute UMAP
            adata = extract_latent_and_compute_umap(
                model_dir=config['model_dir'],
                adata_path=config['adata_path']
            )
            
            # Evaluate clustering
            results = evaluate_scpoli_clustering(
                adata,
                clustering_key=config['label_key'],
                reference_clustering_key=config.get('reference_label', None),
                batch_key=config.get('batch_key', 'batch'),
                embedding_key='X_umap',
                latent_key='X_scpoli'
            )
            
            # Convert to DataFrame
            results_df = pd.DataFrame([results]).T
            results_df.columns = ['Score']
            results_df.reset_index(inplace=True)
            results_df.columns = ['Metric', 'Score']
            results_df['Model'] = config['name']
            
            # Add metadata
            metadata_rows = [
                ['metadata_model_dir', config['model_dir']],
                ['metadata_label_type', config['label_key']],
                ['metadata_n_cells', adata.n_obs],
                ['metadata_n_genes', adata.n_vars],
                ['metadata_latent_dim', adata.obsm['X_scpoli'].shape[1]]
            ]
            
            for meta_key, meta_value in metadata_rows:
                new_row = pd.DataFrame({'Metric': [meta_key], 'Score': [meta_value], 'Model': [config['name']]})
                results_df = pd.concat([results_df, new_row], ignore_index=True)
            
            all_results.append(results_df)
            
            # Save individual results
            individual_output = os.path.join(output_dir, f"{config['name']}_results.csv")
            results_df.to_csv(individual_output, index=False)
            print(f"Individual results saved to: {individual_output}")
            
        except Exception as e:
            print(f"Error processing {config['name']}: {e}")
            continue
    
    # Combine and save comparison results
    if all_results:
        print(f"\n{'='*60}")
        print("Creating comparison results...")
        print(f"{'='*60}")
        
        # Combine all results
        combined_results = pd.concat(all_results, ignore_index=True)
        
        # Filter out metadata for comparison
        comparison_data = combined_results[~combined_results['Metric'].str.startswith('metadata_')]
        
        # Create comparison matrix
        comparison_df = comparison_data.pivot(index='Metric', columns='Model', values='Score')
        
        # Save comparison results
        comparison_path = os.path.join(output_dir, 'all_models_comparison.csv')
        comparison_df.to_csv(comparison_path)
        print(f"Comparison results saved to: {comparison_path}")
        
        # Save detailed results
        detailed_path = os.path.join(output_dir, 'detailed_results.csv')
        combined_results.to_csv(detailed_path, index=False)
        print(f"Detailed results saved to: {detailed_path}")
        
        return comparison_df
    
    return None

# Configuration for your models
def setup_model_configs():
    """
    Setup configuration for all your models
    """
    
    base_path = "/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/final_model"
    
    model_configs = [
        # Lineage models
        {
            'name': 'lineage_hvg2000_dim50_reseed',
            'model_dir': f'{base_path}/lineage_model_hvg2000_dim50_reseed',
            'adata_path': f'{base_path}/lineage_model_hvg2000_dim50_reseed/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_lineage',
            'model_dir': f'{base_path}/enhanced_reference_model_lineage',
            'adata_path': f'{base_path}/enhanced_reference_model_lineage/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_lineage_2ndround',
            'model_dir': f'{base_path}/enhanced_reference_model_lineage_2ndround',
            'adata_path': f'{base_path}/enhanced_reference_model_lineage_2ndround/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method1_balanced_original_lineage',
            'model_dir': f'{base_path}/method1_balanced_original_lineage',
            'adata_path': f'{base_path}/method1_balanced_original_lineage/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method3_balanced_weighted_lineage',
            'model_dir': f'{base_path}/method3_balanced_weighted_lineage',
            'adata_path': f'{base_path}/method3_balanced_weighted_lineage/adata.h5ad',
            'label_key': 'lineage',
            'batch_key': 'orig.ident'
        },
        
        # Reanno models
        {
            'name': 'reanno_hvg4000_dim50',
            'model_dir': f'{base_path}/reanno_model_hvg4000_dim50',
            'adata_path': f'{base_path}/reanno_model_hvg4000_dim50/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_reanno',
            'model_dir': f'{base_path}/enhanced_reference_model_reanno',
            'adata_path': f'{base_path}/enhanced_reference_model_reanno/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'enhanced_reference_reanno_2ndround',
            'model_dir': f'{base_path}/enhanced_reference_model_reanno_2ndround',
            'adata_path': f'{base_path}/enhanced_reference_model_reanno_2ndround/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method1_balanced_original_reanno',
            'model_dir': f'{base_path}/method1_balanced_original_reanno',
            'adata_path': f'{base_path}/method1_balanced_original_reanno/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        },
        {
            'name': 'method3_balanced_weighted_reanno',
            'model_dir': f'{base_path}/method3_balanced_weighted_reanno',
            'adata_path': f'{base_path}/method3_balanced_weighted_reanno/adata.h5ad',
            'label_key': 'reanno',
            'batch_key': 'orig.ident'
        }
    ]
    
    return model_configs

# Main execution function
def main():
    """
    Main function to run the evaluation
    """
    
    print("Setting up model configurations...")
    model_configs = setup_model_configs()
    
    print(f"Total models to evaluate: {len(model_configs)}")
    
    # Run batch evaluation
    results = batch_evaluate_scpoli_models(
        model_configs,
        output_dir='scpoli_evaluation_results'
    )
    
    if results is not None:
        print("\n" + "="*60)
        print("EVALUATION SUMMARY")
        print("="*60)
        
        # Print key metrics comparison
        key_metrics = ['silhouette_latent', 'silhouette_umap', 'n_clusters', 'graph_connectivity', 'kBET']
        
        for metric in key_metrics:
            if metric in results.index:
                print(f"\n{metric}:")
                metric_results = results.loc[metric].dropna().sort_values(ascending=False)
                for model, score in metric_results.items():
                    print(f"  {model}: {score:.4f}")
        
        print(f"\nDetailed results saved in: scpoli_evaluation_results/")
        
    else:
        print("No results generated. Please check for errors.")

# Utility function to evaluate specific model groups
def evaluate_model_group(group_type='lineage'):
    """
    Evaluate only lineage or reanno models
    
    Parameters:
    -----------
    group_type : str
        'lineage' or 'reanno'
    """
    
    all_configs = setup_model_configs()
    group_configs = [config for config in all_configs if group_type in config['label_key']]
    
    print(f"Evaluating {group_type} models only ({len(group_configs)} models)")
    
    results = batch_evaluate_scpoli_models(
        group_configs,
        output_dir=f'scpoli_{group_type}_evaluation_results'
    )
    
    return results

if __name__ == "__main__":
    # Run full evaluation
    main()
    
    # Optionally, run group-specific evaluations
    # lineage_results = evaluate_model_group('lineage')
    # reanno_results = evaluate_model_group('reanno')