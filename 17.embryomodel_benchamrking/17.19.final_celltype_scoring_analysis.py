#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplified Cell Type Correlation Analysis
Calculates correlations based solely on Composite Score within each dataset, then integrates the results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import dendrogram, linkage
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# set font
plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_data(file_path):
    """load data"""
    print("=== load data ===")
    df = pd.read_csv(file_path)
    
    datasets = df['Model'].unique()
    celltypes = df['CellType'].unique()
    
    print(f"Number of datasets: {len(datasets)}")
    print(f"Number of cell types: {len(celltypes)}")
    
    return df, datasets, celltypes

def calculate_dataset_correlation(df, dataset_name, min_score=0.01):
    """Calculate correlations between cell types within a single dataset (based on Composite Score)"""
    print(f"\nAnalyzing dataset: {dataset_name}")
    
    # Filter for current dataset
    dataset_df = df[df['Model'] == dataset_name].copy()
    
    # Create a mapping from cell type to Composite Score
    celltype_scores = dataset_df.set_index('CellType')['Composite_Score']
    
    # Remove cell types with scores below the threshold
    valid_celltypes = celltype_scores[celltype_scores >= min_score].index
    celltype_scores = celltype_scores[valid_celltypes]
    
    print(f"Number of valid cell types: {len(celltype_scores)}")
    
    if len(celltype_scores) < 2:
        print("Insufficient number of valid cell types")
        return []
    
    # Calculate correlations for all cell type pairs
    # Note: In a single dataset, each cell type has only one Composite Score value
    # So we need to calculate correlation based on multiple metrics
    
    # Create the score matrix
    score_columns = ['consistency', 'mean_certainty', 'presence', 'abundance', 
                    'Pearson_r', 'Precision', 'Recall', 'F1_Score', 'Jaccard', 'Composite_Score']
    
    # Calculate Pearson correlation
    score_matrix = dataset_df.set_index('CellType')[score_columns]
    score_matrix = score_matrix.loc[valid_celltypes]
    
    correlations = []
    
    for ct1, ct2 in combinations(valid_celltypes, 2):
        scores1 = score_matrix.loc[ct1].values
        scores2 = score_matrix.loc[ct2].values
        
        # Calculate Pearson correlation
        if np.std(scores1) > 0 and np.std(scores2) > 0:
            corr, p_value = pearsonr(scores1, scores2)
        else:
            corr, p_value = 0, 1
        
        correlations.append({
            'Dataset': dataset_name,
            'CellType1': ct1,
            'CellType2': ct2,
            'Correlation': corr,
            'P_value': p_value,
            'Composite_Score1': celltype_scores[ct1],
            'Composite_Score2': celltype_scores[ct2]
        })
    
    print(f"Calculated correlations for {len(correlations)} pairs")
    return correlations

def aggregate_correlations(all_correlations):
    """Aggregate correlation results from all datasets"""
    print("\n=== Aggregating Correlation Results ===")
    
    # Convert to DataFrame
    corr_df = pd.DataFrame(all_correlations)
    
    # Group by cell type pairs
    grouped = corr_df.groupby(['CellType1', 'CellType2'])
    
    aggregated_results = []
    
    for (ct1, ct2), group in grouped:
        # Calculate frequency and correlations
        frequency = len(group)
        correlations = group['Correlation'].values
        
        # Calculate statistical metrics
        mean_corr = np.mean(correlations)
        std_corr = np.std(correlations)
        median_corr = np.median(correlations)
        
        # Calculate consistency (consistency of correlation direction)
        positive_count = np.sum(correlations > 0)
        negative_count = np.sum(correlations < 0)
        consistency = max(positive_count, negative_count) / frequency
        
        aggregated_results.append({
            'CellType1': ct1,
            'CellType2': ct2,
            'Frequency': frequency,
            'Mean_Correlation': mean_corr,
            'Std_Correlation': std_corr,
            'Median_Correlation': median_corr,
            'Consistency': consistency,
            'Positive_Count': positive_count,
            'Negative_Count': negative_count,
            'Abs_Mean_Correlation': abs(mean_corr)
        })
    
    aggregated_df = pd.DataFrame(aggregated_results)
    
    print(f"Aggregated {len(aggregated_df)} unique cell type pairs")
    
    return aggregated_df

def create_ranking_table(aggregated_df, min_frequency=2):
    """Create a ranking table"""
    print(f"\n=== Creating Ranking Table (minimum frequency: {min_frequency}) ===")
    
    # Filter for pairs with sufficient frequency
    valid_df = aggregated_df[aggregated_df['Frequency'] >= min_frequency].copy()
    
    print(f"Number of valid pairs: {len(valid_df)}")
    
    # Sort by absolute correlation
    ranking_by_abs_corr = valid_df.sort_values('Abs_Mean_Correlation', ascending=False)
    
    # Sort by frequency
    ranking_by_freq = valid_df.sort_values('Frequency', ascending=False)
    
    # Sort by consistency
    ranking_by_consistency = valid_df.sort_values('Consistency', ascending=False)
    
    # Composite ranking (considering correlation strength, frequency, and consistency)
    valid_df['Composite_Rank_Score'] = (
        0.5 * valid_df['Abs_Mean_Correlation'] +
        0.3 * (valid_df['Frequency'] / valid_df['Frequency'].max()) +
        0.2 * valid_df['Consistency']
    )
    
    ranking_composite = valid_df.sort_values('Composite_Rank_Score', ascending=False)
    
    print("Top 10 highest correlation pairs:")
    print(ranking_by_abs_corr[['CellType1', 'CellType2', 'Frequency', 'Mean_Correlation', 'Consistency']].head(10))
    
    return {
        'by_correlation': ranking_by_abs_corr,
        'by_frequency': ranking_by_freq,
        'by_consistency': ranking_by_consistency,
        'composite': ranking_composite
    }

def create_correlation_matrix(aggregated_df, min_frequency=2):
    """Create a correlation matrix for heatmap"""
    print(f"\n=== Creating Correlation Matrix ===")
    
    # Filter for pairs with sufficient frequency
    valid_df = aggregated_df[aggregated_df['Frequency'] >= min_frequency].copy()
    
    # Get all cell types involved
    all_celltypes = set(valid_df['CellType1'].tolist() + valid_df['CellType2'].tolist())
    all_celltypes = sorted(list(all_celltypes))
    
    print(f"Number of cell types involved: {len(all_celltypes)}")
    
     # Create correlation matrix
    corr_matrix = pd.DataFrame(0.0, index=all_celltypes, columns=all_celltypes)
    
    # Fill diagonal with 1
    for ct in all_celltypes:
        corr_matrix.loc[ct, ct] = 1.0
    
    # Fill in correlation values
    for _, row in valid_df.iterrows():
        ct1, ct2 = row['CellType1'], row['CellType2']
        corr = row['Mean_Correlation']
        
        corr_matrix.loc[ct1, ct2] = corr
        corr_matrix.loc[ct2, ct1] = corr  # Symmetric matrix
    
    return corr_matrix

def create_heatmap_with_dendrogram(corr_matrix, save_path='./celltype_comparison_results/correlation_heatmap_with_dendrogram.png'):
    """Create a heatmap with dendrogram"""
    print("\n=== Creating Heatmap with Dendrogram ===")
    
    # Calculate hierarchical clustering
    distance_matrix = 1 - np.abs(corr_matrix.values)
    linkage_matrix = linkage(distance_matrix, method='ward')
    
    # Create heatmap with dendrogram
    plt.figure(figsize=(16, 14))
    
    # Use seaborn's clustermap
    g = sns.clustermap(corr_matrix, 
                       method='ward',
                       metric='euclidean',
                       cmap='RdBu_r',
                       center=0,
                       annot=False,  
                       fmt='.2f',
                       square=True,
                       linewidths=0.1,
                       figsize=(16, 14),
                       cbar_kws={"shrink": 0.8, "label": "Average correlation coefficients"})
    
    # Set title
    g.fig.suptitle('Average Correlation Coefficients among Cell Types (with Hierarchical Clustering)', fontsize=16, fontweight='bold', y=0.98)
    
    # Adjust labels
    g.ax_heatmap.set_xlabel('Cell Type', fontsize=12)
    g.ax_heatmap.set_ylabel('Cell Type', fontsize=12)
    
    # Rotate labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Heatmap saved: {save_path}")

def save_results(all_correlations, aggregated_df, rankings, corr_matrix):
    """Save all results"""
    print("\n=== Saving Results ===")
    
    # 1. Save detailed correlation results for all datasets
    corr_df = pd.DataFrame(all_correlations)
    corr_df.to_csv('./celltype_comparison_results/all_correlations_detailed.csv', index=False)
    print("Detailed correlation results saved: all_correlations_detailed.csv")
    
    # 2. Save aggregated results
    aggregated_df.to_csv('./celltype_comparison_results/aggregated_correlations.csv', index=False)
    print("Aggregated correlation results saved: aggregated_correlations.csv")
    
    # 3. Save various rankings
    for rank_type, rank_df in rankings.items():
        filename = f'./celltype_comparison_results/ranking_{rank_type}.csv'
        rank_df.to_csv(filename, index=False)
        print(f"{rank_type} ranking saved: {filename}")
    
    # 4. Save correlation matrix
    corr_matrix.to_csv('./celltype_comparison_results/mean_correlation_matrix.csv')
    print("Mean correlation matrix saved: mean_correlation_matrix.csv")
    
    # 5. Save separate results for each dataset
    for dataset in corr_df['Dataset'].unique():
        dataset_data = corr_df[corr_df['Dataset'] == dataset].copy()
        dataset_data_sorted = dataset_data.sort_values('Correlation', key=abs, ascending=False)
        filename = f'./celltype_comparison_results/correlations_{dataset.replace("/", "_").replace(" ", "_")}.csv'
        dataset_data_sorted.to_csv(filename, index=False)
        print(f"Results for dataset {dataset} saved: {filename}")

def create_summary_visualizations(aggregated_df, rankings):
    """Create summary visualizations"""
    print("\n=== Creating Summary Visualizations ===")
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Frequency distribution
    axes[0, 0].hist(aggregated_df['Frequency'], bins=range(1, aggregated_df['Frequency'].max()+2),
                   alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 0].set_xlabel('Frequency')
    axes[0, 0].set_ylabel('Number of Pairs')
    axes[0, 0].set_title('Distribution of Cell Type Pair Frequencies')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Correlation distribution
    axes[0, 1].hist(aggregated_df['Mean_Correlation'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    axes[0, 1].set_xlabel('Mean Correlation Coefficient')
    axes[0, 1].set_ylabel('Number of Pairs')
    axes[0, 1].set_title('Distribution of Mean Correlations')
    axes[0, 1].axvline(x=0, color='red', linestyle='--', alpha=0.7)
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Consistency distribution
    axes[1, 0].hist(aggregated_df['Consistency'], bins=20, alpha=0.7, color='orange', edgecolor='black')
    axes[1, 0].set_xlabel('Consistency')
    axes[1, 0].set_ylabel('Number of Pairs')
    axes[1, 0].set_title('Distribution of Correlation Consistency')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Top correlated pairs
    top_10 = rankings['composite'].head(10)
    pair_labels = [f"{row['CellType1'][:15]} - {row['CellType2'][:15]}" for _, row in top_10.iterrows()]
    
    bars = axes[1, 1].barh(range(len(top_10)), top_10['Mean_Correlation'], 
                          color='steelblue', alpha=0.7, edgecolor='black')
    axes[1, 1].set_yticks(range(len(top_10)))
    axes[1, 1].set_yticklabels(pair_labels, fontsize=8)
    axes[1, 1].set_xlabel('Mean Correlation Coefficient')
    axes[1, 1].set_title('Top 10 Highest Correlation Pairs')
    axes[1, 1].grid(True, alpha=0.3, axis='x')
    
    # Add value labels
    for i, (bar, corr) in enumerate(zip(bars, top_10['Mean_Correlation'])):
        axes[1, 1].text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                        f'{corr:.3f}', ha='left', va='center', fontweight='bold', fontsize=8)
    
    plt.tight_layout()
    plt.savefig('./celltype_comparison_results/correlation_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Summary visualization saved: correlation_summary.png")

def main():
    """Main function"""
    print("=== Simplified Cell Type Correlation Analysis ===")
    print("Analysis workflow:")
    print("1. Calculate cell type correlations within each dataset based on multiple metrics")
    print("2. Integrate results from 9 datasets and calculate frequency")
    print("3. Compute mean correlation and create heatmap")
    
    # 1. Load data
    df, datasets, celltypes = load_data('./celltype_comparison_results/all_models_celltype_metrics.csv')
    
    # 2. Calculate correlations for each dataset
    all_correlations = []
    
    for dataset in datasets:
        dataset_correlations = calculate_dataset_correlation(df, dataset)
        all_correlations.extend(dataset_correlations)
    
    if not all_correlations:
        print("No valid correlation data found")
        return
    
    # 3. Aggregate results
    aggregated_df = aggregate_correlations(all_correlations)
    
    # 4. Create rankings
    rankings = create_ranking_table(aggregated_df)
    
    # 5. Create correlation matrix
    corr_matrix = create_correlation_matrix(aggregated_df)
    
    # 6. Create heatmap
    create_heatmap_with_dendrogram(corr_matrix)
    
    # 7. Save results
    save_results(all_correlations, aggregated_df, rankings, corr_matrix)
    
    # 8. Create summary visualizations
    create_summary_visualizations(aggregated_df, rankings)
    
    print("\n=== Analysis Complete ===")
    print("Main output files:")
    print("1. all_correlations_detailed.csv - Detailed correlations for all datasets")
    print("2. aggregated_correlations.csv - Aggregated correlation results")
    print("3. ranking_*.csv - Various ranking results")
    print("4. mean_correlation_matrix.csv - Mean correlation matrix")
    print("5. correlation_heatmap_with_dendrogram.png - Heatmap with dendrogram")
    print("6. correlation_summary.png - Summary visualization")
    
    return {
        'all_correlations': all_correlations,
        'aggregated_df': aggregated_df,
        'rankings': rankings,
        'corr_matrix': corr_matrix
    }

if __name__ == "__main__":
    results = main()