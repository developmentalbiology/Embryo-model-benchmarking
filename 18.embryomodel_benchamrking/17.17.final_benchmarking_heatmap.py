#!/usr/bin/env python
# coding: utf-8
"""
Model Comparison Analysis Pipeline
This script processes pre-computed lineage metrics from multiple scRNA-seq models,
generates a long-format table for plotting, and creates heatmap visualizations.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ----------------------------
# Step 1: Process Metrics and Reshape for Heatmap
# ----------------------------
def process_metrics(input_csv):
    """
    Load the wide-format metrics CSV and reshape it into long-format 
    for flexible plotting (e.g., heatmaps).
    
    Parameters:
    -----------
    input_csv : str
        Path to the input CSV file (e.g., 'all_models_lineage_metrics.csv').
        
    Returns:
    --------
    str
        Path to the saved long-format CSV.
    """
    # Load the wide-format CSV
    df = pd.read_csv(input_csv)
    print(f"Loaded data from {input_csv}. Shape: {df.shape}")

    # Define the metrics to include in the final heatmap
    # ✅ 'Composite_Score' is included here as a key metric
    metrics_to_keep = [
        'consistency', 'mean_certainty', 'presence', 
        'abundance', 'Pearson_r', 'Precision', 'Recall', 
        'F1_Score', 'Jaccard', 'Composite_Score'  # ← Composite_Score is included
    ]
    
    # Define identifier columns
    id_columns = ['Model', 'Lineage']
    
    # Validate required columns
    missing_id_cols = [col for col in id_columns if col not in df.columns]
    missing_metric_cols = [col for col in metrics_to_keep if col not in df.columns]
    
    if missing_id_cols:
        raise ValueError(f"Missing required identifier columns: {missing_id_cols}")
    if missing_metric_cols:
        print(f"Warning: Missing metrics: {missing_metric_cols}")
        metrics_to_keep = [col for col in metrics_to_keep if col in df.columns]

    # Select and reshape
    cols_to_use = id_columns + metrics_to_keep
    df_selected = df[cols_to_use].copy()

    # Melt: Convert from wide to long format
    df_long = df_selected.melt(
        id_vars=id_columns,          
        value_vars=metrics_to_keep,  
        var_name='Metric',           
        value_name='Value'           
    )
    
    # Sort for better readability
    df_long = df_long.sort_values(['Model', 'Lineage', 'Metric']).reset_index(drop=True)

    # Save the processed data
    output_file_path = "processed_metrics_for_heatmap.csv"
    df_long.to_csv(output_file_path, index=False)
    print(f"Long-format data saved to {output_file_path}")
    print(f"Final shape: {df_long.shape}")
    print(f"Unique Models: {df_long['Model'].nunique()}")
    print(f"Unique Lineages: {df_long['Lineage'].nunique()}")
    
    return output_file_path


# ----------------------------
# Step 2: Generate Heatmaps with Consistent Color Scale
# ----------------------------
def generate_heatmaps(input_csv):
    """
    Generate one heatmap per lineage, showing all models and metrics.
    Uses a consistent color scale across all heatmaps.
    """
    # Load the long-format CSV
    df = pd.read_csv(input_csv)
    
    # Define metrics to plot
    metrics_to_plot = [
        'consistency', 'mean_certainty', 'presence', 
        'abundance', 'Pearson_r', 'Precision', 'Recall', 
        'F1_Score', 'Jaccard', 'Composite_Score'
    ]
    
    # Filter to existing metrics
    available_metrics = [m for m in metrics_to_plot if m in df.columns]
    if not available_metrics:
        print("No metrics available for heatmap.")
        return

    # Calculate global min and max for consistent color scaling
    global_min = df[available_metrics].min().min()
    global_max = df[available_metrics].max().max()

    # Create output directory
    output_dir = "./lineage_heatmaps_pdf_consistent_color_scale"
    os.makedirs(output_dir, exist_ok=True)

    # Generate a heatmap for each unique lineage
    for lineage in df['Lineage'].unique():
        # Filter data for this lineage
        lineage_df = df[df['Lineage'] == lineage]
        if lineage_df.empty:
            print(f"No data for lineage: {lineage}")
            continue

        # Pivot to create the heatmap matrix
        heatmap_data = lineage_df.set_index('Model')[available_metrics].T
        
        # Ensure all values are numeric
        heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce').fillna(0)

        if heatmap_data.empty:
            continue

        # Create the heatmap
        plt.figure(figsize=(len(heatmap_data.columns) * 0.8, len(available_metrics) * 0.5))
        sns.heatmap(
            heatmap_data,
            annot=True, 
            fmt=".2f", 
            cmap="viridis",
            cbar_kws={'label': 'Metric Value'},
            linewidths=0.5,
            vmin=global_min, 
            vmax=global_max
        )
        plt.title(f"{lineage} Lineage - Model Performance", fontsize=14)
        plt.xlabel("Models", fontsize=12)
        plt.ylabel("Metrics", fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        # Save as PDF
        output_path = os.path.join(output_dir, f"{lineage.replace(' ', '_').replace('/', '_')}_heatmap.pdf")
        plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()

    print(f"All lineage heatmaps saved to {output_dir}")


# ----------------------------
# Step 3: Generate Final Ranking Heatmap (FIXED)
# ----------------------------
def generate_ranking_heatmap(input_csv):
    """
    Generate a ranking heatmap based on Composite_Score for each lineage.
    Uses a custom 9-color palette. Handles missing (Model, Lineage) combinations.
    """
    # Load the wide-format CSV
    df = pd.read_csv(input_csv)
    
    if 'Composite_Score' not in df.columns:
        print("Error: 'Composite_Score' column not found.")
        return

    # Get unique lineages and models
    lineages = df['Lineage'].unique()
    models = df['Model'].unique()

    # Create a complete pivot table with all (Model, Lineage) combinations
    # Fill missing values with 0.0 (worst possible score)
    pivot_df = df.pivot_table(
        index='Lineage',
        columns='Model',
        values='Composite_Score',
        fill_value=0.0  # Critical: Fill missing combinations with 0
    )

    # Now pivot_df has a row for every lineage and a column for every model
    # Missing combinations (e.g., Model_A has no Notochord) are filled with 0.0

    # Create a DataFrame to store rankings
    ranking_df = pd.DataFrame(index=lineages, columns=models)

    # Calculate dense ranking for each lineage
    for lineage in lineages:
        # Get the Composite_Score for all models in this lineage
        scores = pivot_df.loc[lineage, :]
        
        # Rank models by Composite_Score, descending (higher score = better rank)
        # 'min' or 'dense' both work; 'min' gives ties the same rank
        ranking_df.loc[lineage] = scores.rank(method='dense', ascending=False).astype(int)

    # Convert to numeric and handle any remaining NaN (shouldn't be any)
    ranking_df = ranking_df.apply(pd.to_numeric, errors='coerce').fillna(len(models)).astype(int)

    # Define a custom 9-color palette
    custom_palette = [
        "#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#e6f598", 
        "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"
    ]
    cmap = sns.color_palette(custom_palette, as_cmap=True)

    # Create the ranking heatmap
    plt.figure(figsize=(len(models) * 0.8, len(lineages) * 0.5))
    sns.heatmap(
        ranking_df,
        annot=True, 
        fmt=".0f", 
        cmap=cmap,
        cbar_kws={'label': 'Rank'}, 
        linewidths=0.5,
        vmin=1, 
        vmax=len(models)  # Dynamic max rank based on number of models
    )
    plt.title("Model Rankings by Lineage (Dense Ranking)", fontsize=14)
    plt.xlabel("Models", fontsize=12)
    plt.ylabel("Lineages", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save as PDF
    output_dir = "./final_ranking_heatmaps_dense_ranking"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "final_model_ranking.pdf")
    plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Final ranking heatmap saved to {output_path}")


# ----------------------------
# Main Execution
# ----------------------------
if __name__ == "__main__":
    # Set the path to your combined metrics CSV
    # This file should be generated by your main analysis script
    input_csv = "/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250806_embryo_model_benchmarking/lineage_comparison_results/all_models_lineage_metrics.csv"

    # Step 1: Process the data into long format for plotting
    processed_csv = process_metrics(input_csv)

    # Step 2: Generate individual lineage heatmaps
    generate_heatmaps(input_csv)  # Use the original wide-format CSV

    # Step 3: Generate the final ranking heatmap
    generate_ranking_heatmap(input_csv)  # Use the original wide-format CSV

    print("Analysis completed successfully.")

