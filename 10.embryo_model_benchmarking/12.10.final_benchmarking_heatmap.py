#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Step 1: Process Metrics and Save to CSV
def process_metrics(input_csv):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(input_csv)

    # Define metrics to include in the final table
    metrics_to_plot = [
        'consistency', 'presence', 'percentage_certain_cells', 'mean_certainty',
        'abundance', 'Mean_Pearson_r',
        'Mean_Precision', 'Mean_Recall', 'Mean_F1_Score', 'Mean_Jaccard'
    ]

    processed_data = []

    # Group by Lineage and Model
    grouped = df.groupby(['Lineage', 'Model'])
    for (lineage, model), group in grouped:
        direct_row = group[group['Type'] == 'Direct'].iloc[0] if not group[group['Type'] == 'Direct'].empty else None
        derived_row = group[group['Type'] == 'Derived'].iloc[0] if not group[group['Type'] == 'Derived'].empty else None

        metrics_dict = {'Lineage': lineage, 'Model': model}
        for metric in metrics_to_plot:
            if metric in ['Mean_Pearson_r', 'Mean_Precision', 'Mean_Recall', 'Mean_F1_Score', 'Mean_Jaccard']:
                metrics_dict[metric] = derived_row[metric] if derived_row is not None and metric in derived_row else 0
            elif metric in ['abundance', 'percentage_certain_cells', 'mean_certainty', 'presence']:
                direct_value = direct_row[metric] if direct_row is not None and metric in direct_row else 0
                derived_value = derived_row[metric] if derived_row is not None and metric in derived_row else 0
                metrics_dict[metric] = (direct_value + derived_value) / 2
            elif metric == 'consistency':
                metrics_dict[metric] = direct_row[metric] if direct_row is not None and metric in direct_row else 0

        processed_data.append(metrics_dict)

    # Save processed data to a new CSV file
    processed_df = pd.DataFrame(processed_data)
    output_file_path = "processed_metrics_for_plotting.csv"
    processed_df.to_csv(output_file_path, index=False)
    print(f"Processed data saved to {output_file_path}")
    return output_file_path


# Step 2: Calculate Final Score
def calculate_final_scores(input_csv):
    METRIC_WEIGHTS = {
        'Mean_Precision': 0.10,
        'Mean_Recall': 0.15,
        'Mean_F1_Score': 0.15,
        'Mean_Jaccard': 0.15,
        'Mean_Pearson_r': 0.15,
        'consistency': 0.05,
        'abundance': 0.1,
        'presence': 0.05,
        'percentage_certain_cells': 0.05,
        'mean_certainty': 0.05
    }

    def calculate_final_score(row):
        score = 0
        for metric, weight in METRIC_WEIGHTS.items():
            if metric in row and not pd.isna(row[metric]):
                score += row[metric] * weight
        return score

    # Load processed CSV and calculate final scores
    df = pd.read_csv(input_csv)
    df['Final_Score'] = df.apply(calculate_final_score, axis=1)

    # Save updated DataFrame with final scores
    output_file_path = "processed_metrics_with_final_score.csv"
    df.to_csv(output_file_path, index=False)
    print(f"Processed data with final scores saved to {output_file_path}")
    return output_file_path


# Step 3: Generate Heatmaps with Consistent Color Scale
def generate_heatmaps(input_csv):
    # Load processed CSV
    df = pd.read_csv(input_csv)

    # Define metrics to include in the heatmap
    metrics_to_plot = [
        'consistency', 'presence', 'percentage_certain_cells', 'mean_certainty',
        'abundance', 'Mean_Pearson_r',
        'Mean_Precision', 'Mean_Recall', 'Mean_F1_Score', 'Mean_Jaccard',
        'Final_Score'
    ]

    # Get global min and max for consistent color scaling
    global_min = df[metrics_to_plot].min().min()
    global_max = df[metrics_to_plot].max().max()

    # Create output directory for heatmaps
    output_dir = "./lineage_heatmaps_pdf_consistent_color_scale"
    os.makedirs(output_dir, exist_ok=True)

    # Generate a heatmap for each lineage
    for lineage in df['Lineage'].unique():
        lineage_df = df[df['Lineage'] == lineage]
        heatmap_data = lineage_df.set_index('Model')[metrics_to_plot].T
        heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce').fillna(0)

        if heatmap_data.empty:
            print(f"No data available for lineage: {lineage}")
            continue

        plt.figure(figsize=(len(heatmap_data.columns) * 0.8, len(metrics_to_plot) * 0.5))
        sns.heatmap(
            heatmap_data,
            annot=True, fmt=".2f", cmap="viridis",
            cbar_kws={'label': 'Metric Value'},
            linewidths=0.5,
            vmin=global_min, vmax=global_max
        )
        plt.title(f"{lineage} Lineage - Metrics Across Models", fontsize=14)
        plt.xlabel("Models", fontsize=12)
        plt.ylabel("Metrics", fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        # Save heatmap as PDF
        output_path = os.path.join(output_dir, f"{lineage}_heatmap.pdf")
        plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()

    print("All heatmaps generated and saved as PDF files with consistent color scale.")


# Step 4: Generate Ranking Heatmap
def generate_ranking_heatmap(input_csv):
    # Load processed CSV
    df = pd.read_csv(input_csv)

    # Define a nine-color palette for ranking
    custom_palette = [
         "#9e0142","#d53e4f","#f46d43","#fdae61","#e6f598", 
         "#abdda4","#66c2a5","#3288bd", "#5e4fa2",  
    ]
    cmap = sns.color_palette(custom_palette, as_cmap=True)

    # Get unique lineages and models
    lineages = df['Lineage'].unique()
    models = df['Model'].unique()

    # Create an empty DataFrame to store rankings
    ranking_df = pd.DataFrame(index=lineages, columns=models)

    # Calculate dense rankings for each lineage
    for lineage in lineages:
        lineage_df = df[df['Lineage'] == lineage]
        ranked_models = lineage_df.sort_values(by='Final_Score', ascending=False)
        ranked_models['Rank'] = ranked_models['Final_Score'].rank(method='dense', ascending=False).astype(int)

        for _, row in ranked_models.iterrows():
            ranking_df.loc[lineage, row['Model']] = row['Rank']

    # Convert ranking DataFrame to numeric
    ranking_df = ranking_df.apply(pd.to_numeric, errors='coerce')

    # Generate the ranking heatmap
    plt.figure(figsize=(len(models) * 0.8, len(lineages) * 0.5))
    sns.heatmap(
        ranking_df,
        annot=True, fmt=".0f", cmap=cmap,
        cbar_kws={'label': 'Rank'}, linewidths=0.5,
        vmin=1, vmax=9
    )
    plt.title("Final Rankings of Models Across Lineages", fontsize=14)
    plt.xlabel("Models", fontsize=12)
    plt.ylabel("Lineages", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save ranking heatmap as PDF
    output_dir = "./final_ranking_heatmaps_dense_ranking"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "final_ranking_heatmap_dense_ranking.pdf")
    plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Final ranking heatmap saved to {output_path}")


# Main Execution
if __name__ == "__main__":
    # Input CSV file path
    input_csv = "/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/lineage_comparison_results/all_models_lineage_metrics.csv"

    # Step 1: Process metrics
    processed_csv = process_metrics(input_csv)

    # Step 2: Calculate final scores
    scored_csv = calculate_final_scores(processed_csv)

    # Step 3: Generate heatmaps with consistent color scale
    generate_heatmaps(scored_csv)

    # Step 4: Generate ranking heatmap
    generate_ranking_heatmap(scored_csv)