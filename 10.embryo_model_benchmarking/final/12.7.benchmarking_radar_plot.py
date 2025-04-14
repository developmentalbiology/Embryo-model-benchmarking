#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Model name mapping to handle different naming conventions
MODEL_NAME_MAPPING = {
    'Ai_model': 'corrected_processed_Ai_model',
    'Hislop': 'corrected_processed_Hislop',
    'Liu': 'corrected_processed_Liu',
    'Oldak': 'corrected_processed_Oldak',
    'Pedroza': 'corrected_processed_Pedroza',
    'Rowan': 'corrected_processed_Rowan',
    'Weatherbee': 'corrected_processed_Weatherbee',
    'zheng_2019': 'corrected_processed_zheng_2019',
    'zheng_2022': 'corrected_processed_zheng_2022'
}

# Reverse mapping (for display in plots)
DISPLAY_NAME_MAPPING = {v: k for k, v in MODEL_NAME_MAPPING.items()}


def load_metrics_data(metrics_dir):
    """
    Load all benchmark metrics from different CSV files
    
    Parameters
    ----------
    metrics_dir : str
        Directory containing metrics files
    
    Returns
    -------
    pandas.DataFrame
        Combined metrics dataframe
    """
    # Define the metrics files to load
    files_to_load = {
        'celltype_certainty': 'celltype_certainty_comparison.csv',
        'lineage_certainty': 'lineage_certainty_comparison.csv',
        'celltype_mean_certainty': 'celltype_overall_mean_certainty_certain_cells_comparison.csv',
        'lineage_mean_certainty': 'lineage_overall_mean_certainty_certain_cells_comparison.csv',
        'celltype_correlation': 'celltype_correlation_comparison.csv',
        'lineage_correlation': 'lineage_correlation_comparison.csv',
        'celltype_presence': 'celltype_presence_comparison.csv',
        'lineage_presence': 'lineage_presence_comparison.csv',
        'celltype_marker_overlap': 'combined_celltype_metrics.csv',
        'lineage_marker_overlap': 'combined_lineage_metrics.csv'
    }
    
    all_metrics = {}
    loaded_models = set()
    
    # Load each metrics file
    for key, filename in files_to_load.items():
        filepath = os.path.join(metrics_dir, filename)
        if os.path.exists(filepath):
            print(f"Loading {key} metrics from {filepath}")
            df = pd.read_csv(filepath)
            
            # Standardize Model column name if it's different
            if 'model' in df.columns and 'Model' not in df.columns:
                df.rename(columns={'model': 'Model'}, inplace=True)
                
            # Standardize model names to ensure consistency between files
            if key in ['celltype_marker_overlap', 'lineage_marker_overlap']:
                # For the marker overlap files, map from short names to long names
                df['Model'] = df['Model'].map(MODEL_NAME_MAPPING).fillna(df['Model'])
            else:
                # For other files, ensure consistency by keeping long names
                pass
            
            # Extract relevant columns
            if key == 'celltype_certainty':
                metrics_df = df[['Model', 'Certain Percentage']].copy()
                metrics_df.rename(columns={'Certain Percentage': 'celltype_certainty_percentage'}, inplace=True)
            
            elif key == 'lineage_certainty':
                metrics_df = df[['Model', 'Certain Percentage']].copy()
                metrics_df.rename(columns={'Certain Percentage': 'lineage_certainty_percentage'}, inplace=True)
            
            elif key == 'celltype_mean_certainty':
                if 'Mean_Certainty' in df.columns:
                    metrics_df = df[['Model', 'Mean_Certainty']].copy()
                    metrics_df.rename(columns={'Mean_Certainty': 'celltype_mean_certainty'}, inplace=True)
                    # Convert to percentage
                    metrics_df['celltype_mean_certainty'] *= 100
                else:
                    continue
            
            elif key == 'lineage_mean_certainty':
                if 'Mean_Certainty' in df.columns:
                    metrics_df = df[['Model', 'Mean_Certainty']].copy()
                    metrics_df.rename(columns={'Mean_Certainty': 'lineage_mean_certainty'}, inplace=True)
                    # Convert to percentage
                    metrics_df['lineage_mean_certainty'] *= 100
                else:
                    continue
            
            elif key == 'celltype_correlation':
                if 'Pearson r' in df.columns:
                    metrics_df = df[['Model', 'Pearson r']].copy()
                    metrics_df.rename(columns={'Pearson r': 'celltype_correlation'}, inplace=True)
                    # Scale to 0-100 range
                    metrics_df['celltype_correlation'] = (metrics_df['celltype_correlation'] + 1) * 50
                else:
                    continue
            
            elif key == 'lineage_correlation':
                if 'Pearson r' in df.columns:
                    metrics_df = df[['Model', 'Pearson r']].copy()
                    metrics_df.rename(columns={'Pearson r': 'lineage_correlation'}, inplace=True)
                    # Scale to 0-100 range
                    metrics_df['lineage_correlation'] = (metrics_df['lineage_correlation'] + 1) * 50
                else:
                    continue
            
            elif key == 'celltype_presence':
                # Check for different column names that might contain the presence score
                presence_cols = [
                    'Cell Type Presence (%)', 
                    'Weighted Lineage Presence (%)',
                    'Presence Score (%)',
                    'Major Types Presence (%)'
                ]
                
                # Use first matching column found
                for col in presence_cols:
                    if col in df.columns:
                        print(f"  Using '{col}' column for cell type presence")
                        metrics_df = df[['Model', col]].copy()
                        metrics_df.rename(columns={col: 'celltype_presence_percentage'}, inplace=True)
                        break
                else:
                    # If no standard column found, check all columns
                    print("  Searching for presence columns in celltype_presence file...")
                    print(f"  Available columns: {df.columns.tolist()}")
                    continue
            
            elif key == 'lineage_presence':
                # Check for different column names that might contain the presence score
                presence_cols = [
                    'Weighted Lineage Presence (%)',
                    'Presence Score (%)',
                    'Major Types Presence (%)'
                ]
                
                # Use first matching column found
                for col in presence_cols:
                    if col in df.columns:
                        print(f"  Using '{col}' column for lineage presence")
                        metrics_df = df[['Model', col]].copy()
                        metrics_df.rename(columns={col: 'lineage_presence_percentage'}, inplace=True)
                        break
                else:
                    # If no standard column found, check all columns
                    print("  Searching for presence columns in lineage_presence file...")
                    print(f"  Available columns: {df.columns.tolist()}")
                    continue
                
            elif key == 'celltype_marker_overlap':
                # Check if this is the marker overlap file with different model names
                needed_cols = ['Model', 'mean_precision', 'mean_recall', 'mean_f1', 'mean_jaccard']
                available_cols = [col for col in needed_cols if col in df.columns]
                
                if len(available_cols) < 2:  # Need at least Model and one metric
                    continue
                    
                metrics_df = df[available_cols].copy()
                
                # Rename columns to standard format
                column_mapping = {
                    'mean_precision': 'celltype_marker_precision',
                    'mean_recall': 'celltype_marker_recall',
                    'mean_f1': 'celltype_marker_f1',
                    'mean_jaccard': 'celltype_marker_jaccard'
                }
                
                for old_col, new_col in column_mapping.items():
                    if old_col in metrics_df.columns:
                        metrics_df.rename(columns={old_col: new_col}, inplace=True)
                        # Convert to percentage scale
                        metrics_df[new_col] *= 100
            
            elif key == 'lineage_marker_overlap':
                # Check if this is the marker overlap file with different model names
                needed_cols = ['Model', 'mean_precision', 'mean_recall', 'mean_f1', 'mean_jaccard']
                available_cols = [col for col in needed_cols if col in df.columns]
                
                if len(available_cols) < 2:  # Need at least Model and one metric
                    continue
                    
                metrics_df = df[available_cols].copy()
                
                # Rename columns to standard format
                column_mapping = {
                    'mean_precision': 'lineage_marker_precision',
                    'mean_recall': 'lineage_marker_recall',
                    'mean_f1': 'lineage_marker_f1',
                    'mean_jaccard': 'lineage_marker_jaccard'
                }
                
                for old_col, new_col in column_mapping.items():
                    if old_col in metrics_df.columns:
                        metrics_df.rename(columns={old_col: new_col}, inplace=True)
                        # Convert to percentage scale
                        metrics_df[new_col] *= 100
                
            all_metrics[key] = metrics_df
            loaded_models.update(metrics_df['Model'].tolist())
        else:
            print(f"Warning: {filepath} not found")
    
    # List of all models from all files
    all_model_names = sorted(loaded_models)
    print(f"Found {len(all_model_names)} unique models across all metrics files")
    
    # Merge all metrics on 'Model' column
    combined_df = pd.DataFrame({'Model': all_model_names})
    
    for key, df in all_metrics.items():
        combined_df = pd.merge(combined_df, df, on='Model', how='left')
    
    # Fill NaNs with 0 for safer calculations
    combined_df.fillna(0, inplace=True)
    
    # Add display names for better visualization
    combined_df['Display_Name'] = combined_df['Model'].map(DISPLAY_NAME_MAPPING).fillna(combined_df['Model'])
    
    # Manual assignment for presence scores if not loaded correctly
    # If you have a specific table with presence scores, you can load it directly here
    presence_file = os.path.join(metrics_dir, "lineage_presence_comparison.csv")
    if os.path.exists(presence_file):
        try:
            presence_df = pd.read_csv(presence_file)
            print("Loading presence scores directly from lineage_presence_comparison.csv")
            
            # If your table structure is as shown in the example
            if 'Model' in presence_df.columns and 'Presence Score (%)' in presence_df.columns:
                print("Using 'Presence Score (%)' for lineage presence")
                for idx, row in presence_df.iterrows():
                    model = row['Model']
                    if model in combined_df['Model'].values:
                        combined_df.loc[combined_df['Model'] == model, 'lineage_presence_percentage'] = row['Presence Score (%)']
            
            # The table might have different column names
            elif 'Model' in presence_df.columns:
                print(f"Available columns in presence table: {presence_df.columns.tolist()}")
                
                # Try to find the presence column
                for col in presence_df.columns:
                    if "presence" in col.lower() or "score" in col.lower():
                        print(f"Using '{col}' column as presence score")
                        for idx, row in presence_df.iterrows():
                            model = row['Model']
                            if model in combined_df['Model'].values:
                                combined_df.loc[combined_df['Model'] == model, 'lineage_presence_percentage'] = row[col]
                                # Also use this for cell type presence if not already set
                                if combined_df.loc[combined_df['Model'] == model, 'celltype_presence_percentage'].iloc[0] == 0:
                                    combined_df.loc[combined_df['Model'] == model, 'celltype_presence_percentage'] = row[col]
        except Exception as e:
            print(f"Error loading presence file: {str(e)}")
    
    # Print loaded metrics for verification
    print("\nLoaded metrics:")
    for col in combined_df.columns:
        if col not in ['Model', 'Display_Name']:
            print(f"{col}: {combined_df[col].mean():.2f} (mean), {combined_df[col].min():.2f} (min), {combined_df[col].max():.2f} (max)")
    
    return combined_df

def calculate_composite_scores(df):
    """
    Calculate composite scores from raw metrics
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with raw metrics
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with raw metrics and calculated scores
    """
    result_df = df.copy()
    
    # Create metric groupings for aggregation
    metric_groups = {
        'certainty': {
            'celltype': ['celltype_certainty_percentage', 'celltype_mean_certainty'],
            'lineage': ['lineage_certainty_percentage', 'lineage_mean_certainty']
        },
        'marker_overlap': {
            'celltype': ['celltype_marker_precision', 'celltype_marker_recall', 'celltype_marker_f1', 'celltype_marker_jaccard'],
            'lineage': ['lineage_marker_precision', 'lineage_marker_recall', 'lineage_marker_f1', 'lineage_marker_jaccard']
        },
        'correlation': {
            'celltype': ['celltype_correlation'],
            'lineage': ['lineage_correlation']
        },
        'presence': {
            'celltype': ['celltype_presence_percentage'],
            'lineage': ['lineage_presence_percentage']
        }
    }
    
    # Weights for each metric within groups
    metric_weights = {
        'celltype_certainty_percentage': 0.5,
        'celltype_mean_certainty': 0.5,
        'lineage_certainty_percentage': 0.5,
        'lineage_mean_certainty': 0.5,
        
        'celltype_marker_precision': 0.25,
        'celltype_marker_recall': 0.25,
        'celltype_marker_f1': 0.25,
        'celltype_marker_jaccard': 0.25,
        
        'lineage_marker_precision': 0.25,
        'lineage_marker_recall': 0.25,
        'lineage_marker_f1': 0.25,
        'lineage_marker_jaccard': 0.25,
        
        'celltype_correlation': 1.0,
        'lineage_correlation': 1.0,
        
        'celltype_presence_percentage': 1.0,
        'lineage_presence_percentage': 1.0
    }
    
    # Normalize metrics before calculating scores
    normalized_df = result_df.copy()
    
    # Create a dataframe to store normalized values
    for col in result_df.columns:
        if col != 'Model' and col != 'Display_Name' and col in metric_weights:
            # Get min and max for normalization
            col_min = result_df[col].min()
            col_max = result_df[col].max()
            
            # Normalize to [0, 1] range
            if col_max > col_min:
                normalized_df[col] = (result_df[col] - col_min) / (col_max - col_min)
            else:
                normalized_df[col] = 0
    
    # Calculate composite scores for each group
    for group, type_metrics in metric_groups.items():
        # Calculate celltype and lineage scores for this group
        for type_name, metrics in type_metrics.items():
            score_name = f"{type_name}_{group}_score"
            
            # Check which metrics are available
            available_metrics = [m for m in metrics if m in normalized_df.columns]
            
            if not available_metrics:
                result_df[score_name] = 0
                continue
                
            # Calculate weighted sum
            weighted_sum = 0
            total_weight = 0
            
            for metric in available_metrics:
                weight = metric_weights.get(metric, 1.0)
                weighted_sum += normalized_df[metric] * weight
                total_weight += weight
            
            if total_weight > 0:
                result_df[score_name] = 100 * weighted_sum / total_weight
            else:
                result_df[score_name] = 0
                
    # Create overall aggregated scores
    # First aggregate by type
    for type_name in ['celltype', 'lineage']:
        scores = [f"{type_name}_{group}_score" for group in metric_groups.keys()]
        available_scores = [s for s in scores if s in result_df.columns]
        
        if available_scores:
            result_df[f"{type_name}_composite_score"] = result_df[available_scores].mean(axis=1)
    
    # Then create an overall score
    if 'celltype_composite_score' in result_df.columns and 'lineage_composite_score' in result_df.columns:
        result_df['overall_composite_score'] = (result_df['celltype_composite_score'] + result_df['lineage_composite_score']) / 2
    elif 'celltype_composite_score' in result_df.columns:
        result_df['overall_composite_score'] = result_df['celltype_composite_score']
    elif 'lineage_composite_score' in result_df.columns:
        result_df['overall_composite_score'] = result_df['lineage_composite_score'] 
    
    # Add ranks
    if 'celltype_composite_score' in result_df.columns:
        result_df['celltype_rank'] = result_df['celltype_composite_score'].rank(ascending=False)
    
    if 'lineage_composite_score' in result_df.columns:
        result_df['lineage_rank'] = result_df['lineage_composite_score'].rank(ascending=False)
    
    if 'overall_composite_score' in result_df.columns:
        result_df['overall_rank'] = result_df['overall_composite_score'].rank(ascending=False)
    
    return result_df


def create_radar_chart(df, metrics, metric_labels, title, score_column, output_file, model_filter=None, color_mapping=None):
    """
    Create a radar chart for specified metrics with filled areas,
    with values scaled using Z-score standardization
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with metrics
    metrics : list
        List of metrics to plot
    metric_labels : dict
        Dictionary mapping metric names to display labels
    title : str
        Title of the plot
    score_column : str
        Column name for score used in legend
    output_file : str
        Path to save the output file
    model_filter : list, optional
        List of model names to include (use None for all models)
    color_mapping : dict, optional
        Mapping of model names to colors for consistency
    """
    if len(metrics) < 3:
        print(f"Error: Need at least 3 metrics for radar chart. Only found {len(metrics)}")
        return False
    
    # Filter models if specified
    if model_filter:
        plot_df = df[df['Model'].isin(model_filter)].copy()
    else:
        # Use all models but warn if there are too many
        plot_df = df.copy()
        if len(plot_df) > 10:
            print(f"Warning: Plotting {len(plot_df)} models on one radar chart may be crowded")
    
    if len(plot_df) == 0:
        print("Error: No models to plot after filtering")
        return False
    
    # Number of variables
    N = len(metrics)
    
    # Calculate angles for each category
    angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    
    # Close the loop
    angles += angles[:1]
    metrics_loop = metrics + [metrics[0]]
    
    # Create figure - use a larger figure if plotting all models
    if model_filter is None and len(plot_df) > 5:
        fig = plt.figure(figsize=(14, 14))
    else:
        fig = plt.figure(figsize=(10, 10))
    
    ax = plt.subplot(111, polar=True)
    
    # Add labels to chart
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([metric_labels.get(m, m) for m in metrics])
    
    # Create scaled dataframe for radar chart
    # Scale each metric so that max value = 100
    scaled_df = plot_df.copy()
    for metric in metrics:
        if metric in df.columns:
            max_value = df[metric].max()
            if max_value > 0:  # Avoid division by zero
                scaled_df[f"{metric}_scaled"] = (plot_df[metric] / max_value) * 100
            else:
                scaled_df[f"{metric}_scaled"] = 0
    
    # Set y ticks with Z-score equivalents
    y_ticks = [0, 25, 50, 75, 100]
    z_equiv = [-2, -1, 0, 1, 2]  # Corresponding Z-scores
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f"{y} (Z={z})" for y, z in zip(y_ticks, z_equiv)])
    
    # Set y limits
    ax.set_ylim(0, 105)  # Add a little extra for labels
    
    # Use provided color mapping or create a new one if not provided
    if color_mapping is None:
        # Define a standard set of colors for consistency
        if len(df['Model'].unique()) > 10:
            # Use tab20 for more models
            colors = plt.cm.tab20(np.linspace(0, 1, len(df['Model'].unique())))
        else:
            # Use tab10 for fewer models
            colors = plt.cm.tab10(np.linspace(0, 1, len(df['Model'].unique())))
        
        # Create color mapping for all models
        color_mapping = {model: colors[i] for i, model in enumerate(df['Model'].unique())}
    
    # Sort by score for better visualization (best models on top)
    scaled_df = scaled_df.sort_values(score_column, ascending=True)
    
    # Plot each model with filled area
    for _, row in scaled_df.iterrows():
        model = row['Model']
        
        # Get color for this model
        color = color_mapping.get(model, 'gray')  # Use gray as fallback
        
        # Get scaled values for this model
        values = []
        for m in metrics:
            scaled_metric = f"{m}_scaled"
            if scaled_metric in row:
                values.append(row[scaled_metric])
            else:
                values.append(50)  # Default to middle (Z-score of 0)
        
        # Close the loop - add first value at the end
        values.append(values[0])
        
        # Get display name for model
        display_name = row['Display_Name'] if 'Display_Name' in row else row['Model']
        
        # Plot values with filled area
        ax.plot(angles, values, 'o-', linewidth=2, markersize=4,
                label=f"{display_name} ({row[score_column]:.1f})", color=color)
        
        # Fill with semi-transparent color - reduce alpha for many models
        alpha = 0.3 if len(plot_df) <= 5 else 0.2
        ax.fill(angles, values, color=color, alpha=alpha)
    
    # Add legend - adjust position based on number of models
    if len(plot_df) > 5:
        # For many models, place legend below the chart
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                  ncol=min(3, len(plot_df)), fontsize=8)
    else:
        # For few models, place legend in bottom right
        plt.legend(loc='lower right', bbox_to_anchor=(0.15, 0.1))
    
    # Add title
    plt.title(title + "\n(Z-score Standardized)", size=15, y=1.1)
    
    # Add grid
    ax.grid(True)
    
    # Add explanation of Z-score scaling
    plt.figtext(0.5, -0.05, 
                "Z-score scaling: 50 = mean, 75 = +1 std dev, 25 = -1 std dev",
                ha='center', fontsize=10)
    
    # Save as PNG
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    # Also save as PDF
    pdf_output = output_file.replace('.png', '.pdf')
    plt.savefig(pdf_output, dpi=300, bbox_inches='tight', format='pdf')
    print(f"Saved radar chart as PNG and PDF: {output_file} and {pdf_output}")
    
    plt.close()
    
    # Return the color mapping for reuse
    return color_mapping

def create_plots(df, output_dir):
    """Create all visualization plots with specified metrics"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define readable metric labels
    metric_labels = {
        # Cell Type metrics
        'celltype_presence_percentage': 'Cell Type\nPresence %',
        'celltype_certainty_percentage': 'Cell Type\nCertainty %',
        'celltype_mean_certainty': 'Cell Type\nMean Certainty',
        'celltype_correlation': 'Cell Type\nExpression Correlation',
        'celltype_marker_recall': 'Cell Type\nMarker Recall',
        'celltype_marker_precision': 'Cell Type\nMarker Precision',
        'celltype_marker_f1': 'Cell Type\nMarker F1',
        'celltype_marker_jaccard': 'Cell Type\nMarker Jaccard',
        
        
        # Lineage metrics
        'lineage_presence_percentage': 'Lineage\nPresence %',
        'lineage_certainty_percentage': 'Lineage\nCertainty %',
        'lineage_mean_certainty': 'Lineage\nMean Certainty',
        'lineage_correlation': 'Lineage\nExpression Correlation',
        'lineage_marker_recall': 'Lineage\nMarker Recall',
        'lineage_marker_precision': 'Lineage\nMarker Precision',
        'lineage_marker_f1': 'Lineage\nMarker F1',
        'lineage_marker_jaccard': 'Lineage\nMarker Jaccard',
        
        
        # Aggregated scores
        'celltype_certainty_score': 'Cell Type\nCertainty',
        'lineage_certainty_score': 'Lineage\nCertainty',
        'celltype_marker_overlap_score': 'Cell Type\nMarker Overlap',
        'lineage_marker_overlap_score': 'Lineage\nMarker Overlap',
        'celltype_correlation_score': 'Cell Type\nCorrelation',
        'lineage_correlation_score': 'Lineage\nCorrelation',
        'celltype_presence_score': 'Cell Type\nPresence',
        'lineage_presence_score': 'Lineage\nPresence'
    }

    # Create consistent color mapping for all models
    all_models = df['Model'].unique()
    if len(all_models) > 10:
        colors = plt.cm.tab20(np.linspace(0, 1, len(all_models)))
    else:
        colors = plt.cm.tab10(np.linspace(0, 1, len(all_models)))
    
    color_mapping = {model: colors[i] for i, model in enumerate(all_models)}
    
    # Reordered to match your image (mean certainty first, correlation next to marker recall)
    celltype_metrics = [
        'celltype_mean_certainty',
        'celltype_certainty_percentage', 
        'celltype_presence_percentage',
        'celltype_marker_precision',
        'celltype_marker_jaccard',
        'celltype_marker_f1',
        'celltype_marker_recall',
        'celltype_correlation'
    ]
    
    lineage_metrics = [
        'lineage_mean_certainty',
        'lineage_certainty_percentage',
        'lineage_presence_percentage',
        'lineage_marker_precision',
        'lineage_marker_jaccard',
        'lineage_marker_f1',
        'lineage_marker_recall',
        'lineage_correlation'
    ]
    
    # Filter to only include metrics that exist in the dataframe
    celltype_metrics = [m for m in celltype_metrics if m in df.columns]
    lineage_metrics = [m for m in lineage_metrics if m in df.columns]
    
    # Show which metrics are being used
    print(f"Cell type metrics being used: {celltype_metrics}")
    print(f"Lineage metrics being used: {lineage_metrics}")
    
    # Create radar charts for ALL models
    
    # Cell type metrics radar chart for ALL models
    if len(celltype_metrics) >= 3:
        create_radar_chart(
            df,
            celltype_metrics,
            metric_labels,
            'Cell Type Metrics Comparison (All Models)',
            'celltype_composite_score' if 'celltype_composite_score' in df.columns else 'overall_composite_score',
            os.path.join(output_dir, 'celltype_metrics_radar_all.png'),
            model_filter=None,  # No filter = all models
            color_mapping=color_mapping  # Use consistent colors
        )
    else:
        print(f"Warning: Not enough cell type metrics to create radar chart (need at least 3, found {len(celltype_metrics)})")
    
    # Lineage metrics radar chart for ALL models
    if len(lineage_metrics) >= 3:
        create_radar_chart(
            df,
            lineage_metrics,
            metric_labels,
            'Lineage Metrics Comparison (All Models)',
            'lineage_composite_score' if 'lineage_composite_score' in df.columns else 'overall_composite_score',
            os.path.join(output_dir, 'lineage_metrics_radar_all.png'),
            model_filter=None,  # No filter = all models
            color_mapping=color_mapping  # Use consistent colors
        )
    else:
        print(f"Warning: Not enough lineage metrics to create radar chart (need at least 3, found {len(lineage_metrics)})")
    
    # Also create top 5 model charts for better readability
    top_models = df.sort_values('overall_composite_score', ascending=False).head(5)['Model'].tolist()
    
    # Cell type metrics radar chart for TOP 5
    if len(celltype_metrics) >= 3:
        create_radar_chart(
            df,
            celltype_metrics,
            metric_labels,
            'Cell Type Metrics Comparison (Top 5)',
            'celltype_composite_score' if 'celltype_composite_score' in df.columns else 'overall_composite_score',
            os.path.join(output_dir, 'celltype_metrics_radar_top5.png'),
            model_filter=top_models,
            color_mapping=color_mapping  # Use consistent colors
        )
    
    # Lineage metrics radar chart for TOP 5
    if len(lineage_metrics) >= 3:
        create_radar_chart(
            df,
            lineage_metrics,
            metric_labels,
            'Lineage Metrics Comparison (Top 5)',
            'lineage_composite_score' if 'lineage_composite_score' in df.columns else 'overall_composite_score',
            os.path.join(output_dir, 'lineage_metrics_radar_top5.png'),
            model_filter=top_models,
            color_mapping=color_mapping  # Use consistent colors
        )
    
    # Create radar charts showing the top 3 models' performance
    top3_models = df.sort_values('overall_composite_score', ascending=False).head(3)['Model'].tolist()
    
    if len(celltype_metrics) >= 3:
        create_radar_chart(
            df,
            celltype_metrics,
            metric_labels,
            'Top 3 Models - Cell Type Metrics',
            'celltype_composite_score' if 'celltype_composite_score' in df.columns else 'overall_composite_score',
            os.path.join(output_dir, 'top3_celltype_metrics_radar.png'),
            model_filter=top3_models,
            color_mapping=color_mapping  # Use consistent colors
        )
    
    if len(lineage_metrics) >= 3:
        create_radar_chart(
            df,
            lineage_metrics,
            metric_labels,
            'Top 3 Models - Lineage Metrics',
            'lineage_composite_score' if 'lineage_composite_score' in df.columns else 'overall_composite_score',
            os.path.join(output_dir, 'top3_lineage_metrics_radar.png'),
            model_filter=top3_models,
            color_mapping=color_mapping  # Use consistent colors
        )
    
    # 2. Create bar charts for final scores
    score_types = [
        ('celltype_composite_score', 'Cell Type Composite Score'),
        ('lineage_composite_score', 'Lineage Composite Score'),
        ('overall_composite_score', 'Overall Composite Score')
    ]
    
    for score_col, title in score_types:
        if score_col not in df.columns:
            continue
            
        plt.figure(figsize=(12, 6))
        score_df = df.sort_values(score_col, ascending=False)
        
        # Color based on rank
        bar_colors = plt.cm.viridis(np.linspace(0, 1, len(score_df)))
        
        # Plot bars using display names for better readability
        if 'Display_Name' in score_df.columns:
            bars = plt.bar(score_df['Display_Name'], score_df[score_col], color=bar_colors)
        else:
            bars = plt.bar(score_df['Model'], score_df[score_col], color=bar_colors)
        
        # Add exact scores on top of bars
        for i, bar in enumerate(bars):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                    f"{score_df[score_col].iloc[i]:.1f}", 
                    ha='center', va='bottom')
        
        plt.title(f"Embryo Model Benchmark - {title}", size=15)
        plt.ylabel("Score (0-100)")
        plt.ylim(0, 105)
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, f"{score_col.replace('_score', '')}_barchart.png"), 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Create detailed heatmap
    create_detailed_heatmap(df, output_dir)

def create_detailed_heatmap(df, output_dir):
    """Create a detailed heatmap of all metrics and scores"""
    # Separate base metrics from scores
    base_metrics = [col for col in df.columns 
                   if col != 'Model' and col != 'Display_Name' and 
                   not col.endswith('_score') and 
                   not col.endswith('_rank') and
                   not col.endswith('_composite_score')]
    
    # Get aggregate scores
    agg_scores = [col for col in df.columns 
                 if col.endswith('_score') and 
                 not col.endswith('_composite_score')]
    
    # Get composite scores
    comp_scores = [col for col in df.columns 
                  if col.endswith('_composite_score')]
    
    # Select metrics to show in the heatmap
    all_metrics = base_metrics + agg_scores + comp_scores
    
    # Create a pivot table
    pivot_df = df.set_index('Model' if 'Display_Name' not in df.columns else 'Display_Name')[all_metrics].copy()
    
    # Sort by overall composite score
    if 'overall_composite_score' in pivot_df.columns:
        if 'Display_Name' in df.columns:
            pivot_df = pivot_df.reindex(df.sort_values('overall_composite_score', ascending=False)['Display_Name'])
        else:
            pivot_df = pivot_df.reindex(df.sort_values('overall_composite_score', ascending=False)['Model'])
    
    # Define more readable metric labels
    metric_labels = {
        # Individual metrics
        'celltype_certainty_percentage': 'CT Certainty %',
        'celltype_mean_certainty': 'CT Mean Certainty',
        'lineage_certainty_percentage': 'LIN Certainty %',
        'lineage_mean_certainty': 'LIN Mean Certainty',
        
        'celltype_marker_precision': 'CT Marker Precision',
        'celltype_marker_recall': 'CT Marker Recall',
        'celltype_marker_f1': 'CT Marker F1',
        'celltype_marker_jaccard': 'CT Marker Jaccard',
        
        'lineage_marker_precision': 'LIN Marker Precision',
        'lineage_marker_recall': 'LIN Marker Recall',
        'lineage_marker_f1': 'LIN Marker F1',
        'lineage_marker_jaccard': 'LIN Marker Jaccard',
        
        'celltype_correlation': 'CT Expression Corr',
        'lineage_correlation': 'LIN Expression Corr',
        
        'celltype_presence_percentage': 'CT Presence %',
        'lineage_presence_percentage': 'LIN Presence %',
        
        # Aggregated scores
        'celltype_certainty_score': 'CT Certainty Score',
        'lineage_certainty_score': 'LIN Certainty Score',
        
        'celltype_marker_overlap_score': 'CT Marker Score',
        'lineage_marker_overlap_score': 'LIN Marker Score',
        
        'celltype_correlation_score': 'CT Correlation Score',
        'lineage_correlation_score': 'LIN Correlation Score',
        
        'celltype_presence_score': 'CT Presence Score',
        'lineage_presence_score': 'LIN Presence Score',
        
        # Composite scores
        'celltype_composite_score': 'CT Composite',
        'lineage_composite_score': 'LIN Composite',
        'overall_composite_score': 'Overall Score'
    }
    
    # Create normalized version for colors
    norm_pivot = pivot_df.copy()
    for col in pivot_df.columns:
        if col in norm_pivot.columns:
            col_min = pivot_df[col].min()
            col_max = pivot_df[col].max()
            
            if col_max > col_min:
                norm_pivot[col] = (pivot_df[col] - col_min) / (col_max - col_min)
            else:
                norm_pivot[col] = 0
    
    # Create heatmap
    plt.figure(figsize=(18, len(pivot_df) * 0.8))
    
    # Plot heatmap
    ax = sns.heatmap(
        norm_pivot, annot=pivot_df.round(1), 
        fmt='.1f', cmap='viridis',
        cbar_kws={'label': 'Normalized Score'}
    )
    
    # Update tick labels
    ax.set_xticklabels(
        [metric_labels.get(m, m) for m in pivot_df.columns], 
        rotation=45, ha='right'
    )
    
    plt.title("Detailed Metrics Comparison Across Models", size=15)
    plt.tight_layout()
    
    # Save figure
    plt.savefig(os.path.join(output_dir, 'model_metrics_detailed_heatmap.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()


def main():
    # Set paths
    metrics_dir = '/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_embryomodel_umap_20250323/embryo_model_benchmarking_metrics'  # Directory containing metrics files
    output_dir = './benchmarking_results'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load all metrics data
    all_metrics = load_metrics_data(metrics_dir)
    
    if all_metrics is None or len(all_metrics) == 0:
        print("Error: No metrics data found")
        return
    
    print(f"Loaded metrics for {len(all_metrics)} models")
    print(f"Available metrics: {[col for col in all_metrics.columns if col not in ['Model', 'Display_Name']]}")
    
    # Calculate composite scores
    scored_metrics = calculate_composite_scores(all_metrics)
    
    # Create all plots
    create_plots(scored_metrics, output_dir)
    
    # Save comprehensive results to CSV
    result_file = os.path.join(output_dir, 'model_benchmark_final_scores.csv')
    scored_metrics.to_csv(result_file, index=False)
    
    # Print top models based on overall score
    print("\n=== Top Models by Overall Score ===")
    if 'overall_composite_score' in scored_metrics.columns:
        top_overall = scored_metrics.sort_values('overall_composite_score', ascending=False)
        for i, (_, row) in enumerate(top_overall.iterrows()):
            if i < 5:  # Show top 5
                display_name = row['Display_Name'] if 'Display_Name' in row else row['Model']
                print(f"{i+1}. {display_name} - Score: {row['overall_composite_score']:.2f}")
    
    print(f"\nBenchmark results saved to {output_dir}")


if __name__ == "__main__":
    main()