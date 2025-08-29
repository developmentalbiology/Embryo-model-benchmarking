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

# CS7 specific mapping (handles the exact names in CS7 file)
CS7_MODEL_NAME_MAPPING = {
    'Liu': 'corrected_processed_Liu_scPoli_annotated',
    'Rowan': 'corrected_processed_Rowan_scPoli_annotated', 
    'Hislop': 'corrected_processed_Hislop_scPoli_annotated',
    'Oldak': 'corrected_processed_Oldak_scPoli_annotated',
    'Pedroza': 'corrected_processed_Pedroza_scPoli_annotated',
    'Weatherbee': 'corrected_processed_Weatherbee_scPoli_annotated',
    'Ai_model': 'corrected_processed_Ai_model_scPoli_annotated',
    'model_Zheng_2019': 'corrected_processed_zheng_2019_scPoli_annotated',
    'model_zheng_2022': 'corrected_processed_zheng_2022_scPoli_annotated'
}

# Reverse mapping (for display in plots)
DISPLAY_NAME_MAPPING = {v: k for k, v in MODEL_NAME_MAPPING.items()}
# Add CS7 mappings to display mapping
for short_name, long_name in CS7_MODEL_NAME_MAPPING.items():
    if long_name not in DISPLAY_NAME_MAPPING:
        # Extract the base name for display
        if short_name.startswith('model_'):
            display_name = short_name.replace('model_', '').replace('_', '_')
        else:
            display_name = short_name
        DISPLAY_NAME_MAPPING[long_name] = display_name


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
        'pred_consistency': 'pred_consistency_comparison.csv',
        'celltype_mean_certainty': 'celltype_overall_mean_certainty_consistent_cells_comparison.csv',
        'lineage_mean_certainty': 'lineage_overall_mean_certainty_consistent_cells_comparison.csv',
        'celltype_correlation': 'celltype_correlation_comparison.csv',
        'lineage_correlation': 'lineage_correlation_comparison.csv',
        'celltype_presence': 'celltype_presence_comparison.csv',
        'lineage_presence': 'lineage_presence_comparison.csv',
        'celltype_marker_overlap': 'combined_celltype_metrics.csv',
        'lineage_marker_overlap': 'combined_lineage_metrics.csv',
        'cs7_correlation': 'stage_composition_correlation.csv'
    }
    
    all_metrics = {}
    loaded_models = set()
    
    # Load each metrics file
    for key, filename in files_to_load.items():
        filepath = os.path.join(metrics_dir, filename)
        if os.path.exists(filepath):
            print(f"Loading {key} metrics from {filepath}")
            try:
                df = pd.read_csv(filepath)
                print(f"  Loaded {len(df)} rows, columns: {df.columns.tolist()}")
            except Exception as e:
                print(f"  Error loading {filename}: {str(e)}")
                continue
            
            # Standardize Model column name if it's different
            model_col_candidates = ['Model', 'model', 'dataset', 'Dataset', 'Unnamed: 0']  # Add 'Unnamed: 0' for CS7 file
            model_col = None
            for candidate in model_col_candidates:
                if candidate in df.columns:
                    model_col = candidate
                    break

            if model_col is None:
                print(f"  Warning: No model column found in {filename}. Available columns: {df.columns.tolist()}")
                continue
                
            # Rename to standard 'Model' column
            if model_col != 'Model':
                df.rename(columns={model_col: 'Model'}, inplace=True)
                print(f"  Renamed '{model_col}' to 'Model' in {filename}")
                
            # Standardize model names to ensure consistency between files
            original_models = df['Model'].tolist()
            
            if key in ['celltype_marker_overlap', 'lineage_marker_overlap']:
                # For the marker overlap files, map from short names to long names
                df['Model'] = df['Model'].map(MODEL_NAME_MAPPING).fillna(df['Model'])
                # Ensure they have the _scPoli_annotated suffix
                df['Model'] = df['Model'].apply(lambda x: x + '_scPoli_annotated' if not x.endswith('_scPoli_annotated') else x)
            elif key == 'cs7_correlation':
                # For CS7 file, use the specialized mapping (already includes _scPoli_annotated)
                df['Model'] = df['Model'].map(CS7_MODEL_NAME_MAPPING).fillna(df['Model'])
            elif key == 'pred_consistency':
                # For consistency file, ensure _scPoli_annotated suffix
                df['Model'] = df['Model'].apply(lambda x: x + '_scPoli_annotated' if not x.endswith('_scPoli_annotated') else x)
            else:
                # For other files, ensure they have the _scPoli_annotated suffix if they don't already
                df['Model'] = df['Model'].apply(lambda x: x + '_scPoli_annotated' if not x.endswith('_scPoli_annotated') and not x.startswith('corrected_processed_') else x)
                # If they start with corrected_processed_ but don't have the suffix, add it
                df['Model'] = df['Model'].apply(lambda x: x + '_scPoli_annotated' if x.startswith('corrected_processed_') and not x.endswith('_scPoli_annotated') else x)
            
            standardized_models = df['Model'].tolist()
            if original_models != standardized_models:
                print(f"  Standardized model names in {filename}:")
                for orig, std in zip(original_models, standardized_models):
                    if orig != std:
                        print(f"    {orig} -> {std}")
            
            # Extract relevant columns
            if key == 'pred_consistency':
                metrics_df = df[['Model', 'consistency_percentage']].copy()
                # Create separate columns for celltype and lineage consistency
                metrics_df['celltype_pred_consistency'] = metrics_df['consistency_percentage']
                metrics_df['lineage_pred_consistency'] = metrics_df['consistency_percentage']
                # Keep the original column for backward compatibility
                metrics_df.rename(columns={'consistency_percentage': 'pred_consistency'}, inplace=True)
            
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
                    'Presence Score (%)'
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
                    'Presence Score (%)'
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
                
            elif key == 'cs7_correlation':
                # Handle CS7 correlation file - model names should already be mapped above
                print(f"  Processing CS7 file with columns: {df.columns.tolist()}")
                print(f"  Models before mapping: {df['Model'].tolist()}")
                
                if 'CS7_correlation' in df.columns:
                    metrics_df = df[['Model', 'CS7_correlation']].copy()
                    
                    # Scale correlation to 0-100 range (assuming correlation is between -1 and 1)
                    metrics_df['CS7_correlation'] = (metrics_df['CS7_correlation'] + 1) * 50
                    
                    print(f"  Successfully loaded CS7 correlation for models: {metrics_df['Model'].tolist()}")
                    print(f"  CS7 correlation values (scaled): {metrics_df['CS7_correlation'].tolist()}")
                else:
                    print(f"  ERROR: CS7_correlation column not found in {filename}")
                    print(f"  Available columns: {df.columns.tolist()}")
                    continue
                
            all_metrics[key] = metrics_df
            loaded_models.update(metrics_df['Model'].tolist())
            print(f"  Added {len(metrics_df)} models from {key}: {metrics_df['Model'].tolist()}")
        else:
            print(f"Warning: {filename} not found, skipping {key}")
    
    # List of all models from all files
    all_model_names = sorted(loaded_models)
    print(f"Found {len(all_model_names)} unique models across all metrics files")
    
    # Merge all metrics on 'Model' column
    combined_df = pd.DataFrame({'Model': all_model_names})
    
    for key, df in all_metrics.items():
        combined_df = pd.merge(combined_df, df, on='Model', how='left')
    
    # Fill NaNs with 0 for safer calculations
    combined_df.fillna(0, inplace=True)
    
    # Remove duplicate pred_consistency column (keep only celltype_pred_consistency and lineage_pred_consistency)
    if 'pred_consistency' in combined_df.columns:
        print("Removing duplicate pred_consistency column")
        combined_df = combined_df.drop(columns=['pred_consistency'])
    
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
            non_zero_count = (combined_df[col] != 0).sum()
            print(f"{col}: {combined_df[col].mean():.2f} (mean), {combined_df[col].min():.2f} (min), {combined_df[col].max():.2f} (max), {non_zero_count}/{len(combined_df)} non-zero")
    
    # Check for consistency data specifically
    consistency_cols = ['celltype_pred_consistency', 'lineage_pred_consistency', 'pred_consistency']
    for col in consistency_cols:
        if col in combined_df.columns:
            print(f"\nConsistency data check for {col}:")
            print(f"  Values: {combined_df[col].tolist()}")
            print(f"  Non-zero models: {combined_df[combined_df[col] != 0]['Model'].tolist()}")
    
    return combined_df

def calculate_composite_scores(df):
    """
    Calculate composite scores from raw metrics using equal weights for all metrics in radar charts
    
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
    
    # Define the metrics used in radar charts (equal weights)
    celltype_radar_metrics = [
        'celltype_mean_certainty',
        'celltype_pred_consistency',
        'celltype_presence_percentage',
        'celltype_marker_precision',
        'celltype_marker_jaccard',
        'celltype_marker_f1',
        'celltype_marker_recall',
        'celltype_correlation',
        'CS7_correlation'  # Added CS7_correlation to celltype composite
    ]
    
    lineage_radar_metrics = [
        'lineage_mean_certainty',
        'lineage_pred_consistency',
        'lineage_presence_percentage',
        'lineage_marker_precision',
        'lineage_marker_jaccard',
        'lineage_marker_f1',
        'lineage_marker_recall',
        'lineage_correlation'
    ]
    
    # Normalize metrics before calculating scores
    normalized_df = result_df.copy()
    
    # Normalize all metrics to [0, 1] range
    all_metrics = celltype_radar_metrics + lineage_radar_metrics
    for col in all_metrics:
        if col in result_df.columns:
            # Get min and max for normalization
            col_min = result_df[col].min()
            col_max = result_df[col].max()
            
            # Normalize to [0, 1] range
            if col_max > col_min:
                normalized_df[col] = (result_df[col] - col_min) / (col_max - col_min)
            else:
                normalized_df[col] = 0
    
    # Calculate celltype composite score (equal weights)
    available_celltype_metrics = [m for m in celltype_radar_metrics if m in normalized_df.columns]
    if available_celltype_metrics:
        result_df['celltype_composite_score'] = normalized_df[available_celltype_metrics].mean(axis=1) * 100
        print(f"Celltype composite score calculated from {len(available_celltype_metrics)} metrics: {available_celltype_metrics}")
    else:
        result_df['celltype_composite_score'] = 0
    
    # Calculate lineage composite score (equal weights)
    available_lineage_metrics = [m for m in lineage_radar_metrics if m in normalized_df.columns]
    if available_lineage_metrics:
        result_df['lineage_composite_score'] = normalized_df[available_lineage_metrics].mean(axis=1) * 100
        print(f"Lineage composite score calculated from {len(available_lineage_metrics)} metrics: {available_lineage_metrics}")
    else:
        result_df['lineage_composite_score'] = 0
    
    # Calculate overall composite score
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
    # Scale each metric so that max value = 100, but handle zero/missing data
    scaled_df = plot_df.copy()
    for metric in metrics:
        if metric in df.columns:
            max_value = df[metric].max()
            min_value = df[metric].min()
            
            # Check if all values are zero or very close to zero
            if max_value <= 0.001:
                print(f"Warning: All values for {metric} are zero or near zero. Setting to 0.")
                scaled_df[f"{metric}_scaled"] = 0
            elif max_value > 0:  # Normal scaling
                scaled_df[f"{metric}_scaled"] = (plot_df[metric] / max_value) * 100
            else:
                scaled_df[f"{metric}_scaled"] = 0
                
            # Debug info for consistency metrics
            if 'consistency' in metric.lower():
                print(f"Debug {metric}: min={min_value:.3f}, max={max_value:.3f}")
                print(f"  Original values: {plot_df[metric].tolist()}")
                print(f"  Scaled values: {scaled_df[f'{metric}_scaled'].tolist()}")
    
    # Additional check: ensure no NaN or infinite values
    for metric in metrics:
        scaled_col = f"{metric}_scaled"
        if scaled_col in scaled_df.columns:
            # Replace NaN and infinite values with 0
            scaled_df[scaled_col] = scaled_df[scaled_col].fillna(0)
            scaled_df[scaled_col] = scaled_df[scaled_col].replace([np.inf, -np.inf], 0)
    
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
        'celltype_pred_consistency': 'Cell Type\nConsistency %',
        'celltype_mean_certainty': 'Cell Type\nMean Certainty',
        'celltype_correlation': 'Cell Type\nExpression Correlation',
        'celltype_marker_recall': 'Cell Type\nMarker Recall',
        'celltype_marker_precision': 'Cell Type\nMarker Precision',
        'celltype_marker_f1': 'Cell Type\nMarker F1',
        'celltype_marker_jaccard': 'Cell Type\nMarker Jaccard',
        'CS7_correlation': 'CS7\nCorrelation',
        
        # Lineage metrics
        'lineage_presence_percentage': 'Lineage\nPresence %',
        'lineage_pred_consistency': 'Lineage\nConsistency %',
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
    
    # Modified metrics lists to include pred_consistency in both celltype and lineage
    celltype_metrics = [
        'celltype_mean_certainty',
        'celltype_pred_consistency',  # Added pred_consistency for celltype
        'celltype_presence_percentage',
        'celltype_marker_precision',
        'celltype_marker_jaccard',
        'celltype_marker_f1',
        'celltype_marker_recall',
        'celltype_correlation',
        'CS7_correlation'  # Added CS7_correlation to celltype radar
    ]
    
    lineage_metrics = [
        'lineage_mean_certainty',
        'lineage_pred_consistency',  # Added pred_consistency for lineage
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
    
    # Debug: Check data availability for each metric
    print("\nData availability check:")
    for metric in celltype_metrics + lineage_metrics:
        if metric in df.columns:
            non_zero_count = (df[metric] != 0).sum()
            print(f"  {metric}: {non_zero_count}/{len(df)} models have non-zero values")
            if non_zero_count == 0:
                print(f"    WARNING: {metric} has all zero values!")
        else:
            print(f"  {metric}: NOT FOUND in dataframe")
    
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
    

    # Create bar charts for final scores
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
    
    # Create detailed heatmap
    create_detailed_heatmap(df, output_dir)

def create_detailed_heatmap(df, output_dir):
    """Create a detailed heatmap of all metrics and scores"""
    # Separate base metrics from scores
    base_metrics = [col for col in df.columns 
                   if col not in ['Model', 'Display_Name'] and not col.endswith('_score') and not col.endswith('_rank')]
    
    score_metrics = [col for col in df.columns if col.endswith('_score')]
    
    # Create heatmap for base metrics
    if base_metrics:
        plt.figure(figsize=(16, 10))
        
        # Prepare data for heatmap
        heatmap_data = df[['Display_Name'] + base_metrics].set_index('Display_Name')
        
        # Normalize each column to 0-100 scale for better visualization
        normalized_data = heatmap_data.copy()
        for col in base_metrics:
            col_min = heatmap_data[col].min()
            col_max = heatmap_data[col].max()
            if col_max > col_min:
                normalized_data[col] = ((heatmap_data[col] - col_min) / (col_max - col_min)) * 100
            else:
                normalized_data[col] = 50  # If all values are the same, set to middle
        
        # Create heatmap
        sns.heatmap(normalized_data, 
                   annot=True, 
                   fmt='.1f',
                   cmap='RdYlBu_r',
                   center=50,
                   cbar_kws={'label': 'Normalized Score (0-100)'},
                   linewidths=0.5)
        
        plt.title('Embryo Model Benchmark - All Metrics Heatmap\n(Normalized to 0-100 scale)', size=15)
        plt.xlabel('Metrics')
        plt.ylabel('Models')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, 'all_metrics_heatmap.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    # Create heatmap for composite scores
    if score_metrics:
        plt.figure(figsize=(12, 8))
        
        # Prepare data for heatmap
        score_data = df[['Display_Name'] + score_metrics].set_index('Display_Name')
        
        # Create heatmap
        sns.heatmap(score_data, 
                   annot=True, 
                   fmt='.1f',
                   cmap='RdYlBu_r',
                   center=50,
                   cbar_kws={'label': 'Composite Score (0-100)'},
                   linewidths=0.5)
        
        plt.title('Embryo Model Benchmark - Composite Scores Heatmap', size=15)
        plt.xlabel('Score Types')
        plt.ylabel('Models')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, 'composite_scores_heatmap.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()

def main():
    """Main function to run the benchmarking analysis"""
    # Configuration
    metrics_dir = '/storage2/liuxiaodongLab/fanxueying/embryo_benchmarking_rebuttal/code/20250806_embryo_model_benchmarking/metrics_for_radar_plot'  # Change this to your metrics directory
    output_dir = './benchmarking_plots'
    
    print("="*60)
    print("EMBRYO MODEL BENCHMARKING ANALYSIS")
    print("="*60)
    
    # Load all metrics data
    print("\nStep 1: Loading metrics data...")
    df = load_metrics_data(metrics_dir)
    
    if df.empty:
        print("Error: No metrics data loaded. Please check your metrics directory and files.")
        return
    
    print(f"Loaded data for {len(df)} models with {len(df.columns)-2} metrics")
    
    # Calculate composite scores
    print("\nStep 2: Calculating composite scores...")
    df_with_scores = calculate_composite_scores(df)
    
    # Save the complete metrics table
    output_csv = os.path.join(output_dir, 'complete_benchmarking_metrics.csv')
    os.makedirs(output_dir, exist_ok=True)
    df_with_scores.to_csv(output_csv, index=False)
    print(f"Saved complete metrics table to {output_csv}")
    
    # Create all plots
    print("\nStep 3: Creating visualization plots...")
    create_plots(df_with_scores, output_dir)
    
    # Print summary
    print("\n" + "="*60)
    print("BENCHMARKING ANALYSIS COMPLETE")
    print("="*60)
    
    if 'overall_composite_score' in df_with_scores.columns:
        top_models = df_with_scores.nlargest(3, 'overall_composite_score')
        print("\nTop 3 models by overall composite score:")
        for i, (_, row) in enumerate(top_models.iterrows(), 1):
            display_name = row['Display_Name'] if 'Display_Name' in row else row['Model']
            print(f"{i}. {display_name}: {row['overall_composite_score']:.2f}")
    
    print(f"\nAll plots saved to: {output_dir}")
    print("Generated files:")
    print("- celltype_metrics_radar_all.png/pdf")
    print("- lineage_metrics_radar_all.png/pdf")
    print("- *_barchart.png (composite score comparisons)")
    print("- all_metrics_heatmap.png")
    print("- composite_scores_heatmap.png")
    print("- complete_benchmarking_metrics.csv")

if __name__ == "__main__":
    main()

