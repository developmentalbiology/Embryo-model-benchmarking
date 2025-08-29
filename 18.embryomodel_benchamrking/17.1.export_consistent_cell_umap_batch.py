import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

def create_reanno_to_lineage_mapping():
    """
    Create the mapping from reanno_pred to lineage based on the provided R code.
    
    Returns:
    --------
    dict: mapping from reanno prediction to lineage
    """
    
    # Define the lineage to reanno mapping (converted from R list)
    lineages = {
        'TE_TrB': ['TE', 'CTB_1', 'CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2'],
        'epi': ['Epiblast_1', 'Epiblast_2', 'Epiblast_3', 'Ectoderm'],
        'Primitive.streak': ['Primitive.streak'],
        'NMP': ['Neuromesodermal.progenitor'],
        'Notochord': ['Notochord'],
        'PGC': ['PGC'],
        'ExE_endo': ['Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm'],
        'Amniotic_ecto': ['Amniontic.epi', 'Amniontic.ectoderm'],
        'neural_ecto': ['Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 
                        'Neural.ectoderm.midbrain', 'Spinal.cord'],
        'Endoderm': ['DE', 'DE_1', 'DE_2', 'Gut'],
        'meso_Exe.meso': ['Paraxial.mesoderm', 'Emergent.mesoderm', 'Pre-somatic.mesoderm', 'Somite', 
                          'Rostral.mesoderm', 'Lateral.plate.mesoderm_1', 'Lateral.plate.mesoderm_2', 
                          'Lateral.plate.mesoderm_3', 'Cardiac.mesoderm', 'Amniotic.mesoderm', 
                          'Exe.meso.progenitor', 'YS.mesoderm_1', 'YS.mesoderm_2'],
        'hemogenic': ['Hemogenic.endothelial.progenitor', 'Endothelium', 'Erythroid', 
                      'Primitive.megakaryocyte', 'Myeloid.progenitor']
    }
    
    # Create the reverse mapping from reanno to lineage
    reanno_to_lineage = {}
    
    for lineage_key, reanno_list in lineages.items():
        for reanno_value in reanno_list:
            # Check for mapping conflicts
            if reanno_value in reanno_to_lineage:
                warnings.warn(f"Warning: reanno value '{reanno_value}' is mapped to multiple lineages: "
                            f"'{reanno_to_lineage[reanno_value]}' and '{lineage_key}'. Using the last one.")
            reanno_to_lineage[reanno_value] = lineage_key
    
    return reanno_to_lineage

def add_reanno_pred_lineage_column(adata, reanno_to_lineage_mapping, 
                                   reanno_col='reanno_pred', 
                                   new_col='reanno_pred_lineage'):
    """
    Add a new column mapping reanno_pred to lineage.
    
    Parameters:
    -----------
    adata : AnnData object
        The AnnData object to modify
    reanno_to_lineage_mapping : dict
        Mapping from reanno prediction to lineage
    reanno_col : str
        Name of the reanno prediction column (default: 'reanno_pred')
    new_col : str
        Name of the new lineage column (default: 'reanno_pred_lineage')
    
    Returns:
    --------
    AnnData: Modified AnnData object with new lineage column
    """
    
    if reanno_col not in adata.obs.columns:
        raise ValueError(f"Column '{reanno_col}' not found in adata.obs")
    
    # Initialize the new column with NaN (equivalent to NA in R)
    adata.obs[new_col] = pd.NA
    
    # Create a pandas Series for efficient mapping
    reanno_values = adata.obs[reanno_col].astype(str)
    
    # Map reanno predictions to lineages
    mapped_lineages = reanno_values.map(reanno_to_lineage_mapping)
    
    # Assign the mapped values to the new column
    adata.obs[new_col] = mapped_lineages
    
    # Report mapping statistics
    mapped_count = mapped_lineages.notna().sum()
    unmapped_count = mapped_lineages.isna().sum()
    
    print(f"Mapping results:")
    print(f"  - Successfully mapped: {mapped_count} cells")
    print(f"  - Unmapped cells: {unmapped_count} cells")
    
    if unmapped_count > 0:
        unmapped_values = set(reanno_values[mapped_lineages.isna()])
        print(f"  - Unmapped reanno_pred values: {unmapped_values}")
    
    # Print lineage distribution
    print(f"\nLineage distribution:")
    lineage_counts = adata.obs[new_col].value_counts(dropna=False)
    for lineage, count in lineage_counts.items():
        print(f"  - {lineage}: {count} cells")
    
    return adata

def create_color_palettes():
    """
    Create color palettes for lineage and reanno predictions.
    
    Returns:
    --------
    tuple: (lineage_colors, ordered_reanno_labels, reanno_colors)
    """
    
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
    
    # godsnot_102 color palette
    godsnot_102 = [
        "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
        "#997D87", "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF",
        "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
        "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
        "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
        "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
        "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
        "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
        "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
        "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
        "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
        "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
        "#0089A3", "#CB7E98", "#A4E804", "#324E72"
    ]
    
    # Your ordered labels list for annotation visualization
    ordered_labels = [
        'TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
        'Epiblast_1','Epiblast_2','Epiblast_3','Ectoderm',
        'Amniontic.epi','Amniontic.ectoderm',
        'PGC',
        'Primitive.streak',
        'Neuromesodermal.progenitor',
        'Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 
        'Neural.ectoderm.midbrain','Spinal.cord',
        'Paraxial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 
        'Rostral.mesoderm', 'Lateral.plate.mesoderm_1',
        'Lateral.plate.mesoderm_2','Lateral.plate.mesoderm_3','Cardiac.mesoderm',
        'Amniotic.mesoderm','Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2',
        'Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm',
        'DE','Gut',
        'Notochord',
        'Hemogenic.endothelial.progenitor','Endothelium','Erythroid',
        'Primitive.megakaryocyte','Myeloid.progenitor'
    ]
    
    # Create reanno color palette mapping
    reanno_colors = {}
    for i, label in enumerate(ordered_labels):
        if i < len(godsnot_102):
            reanno_colors[label] = godsnot_102[i]
        else:
            # Fallback to cycling through colors if more labels than colors
            reanno_colors[label] = godsnot_102[i % len(godsnot_102)]
    
    return lineage_color, ordered_labels, reanno_colors

def plot_consistency_umaps(adata, file_stem, output_path, lineage_colors, ordered_reanno_labels, reanno_colors):
    """
    Plot UMAPs showing consistency between reanno_pred_lineage and lineage_pred.
    
    Parameters:
    -----------
    adata : AnnData
        The AnnData object
    file_stem : str
        File stem for naming output files
    output_path : Path
        Output directory path
    lineage_colors : dict
        Color mapping for lineages
    ordered_reanno_labels : list
        Ordered list of reanno labels for consistent coloring
    reanno_colors : dict
        Color mapping for reanno predictions
    """
    
    # Check if consistency can be calculated
    required_cols = ['reanno_pred_lineage', 'lineage_pred', 'reanno_pred']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"Warning: Missing columns {missing_cols} for consistency plotting")
        return
    
    # Define consistent cells (where reanno_pred_lineage == lineage_pred)
    consistent_cells = adata.obs['reanno_pred_lineage'] == adata.obs['lineage_pred']
    consistent_count = consistent_cells.sum()
    total_count = len(adata)
    
    print(f"Consistent cells: {consistent_count}/{total_count} ({consistent_count/total_count:.1%})")
    
    # Plot 1: Lineage predictions with consistency highlighting
    lineage_plot_colors = adata.obs['lineage_pred'].copy().astype(str)
    lineage_plot_colors[~consistent_cells] = 'Inconsistent'
    adata.obs['lineage_consistency_plot'] = lineage_plot_colors
    
    # Create color palette for lineage plot
    lineage_plot_palette = lineage_colors.copy()
    lineage_plot_palette['Inconsistent'] = '#CCCCCC'  # Grey for inconsistent
    
    plt.figure(figsize=(12, 8))
    sc.pl.umap(adata, color='lineage_consistency_plot', palette=lineage_plot_palette,
              title=f'Lineage Predictions (Consistent vs Inconsistent)',
              show=False, frameon=False)
    
    plt.savefig(output_path / f"{file_stem}_lineage_consistency_umap.pdf", 
               bbox_inches='tight', dpi=300)
    plt.close()
    
    # Plot 2: Reanno predictions with consistency highlighting
    reanno_plot_colors = adata.obs['reanno_pred'].copy().astype(str)
    reanno_plot_colors[~consistent_cells] = 'Inconsistent'
    adata.obs['reanno_consistency_plot'] = reanno_plot_colors
    
    # Set categorical levels to ALL ordered_reanno_labels for consistent coloring
    # This ensures consistent colors across all files regardless of which cell types are present
    all_categories = ordered_reanno_labels + ['Inconsistent']
    adata.obs['reanno_consistency_plot'] = pd.Categorical(
        adata.obs['reanno_consistency_plot'],
        categories=all_categories,
        ordered=True
    )
    
    # Create color list in the EXACT same order as the categorical levels
    # This ensures perfect alignment between categories and colors
    color_list = []
    for category in all_categories:
        if category == 'Inconsistent':
            color_list.append('#D3D3D3')  # Light grey for inconsistent
        else:
            color_list.append(reanno_colors[category])
    
    plt.figure(figsize=(12, 8))
    sc.pl.umap(adata, color='reanno_consistency_plot', palette=color_list,
              title=f'Reanno Predictions (Consistent vs Inconsistent)',
              show=False, frameon=False)
    
    plt.savefig(output_path / f"{file_stem}_reanno_consistency_umap.pdf", 
               bbox_inches='tight', dpi=300)
    plt.close()

def process_h5ad_files_add_lineage(folder_path, output_folder=None, 
                                   reanno_col='reanno_pred', 
                                   new_col='reanno_pred_lineage',
                                   generate_plots=True):
    """
    Process all h5ad files in a folder, add the reanno_pred_lineage column,
    generate consistency plots, and save updated files.
    
    Parameters:
    -----------
    folder_path : str
        Path to folder containing h5ad files
    output_folder : str, optional
        Path to output folder for plots (if None, uses 'consistency_plots' in folder_path)
    reanno_col : str
        Name of the reanno prediction column
    new_col : str
        Name of the new lineage column
    generate_plots : bool
        Whether to generate consistency plots
    """
    
    # Create the mapping and color palettes
    reanno_to_lineage_mapping = create_reanno_to_lineage_mapping()
    lineage_colors, ordered_reanno_labels, reanno_colors = create_color_palettes()
    
    print("Reanno to Lineage Mapping:")
    for reanno, lineage in reanno_to_lineage_mapping.items():
        print(f"  {reanno} -> {lineage}")
    print()
    
    # Get all h5ad files in the folder
    folder_path = Path(folder_path)
    h5ad_files = list(folder_path.glob("*.h5ad"))
    
    if not h5ad_files:
        print(f"No h5ad files found in {folder_path}")
        return
    
    # Set output path for plots
    if output_folder:
        plot_output_path = Path(output_folder)
    else:
        plot_output_path = folder_path / "consistency_plots"
    
    if generate_plots:
        plot_output_path.mkdir(exist_ok=True)
    
    print(f"Found {len(h5ad_files)} h5ad files")
    
    for file_path in h5ad_files:
        print(f"\nProcessing: {file_path.name}")
        
        try:
            # Load the h5ad file
            adata = sc.read_h5ad(file_path)
            
            # Check if reanno_pred column exists
            if reanno_col not in adata.obs.columns:
                print(f"Warning: Column '{reanno_col}' not found in {file_path.name}. Skipping.")
                continue
            
            # Add the lineage column
            adata = add_reanno_pred_lineage_column(adata, reanno_to_lineage_mapping, 
                                                   reanno_col, new_col)
            
            # Generate consistency plots if requested
            if generate_plots:
                plot_consistency_umaps(adata, file_path.stem, plot_output_path, 
                                     lineage_colors, ordered_reanno_labels, reanno_colors)
            
            # Save the modified file (overwrite original)
            adata.write_h5ad(file_path)
            print(f"Updated and saved: {file_path.name}")
            
        except Exception as e:
            print(f"Error processing {file_path.name}: {str(e)}")
            continue
    
    print(f"\nProcessing complete!")
    if generate_plots:
        print(f"Consistency plots saved to: {plot_output_path}")

def validate_mapping(folder_path, reanno_col='reanno_pred'):
    """
    Validate the mapping by showing all unique reanno_pred values across all files.
    
    Parameters:
    -----------
    folder_path : str
        Path to folder containing h5ad files
    reanno_col : str
        Name of the reanno prediction column
    """
    
    folder_path = Path(folder_path)
    h5ad_files = list(folder_path.glob("*.h5ad"))
    
    all_reanno_values = set()
    reanno_to_lineage_mapping = create_reanno_to_lineage_mapping()
    
    for file_path in h5ad_files:
        try:
            adata = sc.read_h5ad(file_path)
            if reanno_col in adata.obs.columns:
                file_reanno_values = set(adata.obs[reanno_col].astype(str).unique())
                all_reanno_values.update(file_reanno_values)
        except Exception as e:
            print(f"Error reading {file_path.name}: {str(e)}")
    
    print("All unique reanno_pred values found:")
    mapped_values = []
    unmapped_values = []
    
    for value in sorted(all_reanno_values):
        if value in reanno_to_lineage_mapping:
            mapped_values.append(value)
            print(f"  ✓ {value} -> {reanno_to_lineage_mapping[value]}")
        else:
            unmapped_values.append(value)
            print(f"  ✗ {value} -> NOT MAPPED")
    
    print(f"\nSummary:")
    print(f"  - Mapped values: {len(mapped_values)}")
    print(f"  - Unmapped values: {len(unmapped_values)}")
    
    if unmapped_values:
        print(f"\nUnmapped values that need attention: {unmapped_values}")

# Example usage
if __name__ == "__main__":
    # Specify your folder path here
    folder_path = "/storage2/liuxiaodongLab/fanxueying/mayanalysis/scPoli/20250801_embryomodel_export_plots/output"  # Change this to your actual folder path
    
    # First, validate the mapping to see what reanno_pred values exist
    print("=== VALIDATION ===")
    validate_mapping(folder_path)
    
    print("\n" + "="*50 + "\n")
    
    # Process all files, add the lineage column, generate consistency plots, and save
    print("=== PROCESSING ===")
    process_h5ad_files_add_lineage(
        folder_path=folder_path,
        output_folder=None,  # Will save plots to 'consistency_plots' subfolder
        generate_plots=True  # Generate consistency UMAP plots
    )