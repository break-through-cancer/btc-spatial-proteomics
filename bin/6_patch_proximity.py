#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import spacec as sp
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

# Parse input arguments
parser = argparse.ArgumentParser(description="Patch Proximity Analysis")
parser.add_argument("input_h5ad", help="Path to input h5ad file")
# parser.add_argument("output_dir", help="Path to output directory")
args = parser.parse_args()

output_dir = './'

# Load AnnData object
adata = sc.read(args.input_h5ad)
print(adata)

# Define analysis parameters
pixel_to_um = 0.5085
distances = [5/pixel_to_um, 10/pixel_to_um, 15/pixel_to_um, 20/pixel_to_um, 25/pixel_to_um]
key_names = ['ppa_result_5', 'ppa_result_10', 'ppa_result_15', 'ppa_result_20', 'ppa_result_25']

# Run Patch Proximity Analysis
for distance, key_name in zip(distances, key_names):
    sp.tl.patch_proximity_analysis(
        adata, 
        region_column="unique_region",
        patch_column="CN_k20_n6_annot",
        group="Germinal Center",
        min_cluster_size=50,
        x_column="x", 
        y_column="y",
        radius=distance,
        edge_neighbours=3,
        key_name=key_name,
        plot=False
    )

# Save the processed AnnData object
adata.write(output_dir + "adata_patch_proximity.h5ad")
# print(f"Saved Patch Proximity analysis results to {output_file}")

# Function to plot results
def plot_results(adata, output_dir):
    """Generate and save donut plots for cell types and CNs around Germinal Centers."""

    key_names = ['ppa_result_5', 'ppa_result_10', 'ppa_result_15', 'ppa_result_20', 'ppa_result_25']
    distances = [5, 10, 15, 20, 25]

    # Check if analysis results exist
    if not any(k in adata.uns.keys() for k in key_names):
        print("No patch proximity results found in `adata.uns`. Ensure the analysis ran correctly.")
        return

    # Cell Type Donut Plots
    for condition in ["tonsil", "tonsillitis"]:
        pdf_filename = os.path.join(output_dir, f"patch_proximity_celltype_{condition}.pdf")

        fig, ax = plt.subplots(figsize=(10, 10))

        adata_filtered = adata[adata.obs['condition'] == condition]

        sp.pl.ppa_res_donut(
            adata_filtered,
            cat_col="celltype",
            key_names=key_names,
            radii=distances,
            unit="µm",
            add_guides=True,
            text="Cell Types Around Germinal Center",
            label_color="white",
            title=f"{condition.capitalize()} Patch Proximity Analysis"
        )

        plt.savefig(pdf_filename, format="pdf", dpi=300, bbox_inches="tight")
        plt.close(fig)  
        print(f"Saved: {pdf_filename}")

    # CN Donut Plots
    for condition in ["tonsil", "tonsillitis"]:
        pdf_filename = os.path.join(output_dir, f"patch_proximity_CN_{condition}.pdf")

        fig, ax = plt.subplots(figsize=(10, 10))

        adata_filtered = adata[adata.obs['condition'] == condition]

        sp.pl.ppa_res_donut(
            adata_filtered,
            cat_col="CN_k20_n6_annot",
            key_names=key_names,
            radii=distances,
            text="CNs Around Germinal Center",
            unit="µm",
            add_guides=True,
	    label_color="white",
            title=f"{condition.capitalize()} Patch Proximity Analysis"
        )
        
        plt.savefig(pdf_filename, format="pdf", dpi=300, bbox_inches="tight")
        plt.close(fig)  
        print(f"Saved: {pdf_filename}")

plot_results(adata, output_dir)
