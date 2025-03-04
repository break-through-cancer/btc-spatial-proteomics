#!/usr/bin/env python3

import argparse
import os

# import spacec first
import spacec as sp

#import standard packages
import pandas as pd
import scanpy as sc

# silencing warnings
import warnings
warnings.filterwarnings('ignore')

# Initialize and parse input arguments
parser = argparse.ArgumentParser(description="Take positional args")
parser.add_argument("input_h5ad")
args = parser.parse_args()

output_dir = './'

sc.settings.set_figure_params(dpi=80, facecolor='white')

# Loading the anndata from notebook 3 [cell type or cluster annotation is necessary for the step]
adata = sc.read(args.input_h5ad)
print(adata)

# compute the potential interactions
distance_pvals, results_dict = sp.tl.identify_interactions(
    adata = adata, # AnnData object
    cellid = "index", # column that contains the cell id (set index if the cell id is the index of the dataframe)
    x_pos = "x", # x coordinate column
    y_pos = "y", # y coordinate column
    cell_type = "celltype", # column that contains the cell type information
    region = "unique_region", # column that contains the region information
    num_iterations=1000, # number of iterations for the permutation test
    num_cores=12,  # number of CPU threads to use
    min_observed = 10, # minimum number of observed interactions to consider a cell type pair
    comparison = 'condition', # column that contains the condition information we want to compare
    distance_threshold=128) # distance threshold in px (20 µm)

# the results_dict contains the results of the permutation test as well as the observed and shuffled distances
results_dict.keys()

# the distance_pvals contains the p-values for each cell type pair and is automatically added to the adata.uns
adata.uns['triDist']

# save adata
adata.write(output_dir + "adata_nn_demo_annotated_cn_5.h5ad")


# Filter for most significant results
distance_pvals_filt = sp.tl.remove_rare_cell_types(adata, 
    distance_pvals, 
    cell_type_column="celltype", 
    min_cell_type_percentage=1
)
print(distance_pvals_filt.shape)
print(distance_pvals_filt.head()) 

# Identify significant cell-cell interactions
# dist_table_filt is a simplified table used for plotting
# dist_data_filt contains the filtered raw data with more information about the pairs
#  The function outputs two dataframes:  and dist_data_filt that contains all filtered interactions and  dist_table_filt that contains a table for all interactions that show a significant value in both tissues
dist_table_filt, dist_data_filt = sp.tl.filter_interactions(
    distance_pvals = distance_pvals_filt,
    pvalue = 0.05,
    logfold_group_abs = 0.1,
    comparison = 'condition')

print(dist_table_filt.shape)
print(dist_table_filt.head()) 


# Visualization
sp.pl.plot_top_n_distances(
    dist_table_filt,
    dist_data_filt,
    n=5,
    colors=None,
    dodge=False,
    savefig=True,
    output_fname="top_n_distances",
    output_dir=output_dir,
    figsize=(5, 5),
    unit="px",
    errorbars=True,
)

sp.pl.dumbbell(data = dist_table_filt,
    figsize=(8,12),
    colors = ['#DB444B', '#006BA2'],
    savefig=True,
    output_fname="dumbbell",
    output_dir=output_dir
)

# Adding filtering to account for low cell count in sample data
distance_pvals = distance_pvals.dropna()
# Identify pairs that appear in both conditions
pair_counts = distance_pvals.groupby(['celltype1', 'celltype2'])['condition'].nunique()
valid_pairs = pair_counts[pair_counts > 1].index

# Filter the DataFrame to keep only valid pairs
filtered_df = distance_pvals[distance_pvals.set_index(['celltype1', 'celltype2']).index.isin(valid_pairs)]

sp.pl.distance_graph(dist_table = dist_data_filt, # the (filtered) distance data table you want to plot 
    distance_pvals = filtered_df,
    condition_pair=['tonsil', 'tonsillitis'],
    node_size=1600, font_size=6,
    palette=None,
    dpi = 600,
    savefig=True,
    output_fname="tonsil",
    output_dir=output_dir
)
