process Cell_Neighborhood_Analysis {
    label 'process_single'
    // publishDir "${params.output}/cell_neighborhood", mode: 'copy'
    container "ghcr.io/break-through-cancer/spacec:cpu"

    cpus 4
    memory 8.GB

    input:
    path input_h5ad

    output:
    path "adata_nn_demo_annotated_cn.h5ad", emit: 'annotated_h5ad'
    path "catplot_spatial_plot.pdf"
    path "heatmap.pdf"
    path "normal_tonsil_CNMap.pdf"
    path "tonsilitis_CNMap.pdf"
    path "normal_tonsil_bc_proj.pdf"
    path "tonsilitis_bc_proj.pdf"

    script:
    """
    4_cna.py $input_h5ad
    """
}