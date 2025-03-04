process Distance_Permutation_Analysis {
    label 'process_low'
    // publishDir "${params.output}/distance_permutation", mode: 'copy'
    container "ghcr.io/break-through-cancer/spacec:cpu"

    cpus 4
    memory 8.GB

    input:
    path input_h5ad

    output:
    path "adata_nn_demo_annotated_cn_5.h5ad", emit: 'interaction_h5ad'
    path "top_n_distances.pdf"
    path "dumbbell.pdf"
    path "tonsil_dist_graph.pdf"


    script:
    """
    5_distance_permutation.py $input_h5ad
    """
}