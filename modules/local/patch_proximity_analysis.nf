process Patch_Proximity_Analysis {
    label 'process_single'
    // publishDir "${params.output}/patch_proximity", mode: 'copy'
    container "ghcr.io/break-through-cancer/spacec:cpu"

    cpus 4
    memory 8.GB

    input:
    path input_h5ad

    output:
    path "adata_patch_proximity.h5ad", emit: 'proximity_h5ad'
    path "patch_proximity_celltype_tonsil.pdf"
    path "patch_proximity_celltype_tonsillitis.pdf"
    path "patch_proximity_CN_tonsil.pdf"
    path "patch_proximity_CN_tonsillitis.pdf"

    script:
    """
    6_patch_proximity.py $input_h5ad
    """
}