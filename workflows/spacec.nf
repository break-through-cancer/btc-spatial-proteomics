include { Cell_Neighborhood_Analysis } from '../modules/local/cell_neighborhood_analysis.nf'
include { Distance_Permutation_Analysis } from '../modules/local/distance_permutation_analysis.nf'
include { Patch_Proximity_Analysis } from '../modules/local/patch_proximity_analysis.nf'

workflow run_spacec {
    if (!params.input) {
        error "Please provide a file path using '--input' when running the pipeline"
    }
    
    files_ch = Channel.fromPath(params.input)

    result_ch = Cell_Neighborhood_Analysis(files_ch)
    
    result_ch_2 = Distance_Permutation_Analysis(result_ch.annotated_h5ad)
    
    Patch_Proximity_Analysis(result_ch_2.interaction_h5ad)
}
