nextflow.enable.dsl = 2

// Include the workflow from workflows/spacec.nf
include { run_spacec } from './workflows/spacec.nf'

// Define the main workflow
workflow {
    // Ensure that the input file parameter is provided
    if (!params.input) {
        error "Please provide a file path using '--input' when running the pipeline"
    }

    // Call the run_spacec workflow from workflows/spacec.nf and pass the parameters to it
    run_spacec()
}