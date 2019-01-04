#!/usr/bin/env nextflow

// This is a nextflow for (re)producing the analyses for this guide thing
// we're working on.

// For temporary files
file("./tmp").mkdirs()
// For this script's reports
file("./reports").mkdirs()

script_make_degenerate_barcode_library_py = Channel
    .fromPath("scripts/make_degenerate-barcode-library.py")

// Run parameters. Schema, by position:
//  1. barcode in IUPAC codes
//  2. number of lineages 
//  3. number of barcode clones (number of barcodes sampled)
run_parameters = Channel.from([
    [ "ATCGNNNNATCG" , 100, 1000 ]
    ])

process make_barcoded_libraries {
    input: 
        each script_make_degenerate_barcode_library_py
        each run_parameters
    output: 
        set run_parameters, file("barcode_library.fasta"), 
            file("barcode_library.yaml"), 
            file("barcode_library.csv") into barcoded_libraries
    shell:
    '''
    python3 !{script_make_degenerate_barcode_library_py} !{run_parameters[0]} barcode_library
    '''
}

