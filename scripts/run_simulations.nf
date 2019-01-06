#!/usr/bin/env nextflow

// This is a nextflow for (re)producing the analyses for this guide thing
// we're working on.

// For temporary files
file("./tmp").mkdirs()
// For this script's reports
file("./reports").mkdirs()

script_make_degenerate_barcode_library_py = Channel
    .fromPath("scripts/make_degenerate-barcode-library.py")

script_sample_yaml_as_fasta = Channel
    .fromPath("scripts/sample_yaml_as_fasta.py")

script_eval_bartender_deviations = Channel
    .fromPath("scripts/eval_bartender_deviations.py")

parameters = Channel.fromPath("scripts/parameters_simulations.csv")
    .splitCsv(header: true)

process make_barcoded_libraries {
    publishDir "tmp"
    input: 
        each script_make_degenerate_barcode_library_py
        each parameters
    output: 
        set val(parameters), file("barcode_library.yaml"), file("barcode_library.fasta") into barcoded_libraries
    shell:
    '''
    python3 !{script_make_degenerate_barcode_library_py} \
        --number-lineages !{parameters["num_lineages"]} \
        --number-barcoded-clones !{parameters["num_barcoded_clones"]} \
        --replicate !{parameters["replicate_id"]} \
        !{parameters["barcode_pattern"]} barcodeLibrary
    mv *yaml barcode_library.yaml
    mv *fasta barcode_library.fasta
    '''
}

process sample_abundance_library {
    publishDir "tmp"
    input: 
        each script_sample_yaml_as_fasta
        set val(parameters), file(barcode_yaml), file(barcode_fasta) from barcoded_libraries
    output: 
        set val(parameters), file(barcode_yaml), file("barcode_library.fasta"), file("abundances.txt") into sampled_libraries
    shell:
    '''
    python3 !{script_sample_yaml_as_fasta} !{barcode_yaml} !{parameters["cell_sample_size"]} barcode_library.fasta abundances.txt
    '''
}

process grinder {
    publishDir "tmp"
    input:
        set val(parameters), file(yaml), file(fasta), file(abundances) from sampled_libraries
    output:
        set val(parameters), file(yaml), file(abundances), file("grinder-reads.fastq"), 
            file("grinder-ranks.txt") into grinder_output 
    shell:
    '''
    grinder -reference_file !{fasta} \
        -total_reads !{parameters["read_depth"]} \
        -read_dist 150 \
        -unidirectional 1 \
        -mutation_dist !{parameters["error_rate"]} \
        -mutation_ratio !{parameters["percentage_subs_indels"]} \
        -chimera_perc !{parameters["chimera_percentage"]} \
        -abundance_model !{parameters["abundance_distribution"]} \
        -random_seed 1234 \
        -qual_levels 35 10 \
        -fastq_output 1 \
        -base_name grinder &> errorz
    '''
}

process bartender_extract {
    publishDir "tmp"
    input:
        set parameters, file(yaml), file(abundances), file(input_fastq), file(input_ranks) from grinder_output 
    output:
        set parameters, file(yaml), file(abundances), file(input_ranks), file("extracted_barcode.txt") into bartender_extracted_codes 
    shell:
    '''
    bartender_extractor_com -f !{input_fastq} -o extracted -q ? \
        -p !{parameters["bartender_pattern"]} -m 2
    '''
}

process bartender_quant {
    publishDir "tmp"
    input:
        set parameters, file(yaml), file(abundances), file(input_ranks), file(extracted_barcodes) from bartender_extracted_codes 
    output:
        set parameters, file(yaml), file(abundances), file(input_ranks),
            file("barcode_quant_barcode.csv"),
            file("barcode_quant_cluster.csv"),
            file("barcode_quant_quality.csv") into bartender_quantifications
    shell:
    '''
    bartender_single_com -f !{extracted_barcodes} \
        -o barcode_quant -d 3
    '''
}

process bartender_eval {
    publishDir "tmp"
    input:
        set parameters, file(yaml), file(abundances), file(input_ranks),
            file(bartender_quant_barcode),
            file(bartender_quant_cluster),
            file(bartender_quant_quality) from bartender_quantifications
        each script_eval_bartender_deviations
    output:
        set parameters, file("bartender_deviations.csv") into bartender_deviations
    shell:
    '''
    python3 !{script_eval_bartender_deviations} !{yaml} !{abundances} !{bartender_quant_cluster} !{bartender_quant_cluster} "!{parameters['barcode_pattern']}"
    '''
}
