#!/usr/bin/env nextflow

// This is a nextflow for (re)producing the analyses for this guide thing
// we're working on.

// For temporary files
file("./tmp").mkdirs()
// For this script's reports
file("./reports").mkdirs()

script_make_degenerate_barcode_library_py = Channel
    .fromPath("scripts/make_degenerate-barcode-library.py")

script_sample_fasta_as_fasta = Channel
    .fromPath("scripts/sample_fasta_as_fasta.py")

script_eval_bartender_deviations = Channel
    .fromPath("scripts/eval_bartender_deviations.py")

script_eval_starcoded = Channel
    .fromPath("scripts/eval_starcode.py")

script_analysis_noise = Channel
    .fromPath("scripts/analysis_noise.R")

script_analysis_barcode_graph = Channel
    .fromPath("scripts/analysis_barcode_graph.R")

script_graph_barcode_library_py = Channel
    .fromPath("scripts/graph_barcode_library.py")


make_id_strings = {
        it["generation_id"] = it["barcode_pattern"]+
            "_"+it["base_mix"]+"_"+
            it["num_lineages"]+"lineages_"+it["num_barcodes_per"]+
            "barcodesper_fixedBarcodesPer"+it["fixed_barcodes"]+
            "_replicate"+it["replicate_id"] ; 
        it["sampling_id"] = it["cell_sample_size"] ; \
        it["sequencing_id"] = it["read_depth"]+"reads_"+
            it["abundance_distribution"]+"abundanceDist_"+
            it["error_rate"]+"error_"+
            it["percentage_subs_indels"]+"subsAndIndels_"+
            it["chimera_percentage"]+"chimeras_" // + it["slapchop_pattern"] ; 
        return it ; 
    } 

parameters = Channel.fromPath("scripts/parameters_simulations.tsv")
    .splitCsv(header: true, strip: true, sep: "\t")
    .map(make_id_strings)

process make_barcoded_libraries {
    input: 
        each script_make_degenerate_barcode_library_py
        each parameters
    output: 
        set val(parameters), 
            file({"barcodeLibrary_"+parameters["generation_id"]+".fasta"}) into barcoded_libraries
    shell:
    '''
    python3 !{script_make_degenerate_barcode_library_py} \
        !{parameters["barcode_pattern"]} \
        --mix !{parameters["base_mix"]} \
        --number-lineages !{parameters["num_lineages"]} \
        --barcodes-per-lineage !{parameters["num_barcodes_per"]} \
        --fixed-barcodes-per-clone !{parameters["fixed_barcodes"]} \
        --replicate !{parameters["replicate_id"]} \
        barcodeLibrary
    '''
}

( barcoded_libraries_graph, barcoded_libraries_sim) = barcoded_libraries.into(2)

process measure_barcode_distances {
    publishDir "tmp"
    input: 
        each script_graph_barcode_library_py
        set val(parameters), file(barcode_fasta) from barcoded_libraries_graph
    output:
        set val(parameters), file("barcodes.tsv"), file("barcodes.graph"), 
            file("barcodes_summary.txt") into barcode_graph
        file({"barcode_distances_"+parameters["generation_id"]+".tsv"}) into barcode_graph_tsv
    shell:
    '''
    python3 !{script_graph_barcode_library_py} \
        --store-graph !{barcode_fasta} barcodes
    cp barcodes.tsv barcode_distances_!{parameters["generation_id"]}.tsv
    '''
}


process collect_barcode_graphs_for_analysis {
    publishDir "tmp"
    input:
        file script_analysis_barcode_graph
        file(barcodes_tsv) from barcode_graph_tsv.collect()
    output:
        file("*.png") into summarize_graphs
    shell:
    '''
    Rscript !{script_analysis_barcode_graph} 
    '''
}

//
//
//
//
//

process sample_abundance_library {
    input: 
        each script_sample_fasta_as_fasta
        set val(parameters), file(barcode_fasta) from barcoded_libraries_sim
    output: 
        set val(parameters), file("barcode_library.fasta"), file("counted_abundances.txt") into sampled_libraries
    shell:
    '''
    python3 !{script_sample_fasta_as_fasta} !{barcode_fasta} !{parameters["cell_sample_size"]} barcode_library.fasta
    grep "> " barcode_library.fasta | sed 's/> //' | sort -n | uniq -c > counted_abundances.txt
    '''
}

process grinder {
    input:
        set val(parameters), file(fasta), file(abundances) from sampled_libraries
    output:
        set val(parameters), file(fasta), file(abundances), file("grinder-reads.fastq"), 
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
        -random_seed !{parameters["replicate"]} \
        -qual_levels 35 10 \
        -fastq_output 1 \
        -base_name grinder &> errorz
    '''
}


process slapchop {
    input:
        set val(parameters), file(fasta), file(abundances), file(grinder_fastq), 
            file(grinder_ranks) from grinder_output 
    output:
        set val(parameters), file(fasta), file(abundances), file("basename_pass.fastq"), 
            file(grinder_ranks) into slapchoped_fastq 
    shell: 
    '''
	python3 /slapchop.py !{grinder_fastq} basename \
        --bite-size 10000 --processes !{task.cpus}    \
        -o "GetLineageRead1: input > !{parameters["slapchop_pattern"]}" \
        --output-seq "barcode" \
        --output-id "input.id" \
        --verbose > stdout
    '''
}

process starcode {
    input:
        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
            file(grinder_ranks) from slapchoped_fastq 
    output:
        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
            file(grinder_ranks), file("starcoded") into starcoded
    shell:
    '''
    starcode --threads !{task.cpus} \
        -d 3 --cluster-ratio 10 \
        --input !{slapchop_pass} \
        --print-clusters --seq-id \
        --output starcoded
    '''
}

process eval_starcode {
    input:
        each script_eval_starcoded
        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
            file(grinder_ranks), file(starcoded) from starcoded
    output:
        file({parameters["generation_id"]+"_"+parameters["sampling_id"]+"_"+parameters["sequencing_id"]+"_clone_code_input_output.tsv"}) into eval_tsvs
    shell:
    '''
    python3 !{script_eval_starcoded} !{abundances} !{starcoded} "!{parameters["generation_id"]}_!{parameters["sampling_id"]}_!{parameters["sequencing_id"]}_clone_code_input_output.tsv"
    '''
}


process collect_tsvs_for_analysis {
    publishDir "tmp"
    input:
        file(each_tsv) from eval_tsvs.collect()
        each script_analysis_noise
    output:
        file("*.png") into summarize_noise
    shell:
    '''
    Rscript !{script_analysis_noise} 
    '''
}


// // Special trigger for `onComplete`. I copied this from documentation.
// // Some predefined variables. It somehow mails it. Cool.
//workflow.onComplete {
//    println "Pipeline completed at: $workflow.complete"
//    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
//
//    def subject = 'barnone run'
//    def recipient = email_address
//
//    ['mail', '-s', subject, recipient].execute() << """
//
//    Pipeline execution summary
//    ---------------------------
//    Completed at: ${workflow.complete}
//    Duration        : ${workflow.duration}
//    Success         : ${workflow.success}
//    workDir         : ${workflow.workDir}
//    exit status : ${workflow.exitStatus}
//    Error report: ${workflow.errorReport ?: '-'}
//    """
//}



