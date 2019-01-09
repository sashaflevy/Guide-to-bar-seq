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

script_plot_starcode = Channel
    .fromPath("scripts/plot_eval.R")

script_graph_barcode_library_py = Channel
    .fromPath("scripts/graph_barcode_library.py")

parameters = Channel.fromPath("scripts/parameters_simulations.tsv")
    .splitCsv(header: true, strip: true, sep: "\t")

process make_barcoded_libraries {
    input: 
        each script_make_degenerate_barcode_library_py
        each parameters
    output: 
        set val(parameters), 
            file({"barcodeLibrary_"+parameters["barcode_pattern"]+"_"+parameters["base_mix"].replaceAll(/,/,"-")+"_"+parameters["num_lineages"]+"lineages_"+parameters["num_barcodes_per"]+"barcodesper_fixedBarcodesPer"+parameters["fixed_barcodes"]+"_replicate"+parameters["replicate_id"]+".fasta"}) into barcoded_libraries
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
        set file("barcodes.tsv"), file("barcodes.graph") into barcode_graph
    shell:
    '''
    python3 !{script_graph_barcode_library_py} \
        --store-graph !{barcode_fasta} barcodes
    '''
}

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
        -random_seed 1234 \
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
        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
            file(grinder_ranks), file(starcoded), 
            file("eval_starcoded.txt") into eval_starcoded
    shell:
    '''
    python3 !{script_eval_starcoded} !{abundances} !{starcoded} eval_starcoded.txt
    '''
}

process plot_starcode {
    publishDir "tmp"
    input:
        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
            file(grinder_ranks), file(starcoded), 
            file(eval_starcoded) from eval_starcoded
        each script_plot_starcode
    output:
        set file({"var_"+parameters['barcode_pattern']+"_"+parameters['num_lineages']+"lineages_"+parameters['num_barcoded_clones']+"clones_replicate"+parameters['replicate_id']+"_per_clone.png"}),
            file({"var_"+parameters['barcode_pattern']+"_"+parameters['num_lineages']+"lineages_"+parameters['num_barcoded_clones']+"clones_replicate"+parameters['replicate_id']+"_per_code.png"}) into output
    shell:
    '''
    Rscript !{script_plot_starcode}
    mv var_per_code.png var_!{parameters['barcode_pattern']}_!{parameters['num_lineages']}lineages_!{parameters['num_barcoded_clones']}clones_replicate!{parameters['replicate_id']}_per_code.png
    mv var_per_clone.png var_!{parameters['barcode_pattern']}_!{parameters['num_lineages']}lineages_!{parameters['num_barcoded_clones']}clones_replicate!{parameters['replicate_id']}_per_clone.png
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



