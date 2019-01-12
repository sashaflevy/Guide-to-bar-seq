#!/usr/bin/env nextflow

// This is a nextflow for analyzing the contribution of experimental design to
// theoretical measurements of DNA barcode libraries.

// For temporary files
file("./tmp").mkdirs()
// For this script's reports
file("./reports").mkdirs()

//
//
//

barcode_parameters = Channel.from(
    [ "barcode_name":"six-base", "barcode_pattern":"NNNNNN", "base_mix":"1-1-1-1", "chimera_length":"2", ] 
    )

// "slapchop_pattern":"(?P<barcode>[ATCG]{5,7}){e<=0}", "bartender_pattern":"" ],

library_parameters = Channel.from(
    [ "barcode_name":"six-base", "num_lineages":"10", "num_barcodes_per":"1.0", "fixed_barcodes":"True" ],
    )
    .combine( Channel.from( ["lib_rep":"a"], ["lib_rep":"b"], ["lib_rep":"c"] ) )
    .map{ it[0] + it[1] } 

sampling_parameters = Channel.from(
    [ "cell_sample_size":"100" ],
    )
    .combine( Channel.from( ["sample_rep":"a"] ) )
    .map{ it[0] + it[1] } 

sequencing_parameters = Channel.from(
    [ "read_depth":"1000", "single_error_rate":"uniform0.0", "percentage_substitutions":"80", "percentage_indels":"20", "chimera_percentage":"01" ],
    )
    .combine( Channel.from( ["seq_rep":"a"] ) )
    .map{ it[0] + it[1] } 

//
//
//

parameters_making_libraries = barcode_parameters.combine(library_parameters) 
    .map{ it[0] + it[1] }

process make_barcoded_libraries {
    input: 
        each file("mbl.py") from Channel.fromPath("scripts/make_degenerate-barcode-library.py")
        each parameters from parameters_making_libraries
    output: 
        set val(parameters), file("barcodeLibrary*.fasta") into barcoded_libraries
    shell:
    '''
    python3 mbl.py !{parameters["barcode_pattern"]} \
        --mix !{parameters["base_mix"]} \
        --number-lineages !{parameters["num_lineages"]} \
        --barcodes-per-lineage !{parameters["num_barcodes_per"]} \
        --fixed-barcodes-per-clone !{parameters["fixed_barcodes"]} \
        --replicate !{parameters["replicate"]} \
        barcodeLibrary_!{parameters.each{ key, value -> key+'='+value }.collect().join('_')}
    '''
}

// What is this bullshit? Why, this is some tricky bits to hide variables from
// nextflow, because it seems to copy links and references and not values 
// between channels of values. Sounds slick, but that means you can't reliably
// fork. So, a hack:

( copy1, copy2 ) = barcoded_libraries.into(2)

barcoded_libraries_to_analyze = Channel.create()
copy1.subscribe{ barcoded_libraries_to_analyze.bind(it)  }

barcoded_libraries_to_sample_pre = Channel.create()
copy2.subscribe{ barcoded_libraries_to_sample_pre.bind(it)  }
barcoded_libraries_to_sample = barcoded_libraries_to_sample_pre
    .combine(sampling_parameters) 
    .map{ 
        it[0] = it[0] + it[2] ; 
        it.removeLast() ;
        it 
        }

//barcoded_libraries_to_analyze.subscribe{ println it }
//barcoded_libraries_to_sample.subscribe{ println it }

process measure_barcode_distances {
    publishDir "tmp"
    input: 
        each file("gbl.py") from Channel.fromPath("scripts/graph_barcode_library.py")
        set val(parameters), file(barcode_fasta) from barcoded_libraries_to_analyze
    output:
        file("barcode_distances_*.tsv") into barcode_graph_tsv
    shell:
    '''
    echo barcode_distances_!{parameters}
    python3 gbl.py --store-graph !{barcode_fasta} \
        barcode_distances_!{parameters.each{ key, value -> key+'='+value }.collect().join('_')}
    '''
}

process analyze_barcode_graphs {
    publishDir "tmp"
    input:
        each file("abg.R") from Channel.fromPath("scripts/analysis_barcode_graph.R")
        file(barcodes_tsv) from barcode_graph_tsv.collect()
    output:
        file("*.png") into summarize_graphs
    shell:
    '''
    Rscript abg.R
    '''
}

//
//
//
//
//


// !{parameters.each{ key, value -> key+'='+value }.collect().join('_')}

process sample_abundance_library {
    input: 
        each file("sff.py") from Channel.fromPath("scripts/sample_fasta_as_fasta.py")
        set val(parameters), file(barcode_fasta) from barcoded_libraries_to_sample
    output: 
        set val(parameters), file("barcode_library.fasta"), file("counted_abundances.txt") into sampled_libraries_from_sampler
    shell:
    '''
    python3 sff.py !{barcode_fasta} !{parameters["cell_sample_size"]} barcode_library.fasta
    grep "> " barcode_library.fasta | sed 's/> //' | sort -n | uniq -c > counted_abundances.txt
    '''
}

sampled_libraries_for_grinder = sampled_libraries_from_sampler
    .combine(sequencing_parameters) 
    .map{ 
        it[0] = it[0] + it[3] ; 
        it.removeLast() ;
        it 
        }

process grinder {
    input:
        set val(parameters), file(fasta), file(abundances) from sampled_libraries_for_grinder
    output:
        set val(parameters), file(fasta), file(abundances), file("grinder-reads.fastq"), 
            file("grinder-ranks.txt") into grinder_output 
    shell:
    '''
    grinder -reference_file !{fasta} \
        -total_reads !{parameters["read_depth"]} \
        -read_dist 150 \
        -unidirectional 1 \
        -mutation_dist !{parameters["single_error_rate"]} \
        -mutation_ratio !{parameters["percentage_substitutions"]+" "+parameters["percentage_indels"]} \
        -chimera_perc !{parameters["chimera_percentage"]} \
        -ck !{parameters["chimera_length"]} \
        -abundance_model uniform \
        -random_seed !{parameters["seq_rep"]} \
        -qual_levels 35 10 \
        -fastq_output 1 \
        -base_name grinder
    '''
}

script_eval_bartender_deviations_py = Channel.fromPath("scripts/eval_bartender_deviations.py")
script_eval_starcoded_py = Channel.fromPath("scripts/eval_starcode.py")
script_analysis_noise_r = Channel.fromPath("scripts/analysis_noise.R")

//
//
//



//process slapchop {
//    input:
//        set val(parameters), file(fasta), file(abundances), file(grinder_fastq), 
//            file(grinder_ranks) from grinder_output 
//    output:
//        set val(parameters), file(fasta), file(abundances), file("basename_pass.fastq"), 
//            file(grinder_ranks) into slapchoped_fastq 
//    shell: 
//    '''
//	python3 /slapchop.py !{grinder_fastq} basename \
//        --bite-size 10000 --processes !{task.cpus}    \
//        -o "GetLineageRead1: input > !{parameters["slapchop_pattern"]}" \
//        --output-seq "barcode" \
//        --output-id "input.id" \
//        --verbose > stdout
//    '''
//}
//
//process starcode {
//    input:
//        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
//            file(grinder_ranks) from slapchoped_fastq 
//    output:
//        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
//            file(grinder_ranks), file("starcoded") into starcoded
//    shell:
//    '''
//    starcode --threads !{task.cpus} \
//        -d 3 --cluster-ratio 10 \
//        --input !{slapchop_pass} \
//        --print-clusters --seq-id \
//        --output starcoded
//    '''
//}
//
//process eval_starcode {
//    input:
//        file script_eval_starcoded_py
//        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
//            file(grinder_ranks), file(starcoded) from starcoded
//    output:
//        file({parameters["generation_id"]+"_"+parameters["sampling_id"]+"_"+parameters["sequencing_id"]+"_clone_code_input_output.tsv"}) into eval_tsvs
//    shell:
//    '''
//    python3 !{script_eval_starcoded} !{abundances} !{starcoded} "!{parameters["generation_id"]}_!{parameters["sampling_id"]}_!{parameters["sequencing_id"]}_clone_code_input_output.tsv"
//    '''
//}
//
//
//process collect_tsvs_for_analysis {
//    publishDir "tmp"
//    input:
//        file(each_tsv) from eval_tsvs.collect()
//        file script_analysis_noise_r
//    output:
//        file("*.png") into summarize_noise
//    shell:
//    '''
//    Rscript !{script_analysis_noise} 
//    '''
//}
//
//
//// // Special trigger for `onComplete`. I copied this from documentation.
//// // Some predefined variables. It somehow mails it. Cool.
////workflow.onComplete {
////    println "Pipeline completed at: $workflow.complete"
////    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
////
////    def subject = 'barnone run'
////    def recipient = email_address
////
////    ['mail', '-s', subject, recipient].execute() << """
////
////    Pipeline execution summary
////    ---------------------------
////    Completed at: ${workflow.complete}
////    Duration        : ${workflow.duration}
////    Success         : ${workflow.success}
////    workDir         : ${workflow.workDir}
////    exit status : ${workflow.exitStatus}
////    Error report: ${workflow.errorReport ?: '-'}
////    """
////}
//
//
//
