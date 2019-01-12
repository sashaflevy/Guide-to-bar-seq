#!/usr/bin/env nextflow

import groovyx.gpars.dataflow.operator.PoisonPill

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
    [ "bname":"six-base", "bpattern":"NNNNNN", // barcode name and pattern
        "bmix":"1-1-1-1", "clength":"2", // nucleotide mix and chimeria length (for grinder)
        "slapchopre":"(?P<barcode>[ATCG]{5,7}){e<=0}", 
        "bartenderpat":"" ],
    )
    .concat(Channel.from(PoisonPill.instance))

library_parameters = Channel.from(
    [ "bname":"six-base", "lineages":"10", "codesper":"1.0", "fixedper":"True" ], // name, number lineages, number codes per lineage, fixed number or poisson
    )
    .combine( Channel.from( ["librep":"a"], ["librep":"b"], ["librep":"c"] ) )
    .map{ it[0] + it[1] } 
    .concat(Channel.from(PoisonPill.instance))

sampling_parameters = Channel.from(
    [ "csample":"100" ], // cell sample size
    )
    .combine( Channel.from( ["samrep":"a"] ) )
    .map{ it[0] + it[1] } 

sequencing_parameters = Channel.from(
    [ "rdepth":"1000", "serror":"uniform0.0", "subrate":"80", "indelrate":"20", "crate":"01" ], // read depth, single base error distribution, sub rate percentage for singles, indel rate percentage for singles, chimera rate in percentage
    )
    .combine( Channel.from( ["seqrep":"a"] ) )
    .map{ it[0] + it[1] } 
    .concat(Channel.from(PoisonPill.instance))

//
//
//
 
make_safe_name = { 
    it.each{ key, value -> key+'='+value }
        .collect().join('_')
        .replaceAll("[\\(\\)<>]","").replaceAll(",",".") 
        .replaceAll("[aeiou]","")
    } 

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
    python3 mbl.py !{parameters["bpattern"]} \
        --mix !{parameters["bmix"]} \
        --number-lineages !{parameters["lineages"]} \
        --barcodes-per-lineage !{parameters["codesper"]} \
        --fixed-barcodes-per-clone !{parameters["fixedper"]} \
        --replicate !{parameters["librep"]} \
        barcodeLibrary_!{make_safe_name(parameters)}
    '''
}

// What is this bullshit? Why, this is some tricky bits to hide variables from
// nextflow, because it seems to copy links and references and not values 
// between channels of values. Sounds slick, but that means you can't reliably
// fork. So, a hack:

( copy1, copy2 ) = barcoded_libraries
    .into(2)

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
    python3 gbl.py --store-graph !{barcode_fasta} \
        barcode_distances_!{make_safe_name(parameters)}
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
// !{make_safe_name(parameters)}

process sample_abundance_library {
    input: 
        each file("sff.py") from Channel.fromPath("scripts/sample_fasta_as_fasta.py")
        set val(parameters), file(barcode_fasta) from barcoded_libraries_to_sample
    output: 
        set val(parameters), file("barcode_library.fasta"), file("counted_abundances.txt") into sampled_libraries_from_sampler
    shell:
    '''
    python3 sff.py !{barcode_fasta} !{parameters["csample"]} barcode_library.fasta
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
        set val(parameters), file(fasta), file(abundances), 
            file("grinder-reads.fastq") into sequenced_fastq
    shell:
    '''
    grinder -reference_file !{fasta} \
        -total_reads !{parameters["rdepth"]} \
        -read_dist 150 \
        -unidirectional 1 \
        -mutation_dist !{parameters["single_error_rate"]} \
        -mutation_ratio !{parameters["subrate"]+" "+parameters["indelrate"]} \
        -chimera_perc !{parameters["crate"]} \
        -ck !{parameters["clength"]} \
        -abundance_model uniform \
        -random_seed !{parameters["seqrep"]} \
        -qual_levels 35 10 \
        -fastq_output 1 \
        -base_name grinder
    '''
}

//
//
//

slapchop_input = Channel.create()
sequenced_fastq.into(slapchop_input)

process slapchop {
    input:
        set val(parameters), file(fasta), file(abundances), 
            file(grinder_fastq) from slapchop_input
    output:
        set val(parameters), file(fasta), file(abundances), 
            file("basename_pass.fastq") into slapchopd_fastq
    shell: 
    '''
	python3 /slapchop.py !{grinder_fastq} basename \
        --bite-size 10000 --processes !{task.cpus}    \
        -o "GetLineageRead1: input > !{parameters["slapchopre"]}" \
        --output-seq "barcode" \
        --output-id "input.id" \
        --verbose
    '''
}

process starcode {
    input:
        set val(parameters), file(fasta), file(abundances), 
            file(slapchop_pass) from slapchopd_fastq
    output:
        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
            file("starcoded") into starcoded
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
        each file("es.py") from Channel.fromPath("scripts/eval_starcode.py")
        set val(parameters), file(fasta), file(abundances), file(slapchop_pass), 
            file(starcoded) from starcoded
    output:
        file("*tsv") into eval_tsvs
    shell:
    '''
    python3 es.py !{abundances} !{starcoded} \
        counts_!{make_safe_name(parameters)}.tsv
    '''
}

// script_eval_bartender_deviations_py = Channel.fromPath("scripts/eval_bartender_deviations.py")

process collect_tsvs_for_analysis {
    publishDir "tmp"
    input:
        file(each_tsv) from eval_tsvs.collect()
        each file("an.R") from Channel.fromPath("scripts/analysis_noise.R")
    output:
        file("*.png") into summarize_noise
    shell:
    '''
    Rscript an.R
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



