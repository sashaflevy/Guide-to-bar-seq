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
            file({"barcodeLibrary_"+parameters['barcode_pattern']+"_"+parameters['num_lineages']+"lineages_"+parameters['num_barcoded_clones']+"clones_replicate"+parameters['replicate_id']+".yaml"}),
            file({"barcodeLibrary_"+parameters['barcode_pattern']+"_"+parameters['num_lineages']+"lineages_"+parameters['num_barcoded_clones']+"clones_replicate"+parameters['replicate_id']+".fasta"}) into barcoded_libraries
    shell:
    '''
    python3 !{script_make_degenerate_barcode_library_py} \
        --number-lineages !{parameters["num_lineages"]} \
        --number-barcoded-clones !{parameters["num_barcoded_clones"]} \
        --replicate !{parameters["replicate_id"]} \
        !{parameters["barcode_pattern"]} barcodeLibrary
    '''
}

( barcoded_libraries_graph, barcoded_libraries_sim) = barcoded_libraries.into(2)

process measure_barcode_distances {
    publishDir "tmp"
    input: 
        each script_graph_barcode_library_py
        set val(parameters), file(barcode_yaml), file(barcode_fasta) from barcoded_libraries_graph
    output:
        set file("barcode_graph.tsv"), file("barcode_graph.graph") into barcode_graph
    shell:
    '''
    python3 !{script_graph_barcode_library_py} \
        --store-graph --yaml \
        !{barcode_yaml} barcode_graph
    '''
}

process sample_abundance_library {
    input: 
        each script_sample_yaml_as_fasta
        set val(parameters), file(barcode_yaml), file(barcode_fasta) from barcoded_libraries_sim
    output: 
        set val(parameters), file(barcode_yaml), file("barcode_library.fasta"), file("counted_abundances.txt") into sampled_libraries
    shell:
    '''
    python3 !{script_sample_yaml_as_fasta} !{barcode_yaml} !{parameters["cell_sample_size"]} barcode_library.fasta abundances.txt
    cat abundances.txt | sort | uniq -c > counted_abundances.txt
    '''
}

process grinder {
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


process slapchop {
    input:
        set val(parameters), file(yaml), file(abundances), file(grinder_fastq), 
            file(grinder_ranks) from grinder_output 
    output:
        set val(parameters), file(yaml), file(abundances), file("basename_pass.fastq"), 
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
        set val(parameters), file(yaml), file(abundances), file(slapchop_pass), 
            file(grinder_ranks) from slapchoped_fastq 
    output:
        set val(parameters), file(yaml), file(abundances), file(slapchop_pass), 
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
        set val(parameters), file(yaml), file(abundances), file(slapchop_pass), 
            file(grinder_ranks), file(starcoded) from starcoded
    output:
        set val(parameters), file(yaml), file(abundances), file(slapchop_pass), 
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
        set val(parameters), file(yaml), file(abundances), file(slapchop_pass), 
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

//
//tmp/sample%_clustered_lineage_tags.starclust: \
//	tmp/sample%_lineageF_pass.fastq.conformed tmp/sample%_lineageR_pass.fastq.conformed
//	-starcode \
//		-1 $(word 1,$^) -2 $(word 2,$^) \
//		-s \
//		-t 3 --print-clusters > $@
//
//tmp/sample%_clustered_lineage_tags.counts: \
//	tmp/sample%_clustered_lineage_tags.starclust
//	gawk '{print $3}' $< | sort -rg | uniq -c > $@
//
//
//process label_counts {
//    input:
//        set MM, file(this_counts_file) from results
//    output:
//        file "labeled_counts.tsv" into labeled_results
//    shell: 
//        '''
//        cat !{this_counts_file} | sed "s/^/!{MM}MM_/" > labeled_counts.tsv
//        '''
//}
//
//process concat_counts {
//    publishDir "tmp", mode: 'copy'
//    input:
//        file 'labeled_counts_list' from labeled_results.collect()
//    output:
//        file "all_counts.tsv"
//    shell: 
//        '''
//        head -n 1 labeled_counts_list1 | \
//            sed 's/^[0-9]*MM_//' | sed 's/Strain/Observation/' \
//            > all_counts.tsv
//        for i in $(ls labeled_counts_list*)
//            do cat ${i} | tail -n+2 >> all_counts.tsv
//        done
//        '''
//}

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


//process bartender_extract {
//    publishDir "tmp"
//    input:
//        set parameters, file(yaml), file(abundances), file(input_fastq), file(input_ranks) from grinder_output 
//    output:
//        set parameters, file(yaml), file(abundances), file(input_ranks), file("extracted_barcode.tsv") into bartender_extracted_codes 
//    shell:
//    '''
//    bartender_extractor_com -f !{input_fastq} -o extracted -q ? \
//        -p !{parameters["bartender_pattern"]} -m 2
//    '''
//}
//
//process bartender_quant {
//    publishDir "tmp"
//    input:
//        set parameters, file(yaml), file(abundances), file(input_ranks), file(extracted_barcodes) from bartender_extracted_codes 
//    output:
//        set parameters, file(yaml), file(abundances), file(input_ranks),
//            file("barcode_quant_barcode.csv"),
//            file("barcode_quant_cluster.csv"),
//            file("barcode_quant_quality.csv") into bartender_quantifications
//    shell:
//    '''
//    bartender_single_com -f !{extracted_barcodes} \
//        -o barcode_quant -d 3
//    '''
//}
//
//process bartender_eval {
//    publishDir "tmp"
//    input:
//        set parameters, file(yaml), file(abundances), file(input_ranks),
//            file(bartender_quant_barcode),
//            file(bartender_quant_cluster),
//            file(bartender_quant_quality) from bartender_quantifications
//        each script_eval_bartender_deviations
//    output:
//        set parameters, file("bartender_deviations.csv") into bartender_deviations
//    shell:
//    '''
//    python3 !{script_eval_bartender_deviations} !{yaml} !{abundances} !{bartender_quant_cluster} !{bartender_quant_cluster} "!{parameters['barcode_pattern']}"
//    '''
//}
