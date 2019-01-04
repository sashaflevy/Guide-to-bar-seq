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
//  4. replicate id
//  5. read depth
//  6. errors after sequencing
//  7. percentage of errors that are substitutions and indels
//  8. chimera percentage
//  9. abundance distribution
run_parameters = Channel.from([
    [ "ATCGNNNNCNNNNCNNNNCNNNNCNNNNATCG" , 100, 1000, "a", 10000, "uniform 0.1", "50 50", "10", "powerlaw 0.6" ],
    [ "ATCGNNNNCNNNNCNNNNCNNNNCNNNNATCG" , 100, 1000, "a", 10000, "uniform 0.1", "50 50", "10", "powerlaw 0.3" ]
    ])

process make_barcoded_libraries {
    input: 
        each script_make_degenerate_barcode_library_py
        each run_parameters
    output: 
        set run_parameters, 
            file("barcodeLibrary_${run_parameters[0]}_${run_parameters[1]}lineages_${run_parameters[2]}clones_replicate${run_parameters[3]}.yaml"), 
            file("barcodeLibrary_${run_parameters[0]}_${run_parameters[1]}lineages_${run_parameters[2]}clones_replicate${run_parameters[3]}.fasta") into barcoded_libraries
    shell:
    '''
    python3 !{script_make_degenerate_barcode_library_py} \
        --number-lineages !{run_parameters[1]} \
        --number-barcoded-clones !{run_parameters[2]} \
        --replicate !{run_parameters[3]} \
        !{run_parameters[0]} barcodeLibrary
    '''
}

process grinder {
    input:
        set run_parameters, file(yaml), file(fasta) from barcoded_libraries
    output:
        set run_parameters, file("grinder-reads.fastq"), 
            file("grinder-ranks.txt") into grinder_output
    shell:
    '''
    grinder -reference_file !{fasta} \
        -total_reads !{run_parameters[4]} \
        -read_dist 150 \
        -unidirectional 1 \
        -mutation_dist !{run_parameters[5]} \
        -mutation_ratio !{run_parameters[6]} \
        -chimera_perc !{run_parameters[7]} \
        -abundance_model !{run_parameters[8]} \
        -random_seed 1234 \
        -qual_levels 30 10 \
        -fastq_output 1 &> errorz
    '''
}

//Cli required arguments:
//    -rf <reference_file> | -reference_file <reference_file> | -gf
//    <reference_file> | -genome_file <reference_file>
//        FASTA file that contains the input reference sequences (full
//        genomes, 16S rRNA genes, transcripts, proteins...) or '-' to read
//        them from the standard input. See the README file for examples of
//        databases you can use and where to get them from. Default: -
//
//Cli optional arguments:
//    -tr <total_reads> | -total_reads <total_reads>
//        Number of shotgun or amplicon reads to generate for each library. Do
//        not specify this if you specify the fold coverage. Default: 100
//
//    -rd <read_dist>... | -read_dist <read_dist>...
//        Desired shotgun or amplicon read length distribution specified as:
//        average length, distribution ('uniform' or 'normal') and standard
//        deviation.
//
//        Only the first element is required. Examples:
//
//          All reads exactly 101 bp long (Illumina GA 2x): 101
//          Uniform read distribution around 100+-10 bp: 100 uniform 10
//          Reads normally distributed with an average of 800 and a standard deviation of 100
//            bp (Sanger reads): 800 normal 100
//          Reads normally distributed with an average of 450 and a standard deviation of 50
//            bp (454 GS-FLX Ti): 450 normal 50
//
//        Reference sequences smaller than the specified read length are not
//        used. Default: 100
//
//    -id <insert_dist>... | -insert_dist <insert_dist>...
//        Create paired-end or mate-pair reads spanning the given insert
//        length. Important: the insert is defined in the biological sense,
//        i.e. its length includes the length of both reads and of the stretch
//        of DNA between them: 0 : off, or: insert size distribution in bp, in
//        the same format as the read length distribution (a typical value is
//        2,500 bp for mate pairs) Two distinct reads are generated whether or
//        not the mate pair overlaps. Default: 0
//
//    -fr <forward_reverse> | -forward_reverse <forward_reverse>
//        Use DNA amplicon sequencing using a forward and reverse PCR primer
//        sequence provided in a FASTA file. The reference sequences and their
//        reverse complement will be searched for PCR primer matches. The
//        primer sequences should use the IUPAC convention for degenerate
//        residues and the reference sequences that that do not match the
//        specified primers are excluded. If your reference sequences are full
//        genomes, it is recommended to use <copy_bias> = 1 and <length_bias>
//        = 0 to generate amplicon reads. To sequence from the forward strand,
//        set <unidirectional> to 1 and put the forward primer first and
//        reverse primer second in the FASTA file. To sequence from the
//        reverse strand, invert the primers in the FASTA file and use
//        <unidirectional> = -1. The second primer sequence in the FASTA file
//        is always optional. Example: AAACTYAAAKGAATTGRCGG and
//        ACGGGCGGTGTGTRC for the 926F and 1392R primers that target the V6 to
//        V9 region of the 16S rRNA gene.
//
//    -un <unidirectional> | -unidirectional <unidirectional>
//        Instead of producing reads bidirectionally, from the reference
//        strand and its reverse complement, proceed unidirectionally, from
//        one strand only (forward or reverse). Values: 0 (off, i.e.
//        bidirectional), 1 (forward), -1 (reverse). Use <unidirectional> = 1
//        for amplicon and strand-specific transcriptomic or proteomic
//        datasets. Default: 0
//
//    -cb <copy_bias> | -copy_bias <copy_bias>
//        In amplicon libraries where full genomes are used as input, sample
//        species proportionally to the number of copies of the target gene:
//        at equal relative abundance, genomes that have multiple copies of
//        the target gene contribute more amplicon reads than genomes that
//        have a single copy. 0 = no, 1 = yes. Default: 1
//
//    -md <mutation_dist>... | -mutation_dist <mutation_dist>...
//        Introduce sequencing errors in the reads, under the form of
//        mutations (substitutions, insertions and deletions) at positions
//        that follow a specified distribution (with replacement): model
//        (uniform, linear, poly4), model parameters. For example, for a
//        uniform 0.1% error rate, use: uniform 0.1. To simulate Sanger
//        errors, use a linear model where the errror rate is 1% at the 5' end
//        of reads and 2% at the 3' end: linear 1 2. To model Illumina errors
//        using the 4th degree polynome 3e-3 + 3.3e-8 * i^4 (Korbel et al
//        2009), use: poly4 3e-3 3.3e-8. Use the <mutation_ratio> option to
//        alter how many of these mutations are substitutions or indels.
//        Default: uniform 0 0
//
//    -mr <mutation_ratio>... | -mutation_ratio <mutation_ratio>...
//        Indicate the percentage of substitutions and the number of indels
//        (insertions and deletions). For example, use '80 20' (4
//        substitutions for each indel) for Sanger reads. Note that this
//        parameter has no effect unless you specify the <mutation_dist>
//        option. Default: 80 20
//
//    -cp <chimera_perc> | -chimera_perc <chimera_perc>
//        Specify the percent of reads in amplicon libraries that should be
//        chimeric sequences. The 'reference' field in the description of
//        chimeric reads will contain the ID of all the reference sequences
//        forming the chimeric template. A typical value is 10% for amplicons.
//        This option can be used to generate chimeric shotgun reads as well.
//        Default: 0 %
//
//    -cd <chimera_dist>... | -chimera_dist <chimera_dist>...
//        Specify the distribution of chimeras: bimeras, trimeras, quadrameras
//        and multimeras of higher order. The default is the average values
//        from Quince et al. 2011: '314 38 1', which corresponds to 89% of
//        bimeras, 11% of trimeras and 0.3% of quadrameras. Note that this
//        option only takes effect when you request the generation of chimeras
//        with the <chimera_perc> option. Default: 314 38 1
//
//    -ck <chimera_kmer> | -chimera_kmer <chimera_kmer>
//        Activate a method to form chimeras by picking breakpoints at places
//        where k-mers are shared between sequences. <chimera_kmer> represents
//        k, the length of the k-mers (in bp). The longer the kmer, the more
//        similar the sequences have to be to be eligible to form chimeras.
//        The more frequent a k-mer is in the pool of reference sequences
//        (taking into account their relative abundance), the more often this
//        k-mer will be chosen. For example, CHSIM (Edgar et al. 2011) uses
//        this method with a k-mer length of 10 bp. If you do not want to use
//        k-mer information to form chimeras, use 0, which will result in the
//        reference sequences and breakpoints to be taken randomly on the
//        "aligned" reference sequences. Note that this option only takes
//        effect when you request the generation of chimeras with the
//        <chimera_perc> option. Also, this options is quite memory intensive,
//        so you should probably limit yourself to a relatively small number
//        of reference sequences if you want to use it. Default: 10 bp
//
//    -am <abundance_model>... | -abundance_model <abundance_model>...
//        Relative abundance model for the input reference sequences: uniform,
//        linear, powerlaw, logarithmic or exponential. The uniform and linear
//        models do not require a parameter, but the other models take a
//        parameter in the range [0, infinity). If this parameter is not
//        specified, then it is randomly chosen. Examples:
//
//          uniform distribution: uniform
//          powerlaw distribution with parameter 0.1: powerlaw 0.1
//          exponential distribution with automatically chosen parameter: exponential
//
//        Default: uniform 1
//
//    -mi <multiplex_ids> | -multiplex_ids <multiplex_ids>
//        Specify an optional FASTA file that contains multiplex sequence
//        identifiers (a.k.a MIDs or barcodes) to add to the sequences (one
//        sequence per library, in the order given). The MIDs are included in
//        the length specified with the -read_dist option and can be altered
//        by sequencing errors. See the MIDesigner or BarCrawl programs to
//        generate MID sequences.
//
//    -rs <random_seed> | -random_seed <random_seed>
//        Seed number to use for the pseudo-random number generator.
//
//    -ql <qual_levels>... | -qual_levels <qual_levels>...
//        Generate basic quality scores for the simulated reads. Good residues
//        are given a specified good score (e.g. 30) and residues that are the
//        result of an insertion or substitution are given a specified bad
//        score (e.g. 10). Specify first the good score and then the bad score
//        on the command-line, e.g.: 30 10. Default:
//
//    -fq <fastq_output> | -fastq_output <fastq_output>
//        Whether to write the generated reads in FASTQ format (with
//        Sanger-encoded quality scores) instead of FASTA and QUAL or not (1:
//        yes, 0: no). <qual_levels> need to be specified for this option to
//        be effective. Default: 0
//
//    -bn <base_name> | -base_name <base_name>
//        Prefix of the output files. Default: grinder

//process slapchop_or_not {
//    input:
//        set run_parameters, file("grinder-reads.fastq"), 
//            file("grinder-ranks.txt") into grinder_output
//    output:
//    shell:
//    '''
//    grinder -reference_file !{fasta} \
//        -total_reads !{run_parameters[4]} \
//        -read_dist 150 \
//        -unidirectional 1 \
//        -mutation_dist !{run_parameters[5]} \
//        -mutation_ratio !{run_parameters[6]} \
//        -chimera_perc !{run_parameters[7]} \
//        -abundance_model !{run_parameters[8]} \
//        -random_seed 1234 \
//        -qual_levels 30 10 \
//        -fastq_output 1 &> errorz
//    '''
//}
