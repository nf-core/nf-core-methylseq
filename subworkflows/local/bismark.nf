/*
 * bismark subworkflow
 */

def modules = params.modules.clone()

def bismark_genomepreparation_options   = modules['bismark_genomepreparation']
bismark_genomepreparation_options.args += params.aligner == 'bismark_hisat' ? ' --hisat2' : ' --bowtie2'
bismark_genomepreparation_options.args += params.slamseq ? ' --slam' : ''
if (!params.save_reference) { bismark_genomepreparation_options['publish_files'] = false }

def bismark_align_options   = modules['bismark_align']
known_splices = params.known_splices ? file("${params.known_splices}", checkIfExists: true) : file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

bismark_align_options.args += params.aligner == 'bismark_hisat' ? ' --hisat2' : ' --bowtie2'
bismark_align_options.args += params.aligner == "bismark_hisat" && known_splices.name != 'dummy_file.txt' ? " --known-splicesite-infile <(hisat2_extract_splice_sites.py ${known_splices})" : ''
bismark_align_options.args += params.pbat ? ' --pbat' : ''
bismark_align_options.args += ( params.single_cell || params.non_directional || params.zymo ) ? ' --non_directional' : ''
bismark_align_options.args += params.unmapped ? ' --unmapped' : ''
bismark_align_options.args += params.relax_mismatches ? " --score_min L,0,-${params.num_mismatches}" : ''
bismark_align_options.args += params.local_alignment ? " --local" : ''
bismark_align_options.args += params.minins ? " --minins ${params.minins}" : ''
bismark_align_options.args += params.maxins ? " --minins ${params.maxins}" : ''
if (params.save_align_intermeds)  { bismark_align_options.publish_files.put('bam','') }

def samtools_sort_options = modules['samtools_sort'].clone()
samtools_sort_options['suffix']   = ".deduplicated.sorted"

def bismark_methylationextractor_options   = modules['bismark_methylationextractor']
bismark_methylationextractor_options.args += params.comprehensive   ? ' --comprehensive --merge_non_CpG' : ''
bismark_methylationextractor_options.args += params.cytosine_report ? ' --cytosine_report --genome_folder BismarkIndex' : ''
bismark_methylationextractor_options.args += params.meth_cutoff     ? " --cutoff ${params.meth_cutoff}" : ''

include { BISMARK_GENOMEPREPARATION    } from '../../modules/local/software/bismark/genomepreparation/main'    addParams( options: bismark_genomepreparation_options    )
include { BISMARK_ALIGN                } from '../../modules/local/software/bismark/align/main'                addParams( options: bismark_align_options                )
include { BISMARK_METHYLATIONEXTRACTOR } from '../../modules/local/software/bismark/methylationextractor/main' addParams( options: bismark_methylationextractor_options )
include { SAMTOOLS_SORT                } from '../../modules/nf-core/software/samtools/sort/main'              addParams( options: samtools_sort_options                )
include { BISMARK_DEDUPLICATE          } from '../../modules/nf-core/software/bismark/deduplicate/main'        addParams( options: modules['bismark_deduplicate']       )
include { BISMARK_REPORT               } from '../../modules/nf-core/software/bismark/report/main'             addParams( options: modules['bismark_report']            )
include { BISMARK_SUMMARY              } from '../../modules/nf-core/software/bismark/summary/main'            addParams( options: modules['bismark_summary']           )

workflow BISMARK {
    take:
    genome // channel: [ val(meta), [ genome ] ]
    reads  // channel: [ val(meta), [ reads ] ]

    main:
    /*
     * Generate bismark index if not supplied
     */

    // there might be indices from user input
    // branch them into different channels
    genome
        .branch{ meta, genome ->
            have_index: genome.containsKey('bismark_index')
                return [meta, genome.bismark_index]
            need_index: !genome.containsKey('bismark_index')
                return [meta, genome.fasta] 
        }
        .set{ch_genome}

    // group by unique fastas, so that we only index each genome once
    ch_genome.need_index.groupTuple(by:1) | BISMARK_GENOMEPREPARATION

    // reverse groupTuple to restore the cardinality of the input channel
    // then join back with pre-existing indices
    bismark_index = BISMARK_GENOMEPREPARATION.out.index | transpose | mix(ch_genome.have_index)

    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        reads.join(bismark_index)
    )

    if (params.skip_deduplication || params.rrbs) {
        alignments = BISMARK_ALIGN.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    } else {
        /*
        * Run deduplicate_bismark
        */
        BISMARK_DEDUPLICATE( BISMARK_ALIGN.out.bam )

        alignments = BISMARK_DEDUPLICATE.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
    }

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments.join(bismark_index)
    )

    /*
     * Generate bismark sample reports
     */
    BISMARK_REPORT (
    alignment_reports
        .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
        .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )

    /*
     * Generate bismark summary report
     */
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[2] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] }.ifEmpty([])
    )

    /*
     * MODULE: Run samtools sort
     */
    SAMTOOLS_SORT (
        alignments
    )

    if (!params.skip_multiqc) {
        /*
        * Collect MultiQC inputs
        */
        BISMARK_SUMMARY.out.summary.ifEmpty([])
            .mix(alignment_reports.collect{ it[1] })
            .mix(alignment_reports.collect{ it[2] })
            .mix(BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] })
            .mix(BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] })
            .mix(BISMARK_REPORT.out.report.collect{ it[1] })
            .set{ ch_multiqc_files }
    } else {
        ch_multiqc_files = Channel.empty()
    }

    /*
     * Collect Software Versions
     */
    SAMTOOLS_SORT.out.version
        .mix(BISMARK_ALIGN.out.version)
        .set{ versions }
    

    emit:
    bam              = BISMARK_ALIGN.out.bam          // channel: [ val(meta), [ bam ] ]
    dedup            = SAMTOOLS_SORT.out.bam          // channel: [ val(meta), [ bam ] ]

    mqc              = ch_multiqc_files               // path: *{html,txt}
    versions                                          // path: *.version.txt

}