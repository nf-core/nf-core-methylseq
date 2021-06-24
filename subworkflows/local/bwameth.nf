/*
 * bwameth subworkflow
 */

def modules = params.modules.clone()

def bwameth_index_options   = modules['bwameth_index']
if (!params.save_reference) { bwameth_index_options['publish_files'] = false }
def samtools_faidx_options   = modules['samtools_faidx']
if (!params.save_reference) { samtools_faidx_options['publish_files'] = false }

def bwameth_align_options   = modules['bwameth_align']
if (!params.save_align_intermeds)  { bwameth_align_options['publish_files'] = false }
def samtools_sort_options   = modules['samtools_sort'].clone()
if (!params.save_align_intermeds)  { samtools_sort_options['publish_files'] = false }
samtools_sort_options['suffix']   = ".sorted"


def samtools_index_alignments_options   = modules['samtools_index'].clone()
if (!params.save_align_intermeds)  { samtools_index_alignments_options['publish_files'] = false }
samtools_index_alignments_options['publish_dir']   = "${params.aligner}/alignments"

def samtools_index_deduplicated_options   = modules['samtools_index'].clone()
samtools_index_deduplicated_options['publish_dir']   = "${params.aligner}/deduplicated"

def methyldackel_extract_options   = modules['methyldackel_extract']
methyldackel_extract_options.args += params.comprehensive ? ' --CHG --CHH' : ''
methyldackel_extract_options.args += params.ignore_flags ? " --ignoreFlags" : ''
methyldackel_extract_options.args += params.methyl_kit ? " --methylKit" : ''
methyldackel_extract_options.args += params.min_depth > 0 ? " --minDepth ${params.min_depth}" : ''

def methyldackel_mbias_options   = modules['methyldackel_mbias']
methyldackel_mbias_options.args += params.comprehensive ? ' --CHG --CHH' : ''
methyldackel_mbias_options.args += params.ignore_flags ? " --ignoreFlags" : ''

include { SAMTOOLS_STATS                                } from '../../modules/nf-core/software/samtools/stats/main'        addParams( options: modules['samtools_stats']           )
include { SAMTOOLS_FLAGSTAT                             } from '../../modules/nf-core/software/samtools/flagstat/main'     addParams( options: modules['samtools_flagstat']        )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../modules/nf-core/software/samtools/index/main'        addParams( options: samtools_index_alignments_options   )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../modules/nf-core/software/samtools/index/main'        addParams( options: samtools_index_deduplicated_options )
include { PICARD_MARKDUPLICATES                         } from '../../modules/nf-core/software/picard/markduplicates/main' addParams( options: modules['picard_markduplicates']    )
include { SAMTOOLS_SORT                                 } from '../../modules/nf-core/software/samtools/sort/main'         addParams( options: samtools_sort_options               )
include { SAMTOOLS_FAIDX                                } from '../../modules/local/software/samtools/faidx/main'          addParams( options: samtools_faidx_options              )
include { BWAMETH_INDEX                                 } from '../../modules/local/software/bwameth/index/main'           addParams( options: bwameth_index_options               )
include { BWAMETH_ALIGN                                 } from '../../modules/local/software/bwameth/align/main'           addParams( options: bwameth_align_options               )

include { METHYLDACKEL } from './methyldackel' addParams( mbias_options: methyldackel_mbias_options, extract_options: methyldackel_extract_options )

workflow BWAMETH {
    take:
    genome // channel: [ val(meta), [ genome ] ]
    reads  // channel: [ val(meta), [ reads ] ]

    main:
    /*
     * Generate bwameth index if not supplied
     */

    // there might be indices from iGenomes or user input
    // branch them into different channels
    genome
        .branch{ meta, genome ->
            have_bwameth_index: genome.containsKey('bwameth_index')
                return [meta, genome.bwameth_index]
            need_bwameth_index: !genome.containsKey('bwameth_index')
                return [meta, genome.fasta] 
        }
        .set{ch_genome}

    // group by unique fastas, so that we only index each genome once
    ch_genome.need_bwameth_index.groupTuple(by:1) | BWAMETH_INDEX

    // reverse groupTuple to restore the cardinality of the input channel
    // then join back with pre-existing indices
    bwameth_index = BWAMETH_INDEX.out.index | transpose(by:0) | mix(ch_genome.have_bwameth_index)

    /*
     * Align with bwameth
     */
    BWAMETH_ALIGN (
        reads.join(bwameth_index)
    )

    /*
     * Run samtools sort
     */
    SAMTOOLS_SORT (BWAMETH_ALIGN.out.bam)

    /*
     * Run samtools index on alignment
     */
    SAMTOOLS_INDEX_ALIGNMENTS (SAMTOOLS_SORT.out.bam)

    /*
     * Run samtools flagstat and samtools stats
     */
    BWAMETH_ALIGN.out.bam.join(SAMTOOLS_INDEX_ALIGNMENTS.out.bai) | (SAMTOOLS_FLAGSTAT & SAMTOOLS_STATS)

    if (params.skip_deduplication || params.rrbs) {
        alignments = SAMTOOLS_SORT.out.bam
        bam_index = SAMTOOLS_INDEX_ALIGNMENTS.out.bai
        picard_metrics = Channel.empty()
        picard_version = Channel.empty()
    } else {
        /*
        * Run Picard MarkDuplicates
        */
        PICARD_MARKDUPLICATES (SAMTOOLS_SORT.out.bam)   
        /*
         * Run samtools index on deduplicated alignment
        */
        SAMTOOLS_INDEX_DEDUPLICATED (PICARD_MARKDUPLICATES.out.bam)

        alignments = PICARD_MARKDUPLICATES.out.bam
        bam_index = SAMTOOLS_INDEX_DEDUPLICATED.out.bai
        picard_metrics = PICARD_MARKDUPLICATES.out.metrics
        picard_version = PICARD_MARKDUPLICATES.out.version
    }

    /*
     * Generate fasta index if not supplied
     */
    genome
        .branch{ meta, genome ->
            have_fasta_index: genome.containsKey('fasta_index')
                return [meta, genome.fasta_index]
            need_fasta_index: !genome.containsKey('fasta_index')
                return [meta, genome.fasta] 
        }
        .set{ch_genome}

    ch_genome.need_fasta_index.groupTuple(by:1) | SAMTOOLS_FAIDX

    fasta_index = SAMTOOLS_FAIDX.out.fai | transpose | mix(ch_genome.have_fasta_index)

    /*
     * SUBWORKFLOW: Extract per-base methylation and plot methylation bias
     */
    METHYLDACKEL (
        alignments,
        bam_index,
        genome.map{ meta, genome -> [meta, genome.fasta] },
        fasta_index
    )

    if (!params.skip_multiqc) {
        /*
        * Collect MultiQC inputs
        */
        picard_metrics.collect{ it[1] }
            .mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{ it[1] })
            .mix(SAMTOOLS_STATS.out.stats.collect{ it[1] })
            .mix(METHYLDACKEL.out.mbias.collect{ it[1] })
            .mix(METHYLDACKEL.out.consensus.collect{ it[1] })
            .set{ ch_multiqc_files }
    } else {
        ch_multiqc_files = Channel.empty()
    }

    /*
     * Collect Software Versions
     */
    picard_version
        .mix(METHYLDACKEL.out.version)
        .mix(SAMTOOLS_SORT.out.version)
        .mix(BWAMETH_INDEX.out.version)
        .set{ versions }

    emit:
    bam                  = SAMTOOLS_SORT.out.bam          // channel: [ val(meta), [ bam ] ]
    dedup                = alignments                     // channel: [ val(meta), [ bam ] ]

    mqc                  = ch_multiqc_files               // path: *{html,txt}
    versions                                              // path: *.version.txt
}