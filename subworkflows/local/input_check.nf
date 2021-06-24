/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .set { sample }
        
    reads = sample.map { get_samplesheet_paths(it) }
    genome = sample.map { get_genome_paths(it, params.genomes) }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
    genome // channel: [ val(meta), [ fasta ] ]
}

// Function to get a HashMap with id and single_end fields
def get_metamap(LinkedHashMap sample) {
    def meta = [:]
    meta.id           = sample.sample
    meta.single_end   = sample.single_end.toBoolean()

    return meta
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap sample) {
    def meta = get_metamap(sample)

    if (!file(sample.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${sample.fastq_1}"
    }
    if (meta.single_end) {
        return [ meta, [ file(sample.fastq_1) ] ]
    } else {
        if (!file(sample.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${sample.fastq_2}"
        }
        return [ meta, [ file(sample.fastq_1), file(sample.fastq_2) ] ]
    }
}

// Function to get list of [ meta, genome ]
def get_genome_paths(LinkedHashMap sample, LinkedHashMap genomeMap) {
    def meta = get_metamap(sample)

    def genome = [:]

    if (sample.genome) {
        if (genomeMap && genomeMap.containsKey(sample.genome)) {
            // get fasta and indices from iGenomes
            genome.fasta = file(genomeMap[sample.genome].fasta, checkIfExists: true)
            if (params.aligner == 'bismark') {
                genome.bismark_index = file(genomeMap[sample.genome].bismark, checkIfExists: true)
            }
            else if (params.aligner == 'bwameth') {
                genome.bwa_meth = file(genomeMap[sample.genome].bwa_meth, checkIfExists: true)
                genome.fasta_index = file(genomeMap[sample.genome].fasta_index, checkIfExists: true)
            }
        } else {
            // genome is a fasta file, or not part of iGenomes
            genome.fasta = file(sample.genome, checkIfExists: true)
        }
    } else if ( params.fasta ) {
      // samplesheet does not contain genome column, fall back to params.fasta
      genome.fasta = file(params.fasta, checkIfExists: true)
    } else {
        exit 1, "ERROR: Please either supply a fasta file with --fasta or specify genome column in the samplesheet"
    }

    return [ meta, genome ]
}