////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.known_splices ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Trimming / kit presets
if(params.pbat){
    params.clip_r1 = 9
    params.clip_r2 = 9
    params.three_prime_clip_r1 = 9
    params.three_prime_clip_r2 = 9
}
else if( params.single_cell ){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 6
}
else if( params.epignome ){
    params.clip_r1 = 8
    params.clip_r2 = 8
    params.three_prime_clip_r1 = 8
    params.three_prime_clip_r2 = 8
}
else if( params.accel || params.zymo ){
    params.clip_r1 = 10
    params.clip_r2 = 15
    params.three_prime_clip_r1 = 10
    params.three_prime_clip_r2 = 10
}
else if( params.cegx ){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 2
    params.three_prime_clip_r2 = 2
}
else if( params.em_seq ){
    params.maxins = 1000
    params.clip_r1 = 8
    params.clip_r2 = 8
    params.three_prime_clip_r1 = 8
    params.three_prime_clip_r2 = 8
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

def trimgalore_options   = modules['trimgalore']
trimgalore_options.args += params.rrbs ? ' --rrbs' : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }


// Local: Modules
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['csv':'']] )

// Local: Sub-workflows
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

// nf-core/modules: Modules
include { FASTQC          } from '../modules/nf-core/software/fastqc/main'          addParams( options: modules['fastqc']          )
include { TRIMGALORE      } from '../modules/nf-core/software/trimgalore/main'      addParams( options: trimgalore_options         )
include { QUALIMAP_BAMQC  } from '../modules/nf-core/software/qualimap/bamqc/main'  addParams( options: modules['qualimap_bamqc']  )
include { PRESEQ_LCEXTRAP } from '../modules/nf-core/software/preseq/lcextrap/main' addParams( options: modules['preseq_lcextrap'] )
include { MULTIQC         } from '../modules/nf-core/software/multiqc/main'         addParams( options: multiqc_options            )

// nf-core/modules: Sub-workflows
if( params.aligner =~ /bismark/ ){
    include { BISMARK as ALIGNER } from '../subworkflows/local/bismark'
} else if ( params.aligner == 'bwameth' ){
    include { BWAMETH as ALIGNER } from '../subworkflows/local/bwameth'
}

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow METHYLSEQ {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( 
        ch_input
    )
    
    /*
     * MODULE: Run FastQC
     */
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    
    if (!params.skip_trimming) {
        /*
        * MODULE: Run TrimGalore!
        */
        TRIMGALORE(INPUT_CHECK.out.reads)

        reads = TRIMGALORE.out.reads
        ch_software_versions = ch_software_versions.mix(TRIMGALORE.out.version.first().ifEmpty(null))
    } else {
        reads = INPUT_CHECK.out.reads
    }

    /*
     * SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark 
     */
    ALIGNER (
        INPUT_CHECK.out.genome,
        reads
    )
    ch_software_versions = ch_software_versions.mix(ALIGNER.out.versions.unique{ it.baseName }.ifEmpty(null))

    /*
     * MODULE: Qualimap BamQC
     */
    QUALIMAP_BAMQC (
        ALIGNER.out.dedup,
        ch_dummy_file,
        false
    )
    ch_software_versions = ch_software_versions.mix(QUALIMAP_BAMQC.out.version.first().ifEmpty(null))

    /*
     * MODULE: Run Preseq
     */
    PRESEQ_LCEXTRAP (
        ALIGNER.out.bam
    )
    ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.version.first().ifEmpty(null))

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    /*
     * MultiQC
     */  
    if (!params.skip_multiqc) {
        workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        Channel.from(ch_multiqc_config)
            .mix(ch_multiqc_custom_config.collect())
            .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            .mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
            .mix(FASTQC.out.zip.collect{ it[1] })
            .mix(QUALIMAP_BAMQC.out.results.collect{ it[1] })
            .mix(PRESEQ_LCEXTRAP.out.ccurve.collect{ it[1] })
            .mix(ALIGNER.out.mqc)
            .set { ch_multiqc_files }

        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{ it[1] })
        }

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////