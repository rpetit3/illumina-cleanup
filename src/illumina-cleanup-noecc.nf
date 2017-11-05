#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
params.help = null
params.output = null
params.fq1 = null
params.fq2 = null
params.cpu = 1
params.genome_size = 2814816 // S. aureus, adjust accordingly
params.reference_path = '/opt/references'
params.coverage = 100
params.is_miseq = false
params.force = false
params.clear_cache_on_success = true
params.clear_logs = true

if (params.help) {
    print_usage()
    exit 0
}

check_input_params()

// Set some global variables
sample = params.sample
outdir = params.output ? params.output : './'
is_miseq = params.is_miseq
genome_size = params.genome_size
coverage = params.coverage
reference_path = params.reference_path
cpu = params.cpu

// Output folders
logs_folder= outdir + "/logs"

/* ==== BEGIN FASTQ CLEANUP ==== */
process original_fastq_stats {
    publishDir logs_folder, mode: 'copy', overwrite: true

    input:
        file fq from create_input_channel(params.fq1, params.fq2)
    output:
        file "original.json" into ORIGINAL_JSON
    shell:
        '''
        zcat !{fq[0]} !{fq[1]} | fastq-stats > "original.json" !{genome_size}
        '''
}

process bbduk_adapter_filter {
    publishDir logs_folder, mode: 'copy', overwrite: true, pattern: "*.log"

    input:
        file fq from create_input_channel(params.fq1, params.fq2)
    output:
        file 'bbduk-adapter-R*.fq' into BBDUK, BBDUK_STATS
        file '*.log'
    shell:
        '''
        bbduk.sh -Xmx2g threads=!{cpu} in=!{fq[0]} in2=!{fq[1]} out=bbduk-phix-R1.fq \
        out2=bbduk-phix-R2.fq stats=bbduk-phix.log hdist=1 k=31 overwrite=t \
        ordered=t ref=!{reference_path}/phiX-NC_001422.fasta

        bbduk.sh -Xmx2g threads=!{cpu} in=bbduk-phix-R1.fq in2=bbduk-phix-R2.fq \
        out=bbduk-adapter-R1.fq out2=bbduk-adapter-R2.fq stats=bbduk-adapter.txt \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo qout=33 minlength=36  overwrite=t \
        ordered=t ref=!{reference_path}/adapters.fasta
        cp .command.err bbduk-stderr.log
        cp .command.out bbduk-stdout.log
        '''
}

process bbduk_fastq_stats {
    publishDir logs_folder, mode: 'copy', overwrite: true

    input:
        file fq from BBDUK_STATS
    output:
        file "post-adapter.json" into BBDUK_FQ_STATS, ADAPTER_JSON
    shell:
        '''
        cat !{fq[0]} !{fq[1]} | fastq-stats > post-adapter.json !{genome_size}
        '''
}

process illumina_cleanup {
    input:
        file fq from BBDUK
        file stats from BBDUK_FQ_STATS
    output:
        file '*cleanup.fastq' into FASTQ_STATS, FASTQ_SPLIT
    shell:
        no_length_filter = is_miseq ? '--no_length_filter' : ''
        '''
        fastq-interleave !{fq[0]} !{fq[1]} | illumina-cleanup.py --paired --stats !{stats} \
        --coverage !{params.coverage} --genome_size !{genome_size} !{no_length_filter} > \
        !{sample}.cleanup.fastq
        '''
}

process final_stats {
    publishDir logs_folder, mode: 'copy', overwrite: true

    input:
        file fq from FASTQ_STATS
    output:
        file "cleanup.json" into CLEANUP_JSON
    shell:
        '''
        cat !{fq} | fastq-stats > cleanup.json !{genome_size}
        '''
}

process split_fastq {
    publishDir outdir, mode: 'copy', overwrite: true

    input:
        file fq from FASTQ_STATS
    output:
        file "*.fastq.gz"
    shell:
        '''
        reformat.sh in=!{fq} out1=!{sample}-R1.cleanup.fastq out2=!{sample}-R2.cleanup.fastq
        gzip --best !{sample}-R1.cleanup.fastq
        gzip --best !{sample}-R2.cleanup.fastq
        '''
}

process merge_json {
    publishDir outdir, mode: 'copy', overwrite: true

    input:
        file original from ORIGINAL_JSON
        file adapter from ADAPTER_JSON
        file spades from SPADES_JSON
        file cleanup from CLEANUP_JSON
    output:
        file {"${sample}-illumina-cleanup.json"}
    shell:
        '''
        merge-json.py !{original} !{adapter} !{spades} !{cleanup} > !{sample}-illumina-cleanup.json
        '''
}

/* ==== END FASTQ CLEANUP ==== */

workflow.onComplete {
    if (workflow.success == true && params.clear_cache_on_success) {
        // No need to resume completed run so remove cache.
        file('./work/').deleteDir()
        if (params.clear_logs) {
            file(logs_folder).deleteDir()
        }
    }
    println """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

// Utility Functions
def print_usage() {
    log.info 'Illumina Cleanup Pipeline'
    log.info ''
    log.info 'Required Options:'
    log.info '    --fq1  FASTQ.GZ    Input FASTQ, compressed using GZIP'
    log.info '    --fq2  FASTQ.GZ    Second set of reads for paired end input.'
    log.info '    --sample  STR      A sample name to give the run.'
    log.info ''
    log.info 'Optional:'
    log.info '    --outdir  DIR      Directory to write results to. (Default ./${NAME})'
    log.info '    --genome_size  INT Expected genome size (bp) for coverage estimation.'
    log.info '    --coverage  INT    Reduce samples to a given coverage. (Default: 100x)'
    log.info '    --miseq            For Illumina MiSeq (variable read lengths), reads '
    log.info '                           will not be filtered base on read lengths.'
    log.info '    --help          Show this message and exit'
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow fastq-cleanup.nf --fq1 input.fastq.gz --fq2 input.fastq.gz --sample saureus [more options]'
}

def check_input_params() {
    error = false
    if (!params.sample) {
        log.info('A sample name is required to continue. Please use --sample')
        error = true
    }
    if (!params.fq1) {
        log.info('Compressed FASTQ (gzip) is required. Please use --fq1')
        error = true
    } else if (!file(params.fq1).exists()) {
        log.info('Invailid input (--fq1), please verify "' + params.fq1 + '"" exists.')
        error = true
    }

    if (!params.fq2) {
        log.info('Compressed FASTQ (gzip) is required. Please use --fq2')
        error = true
    } else if (!file(params.fq2).exists()) {
        log.info('Invailid input (--fq2), please verify "' + params.fq2 + '"" exists.')
        error = true
    }

    if (error) {
        log.info('See --help for more information')
        exit 1
    }
}

def create_input_channel(input_1, input_2) {
    return Channel.value([file(input_1), file(input_2)])
}
