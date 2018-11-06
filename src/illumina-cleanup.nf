#!/usr/bin/env nextflow
NAME = 'illumina-cleanup'
VERSION = 0.4
params.help = null
params.version = null
params.example_fastqs = null
params.cpus = 1
params.fastqs = null
params.outdir = null
params.single_end = null
params.genome_size = null
params.coverage = '100'
params.adapters = 'adapters'
params.ktrim = 'r'
params.adapter_k = 23
params.mink = 11
params.hdist = 1
params.tpe = 't'
params.tbo = 't'
params.phix = 'phix'
params.phix_k = 31
params.qtrim = 'rl'
params.trimq = 6
params.minlength = 35
params.maq = 20
params.qout = 33
params.tossjunk = 't'
params.ftm = 5
params.maxcor = 1
params.mash_s = 10000
params.mash_k = 31
params.mash_m = 3

// Parse input parameters
if (workflow.commandLine.endsWith(workflow.scriptName)) print_usage();
if (params.example_fastqs) print_example_fastqs();
if (params.version) print_version();
if (params.help) print_usage();
check_input_params()
check_input_fastqs(params.fastqs)

// Setup output directories
outdir = params.outdir ? params.outdir : './'
log_folder = outdir + '/logs'
summary_folder = outdir + '/summary'

process estimate_genome_size {
    /* Estimate a genome size if one is not given at runtime. */
    tag "${sample}"
    cpus 1
    publishDir "${outdir}/${sample}/sequence-qc/logs", mode: 'rellink', overwrite: true, pattern: "*.log"

    input:
    set val(sample), val(single_end), file(fq) from create_fastq_channel(params.fastqs)

    output:
    file 'mash-sketch-original.log'
    stdout into MASH_OUTPUT

    shell:
    if (params.genome_size)
        """
        echo !{params.genome_size}
        """
    else
        """
        mash sketch -o temp -s 10000 -k 31 -m 3 -r !{fq[0]} 2>&1 1> mash-sketch-original.log
        mash info -t temp.msh | tail -n 1 | cut -f2,2
        """
}

// Pull out genome size from Mash output
MASH_OUTPUT
    .flatMap{ it.tokenize('\n') }
    .last()
    .into { GS_ORIGINAL_SUMMARY; GS_ERROR_CORRECTION; GS_REDUCE_COVERAGE; GS_FINAL_SUMMARY }


process original_summary {
    /* Run FASTQC on the input FASTQ files. */
    cpus params.cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/sequence-qc/summary", mode: 'rellink', overwrite: true, pattern: '*.{html,json,zip}'

    input:
    set val(sample), val(single_end), file(fq) from create_fastq_channel(params.fastqs)
    val genome_size from GS_ORIGINAL_SUMMARY

    output:
    file '*.json'
    file '*fastqc.html'
    file '*fastqc.zip'

    shell:
    fq2 = single_end ? '' : fq[1]
    """
    zcat !{fq[0]} !{fq2} > !{sample}-original
    cat !{sample}-original | fastq-scan !{genome_size} > fastq-scan-original.json
    fastqc --noextract -f fastq -t !{params.cpus} !{sample}-original
    """
}

process adapter_removal {
    /* Remove Illumina related adapters using BBDuk */
    cpus params.cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/sequence-qc/logs", mode: 'rellink', overwrite: true, pattern: "*.log"

    input:
    set val(sample), val(single_end), file(fq) from create_fastq_channel(params.fastqs)

    output:
    set val(sample), val(single_end), file('adapter-r*.fq') into PHIX_REMOVAL
    file '*.log'

    shell:
    in2 = single_end ? '' : 'in2=' + fq[1]
    out2 = single_end ? '' : 'out2=adapter-r2.fq'
    """
    bbduk.sh \
        in=!{fq[0]} !{in2} \
        out=adapter-r1.fq !{out2} \
        ref=!{params.adapters} \
        k=!{params.adapter_k} \
        ktrim=!{params.ktrim} \
        mink=!{params.mink} \
        hdist=!{params.hdist} \
        tpe=!{params.tpe} \
        tbo=!{params.tbo} \
        threads=!{params.cpus} \
        ftm=!{params.ftm} \
        ordered=t \
        stats=bbduk-adapter.log
    """
}

process phix_removal {
    /* Remove contaminant (phiX) and quality-trim using BBDuk */
    cpus params.cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/sequence-qc/logs", mode: 'rellink', overwrite: true, pattern: "*.log"

    input:
    set val(sample), val(single_end), file(fq) from PHIX_REMOVAL

    output:
    set val(sample), val(single_end), file('phix-r*.fq') into ERROR_CORRECTION
    file '*.log'

    shell:
    in2 = single_end ? '' : 'in2=' + fq[1]
    out2 = single_end ? '' : 'out2=phix-r2.fq'
    """
    bbduk.sh \
        in=!{fq[0]} !{in2} \
        out=phix-r1.fq !{out2} \
        ref=!{params.phix} \
        k=!{params.phix_k} \
        hdist=!{params.hdist} \
        tpe=!{params.tpe} \
        tbo=!{params.tbo} \
        qtrim=!{params.qtrim} \
        trimq=!{params.trimq} \
        minlength=!{params.minlength} \
        minavgquality=!{params.maq} \
        qout=!{params.qout} \
        tossjunk=!{params.tossjunk} \
        threads=!{params.cpus} \
        ordered=t \
        stats=bbduk-phix.log
    """
}

process error_correction {
    /* Attempt to correct any base-call errors using Lighter */
    cpus params.cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/sequence-qc/logs", mode: 'rellink', overwrite: true, pattern: "*.log"

    input:
    set val(sample), val(single_end), file(fq) from ERROR_CORRECTION
    val genome_size from GS_ERROR_CORRECTION

    output:
    set val(sample), val(single_end), file('phix-r*.cor.fq') into REDUCE_COVERAGE
    file '*.log'

    shell:
    fq2 = single_end ? '' : '-r ' + fq[1]
    """
    lighter -od . -r !{fq[0]} !{fq2} -K 31 !{genome_size} -maxcor 1 \
            -t !{params.cpus} 2>&1 > lighter-ecc.log
    """
}

process reduce_coverage {
    /* Reduce the sequence coverage using BBDuk (reformat) */
    cpus 1
    tag "${sample}"

    input:
    set val(sample), val(single_end), file(fq) from REDUCE_COVERAGE
    val genome_size from GS_REDUCE_COVERAGE

    output:
    set val(sample), val(single_end), file('subsample-r*.fq') into FINAL_SUMMARY, MASH_SKETCH, FINISH_UP
    file '*.log'

    shell:
    total_bp = params.coverage.toInteger() * genome_size.toInteger()
    in2 = single_end ? '' : 'in2=' + fq[1]
    out2 = single_end ? '' : 'out2=subsample-r2.fq'
    """
    reformat.sh \
        in=!{fq[0]} !{in2} \
        out=subsample-r1.fq !{out2} \
        samplebasestarget=!{total_bp} \
        overwrite=t 2>&1 > reformat.log
    """
}

process final_summary {
    /* Run FASTQC on the processed FASTQ files. */
    cpus params.cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/sequence-qc/summary", mode: 'rellink', overwrite: true, pattern: '*.{html,json,zip}'

    input:
    set val(sample), val(single_end), file(fq) from FINAL_SUMMARY
    val genome_size from GS_FINAL_SUMMARY

    output:
    file '*.json'
    file '*fastqc.html'
    file '*fastqc.zip'

    shell:
    fq2 = single_end ? '' : fq[1]
    """
    cat !{fq[0]} !{fq2} > !{sample}-final
    cat !{sample}-final | fastq-scan !{genome_size} > fastq-scan-final.json
    fastqc --noextract -f fastq -t !{params.cpus} !{sample}-final
    """
}

process mash_sketch {
    /* Create a Mash Sketch of the processed FASTQ. */
    cpus 1
    tag "${sample}"
    publishDir "${outdir}/${sample}/sequence-qc", mode: 'rellink', overwrite: true, pattern: "*.msh"

    input:
    set val(sample), val(single_end), file(fq) from MASH_SKETCH

    output:
    file '*.msh'

    shell:
    fq2 = single_end ? '' : fq[1]
    """
    cat !{fq[0]} !{fq2} | \
    mash sketch -o !{sample} -s !{params.mash_s} -k !{params.mash_k} -m !{params.mash_m} -r -
    """
}

process finish_up {
    /* Compress the processed FASTQ files. */
    cpus 1
    tag "${sample}"
    publishDir "${outdir}/${sample}/sequence-qc", mode: 'rellink', overwrite: true, pattern: '*.fastq.gz'


    input:
    set val(sample), val(single_end), file(fq) from FINISH_UP

    output:
    file '*.fastq.gz'

    shell:
    if (single_end)
        """
        gzip --best -c -n !{fq[0]} > !{sample}.fastq.gz
        """
    else
        """
        gzip --best -c -n !{fq[0]} > !{sample}_R1.fastq.gz
        gzip --best -c -n !{fq[1]} > !{sample}_R2.fastq.gz
        """
}

/* ==== END FASTQ CLEANUP ==== */
workflow.onComplete {
    /*
    if (workflow.success == true && params.clear_cache_on_success) {
        // No need to resume completed run so remove cache.
        file('./work/').deleteDir()
        if (params.clear_logs) {
            file(logs_folder).deleteDir()
        }
    }
    */
    println """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    Launch Dir  : ${workflow.launchDir}
    Working Dir : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Command line: ${workflow.commandLine}
    Resumed?    : ${workflow.resume}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

// Utility Functions
def print_version() {
    println(NAME + ' ' + VERSION)
    exit 0
}

def print_usage() {
    log.info NAME + ' v' + VERSION
    log.info ''
    log.info 'Requirements: bbmap fastqc fastq-scan lighter mash nextflow'
    log.info ''
    log.info 'Required Options:'
    log.info '    --fastqs DIR        An input file containing the sample name and absolute'
    log.info '                            paths to FASTQs to process.'
    log.info ''
    log.info 'Optional Parameters:'
    log.info ''
    log.info '    --example_fastqs    Print an example of expected input for FASTQs file.'
    log.info ''
    log.info '    --check_fastqs      Verify "--fastqs" produces the expected output'
    log.info ''
    log.info '    --coverage INT      Reduce samples to a given coverage. (Default: 100)'
    log.info ''
    log.info '    --genome_size INT   Expected genome size (bp). (Default: mash estimate).'
    log.info ''
    log.info '    --outdir DIR        Directory to write results to. (Default ./)'
    log.info ''
    log.info '    --cpus INT          Number of processors to use for processing. '
    log.info '                           (Default: 1)'
    log.info ''
    log.info 'Mash Sketch Parameters:'
    log.info '    --mash_s INT        Sketch size. Each sketch will have at most this many '
    log.info '                            non-redundant min-hashes. (Default: 10000)'
    log.info ''
    log.info '    --mash_k INT        K-mer size. Hashes will be based on strings of this '
    log.info '                            many nucleotides. (Default: 31)'
    log.info ''
    log.info '    --mash_m INT        Minimum copies of each k-mer required to pass noise '
    log.info '                            filter for reads. (Default 3)'
    log.info ''
    log.info 'BBDuk Parameters:'
    log.info '    --adapters FASTA    Illumina adapters to remove (Default: BBmap adapters)'
    log.info ''
    log.info '    --adapter_k INT     Kmer length used for finding adapters. Adapters'
    log.info '                            shorter than k will not be found. (Default: 23)'
    log.info ''
    log.info '    --phix FASTA        phiX174 reference genome to remove (Default: NC_001422)'
    log.info ''
    log.info '    --phix_k INT        Kmer length used for finding phiX174. Contaminants'
    log.info '                            shorter than k will not be found. (Default: 31)'
    log.info ''
    log.info '    --ktrim STR         Trim reads to remove bases matching reference kmers.'
    log.info '                            Values:'
    log.info '                                f (do not trim),'
    log.info '                                r (trim to the right, Default),'
    log.info '                                l (trim to the left)'
    log.info ''
    log.info '    --mink INT          Look for shorter kmers at read tips down to this '
    log.info '                            length, when k-trimming or masking. 0 means '
    log.info '                            disabled. Enabling this will disable maskmiddle.'
    log.info '                            (Default: 11)'
    log.info ''
    log.info '    --hdist INT         Maximum Hamming distance for ref kmers (subs only).'
    log.info '                            Memory use is proportional to (3*K)^hdist.'
    log.info '                            (Default: 1)'
    log.info ''
    log.info '    --tpe BOOL          When kmer right-trimming, trim both reads to the '
    log.info '                            minimum length of either.'
    log.info '                            Values:'
    log.info '                                f (do not equally trim),'
    log.info '                                t (equally trim to the right, Default)'
    log.info ''
    log.info '    --tbo BOOL          Trim adapters based on where paired reads overlap.'
    log.info '                            Values:'
    log.info '                                f (do not trim by overlap),'
    log.info '                                t (trim by overlap, Default)'
    log.info ''
    log.info '    --qtrim STR         Trim read ends to remove bases with quality below '
    log.info '                            trimq. Performed AFTER looking for kmers. Values:'
    log.info '                                rl (trim both ends, Default),'
    log.info '                                f (neither end),'
    log.info '                                r (right end only),'
    log.info '                                l (left end only),'
    log.info '                                w (sliding window)'
    log.info ''
    log.info '    --trimq FLOAT       Regions with average quality BELOW this will be '
    log.info '                            trimmed if qtrim is set to something other than f.'
    log.info '                            (Default: 6)'
    log.info ''
    log.info '    --maq INT           Reads with average quality (after trimming) below '
    log.info '                            this will be discarded. (Default: 20)'
    log.info ''
    log.info '    --minlength INT     Reads shorter than this after trimming will be '
    log.info '                            discarded. Pairs will be discarded if both are '
    log.info '                            shorter. (Default: 35)'
    log.info ''
    log.info '    --ftm INT           If positive, right-trim length to be equal to zero, '
    log.info '                            modulo this number. (Default: 0)'
    log.info ''
    log.info '    --tossjunk          Discard reads with invalid characters as bases.'
    log.info '                            Values:'
    log.info '                                f (keep all reads),'
    log.info '                                t (toss reads with ambiguous bases, Default)'
    log.info ''
    log.info '    --qout STR          Output quality offset.'
    log.info '                            Values:'
    log.info '                                33 (PHRED33 offset quality scores, Default),'
    log.info '                                64 (PHRED64 offset quality scores)'
    log.info '                                auto (keeps the current input offset)'
    log.info ''
    log.info 'Lighter Parameters:'
    log.info '    --maxcor INT        Max number of corrections within a 20bp window '
    log.info '                            (Default: 1)'
    log.info ''
    log.info '    --version           Print version information'
    log.info '    --help              Show this message and exit'
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow illumina-cleanup.nf --fastqs my_fastqs.txt'

    exit 0
}

def print_example_fastqs() {
    log.info 'Printing example input for "--fastqs"'
    log.info ''
    log.info 'sample\tr1\tr2'
    log.info 'test001\t/path/to/fastqs/test_R1.fastq.gz\t/path/to/fastqs/test_R2.fastq.gz'
    log.info 'test002\t/path/to/fastqs/test.fastq.gz\t'

    exit 0
}

def check_input_params() {
    error = false
    missing_requirement = false
    if (!params.fastqs) {
        log.info('Missing required "--fastqs" input')
        error = true
    } else if (!file(params.fastqs).exists()) {
        log.info('Invalid input (--fastqs), please verify "' + params.fastqs + '"" exists.')
        error = true
    }

    if (error) {
        log.info('See --help for more information')
        exit 1
    }
}

def process_csv(line) {
    /* Parse line and determine if single end or paired reads*/
    if (line.r2) {
        // Paired
        return tuple(line.sample, false, [file(line.r1), file(line.r2)])
    } else {
        // Single End
        return tuple(line.sample, true, [file(line.r1)])
    }
}

def create_fastq_channel(fastq_input) {
    return Channel.fromPath( file(fastq_input) )
            .splitCsv(header: true, sep: '\t')
            .map { row -> process_csv(row) }
}

def check_input_fastqs(fastq_input) {
    /* Read through --fastqs and verify each input exists. */
    samples = [:]
    error = false
    file(fastq_input).splitEachLine('\t') { cols ->
        if (cols[0] != 'sample') {
            if (samples.containsKey(cols[0])) {
                samples[cols[0]] = samples[cols[0]] + 1
            } else {
                samples[cols[0]] = 1
            }
            if (cols[1]) {
                if (!file(cols[1]).exists()) {
                    log.info 'Error: Please verify ' + cols[1]+ ' exists, and try again'
                    error = true
                }
            }
            if (cols[2]) {
                if (!file(cols[2]).exists()) {
                    log.info 'Error: Please verify ' + cols[2]+ ' exists, and try again'
                    error = true
                }
            }
        }
    }

    samples.each{ sample, count ->
        if (count > 1) {
            error = true
            log.info 'Sample name "'+ sample +'" is not unique, please revise sample names'
        }
    }

    if (error) {
        log.info 'Verify sample names are unique and/or FASTQ paths are correct'
        log.info 'See "--example_fastqs" for an example'
        log.info 'Exiting'
        exit 1
    }
}
