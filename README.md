A Nextflow powered pipeline for cleaning up Illumina sequence reads.

# illumina-cleanup
*illumina-cleanup* is a simple Nextflow pipeline to cleanup Illumina FASTQ files. It is essentially a wrapper around [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) (BBDuk, Reformat), [Lighter](https://github.com/mourisl/Lighter) and [fastq-scan](https://github.com/rpetit3/fastq-scan).

# Installation
### Bioconda


### From Souce
# Dependencies
In order to use *illumina-cleanup* you will need to have the following programs installed. The program versions that were tested are given in parentheses.

* __Nextflow (v19.01.0)__  
Used to manage the workflow.  
_Di Tommaso, P., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E., Notredame, C., 2017. [Nextflow enables reproducible computational workflows.](https://www.nature.com/articles/nbt.3820.pdf?origin=ppub) Nat. Biotechnol. 35, 316–319._

* __BBTools (v38.34)__  
Used *BBDuk* and *Reformat* for removing adapters, PhiX contaminants, quality (PHRED-based) filtering, and randomly subsampling.  
_[BBTools Home Page](https://jgi.doe.gov/data-and-tools/bbtools/)_

* __FastQC (v0.11.8)__  
Used to generate visuals of sequencing quality and as input for [MultiQC](https://multiqc.info/)  
_Andrews, S. FastQC: a quality control tool for high throughput sequence data.
(http://www.bioinformatics.babraham.ac.uk/projects/fastqc)._

* __fastq-scan (v0.3)__  
Used to get a quick summary of sequencing quality in JSON format.  
_[fastq-scan Home Page](https://github.com/rpetit3/fastq-scan)_

* __Lighter (v1.1.2)__  
Used for sequence error correction.   
_Song, L., Florea, L. and Langmead, B., [Lighter: Fast and Memory-efficient Sequencing Error Correction without Counting](http://genomebiology.com/2014/15/11/509/). Genome Biol. 2014 Nov 15;15(11):509._

* __pigz (v2.3.4)__  
Used to speed up G-zip steps.  
_Adler, Mark. "pigz: A parallel implementation of gzip for modern multi-processor, multi-core machines." Jet Propulsion Laboratory (2015)._

```
illumina-cleanup
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [distracted_bartik] - revision: 0416ba407c

illumina-cleanup v1.0.0
Usage: illumina-cleanup --fastqs input_fastqs.txt --max_cpus 8

Requirements: bbmap fastqc fastq-scan lighter nextflow pigz

Required Parameters:
    --fastqs STR        An input file containing the sample name and absolute
                            paths to FASTQs to process.

    --max_cpus INT       The maximum number of processors this workflow should
                             have access to at any given moment. The default
                             for Nextflow is to use all available processors.

    --example_fastqs    Print an example of expected input for FASTQs file.

    --check_fastqs      Verify "--fastqs" produces the expected inputs.

     -name              Name for the pipeline run. If not specified, Nextflow
                            will automatically generate a random mnemonic.
Optional Parameters:

    --coverage INT      Reduce samples to a given coverage. (Default: 100)

    --genome_size INT   Expected genome size (bp) for all samples. (Default: 0).

    --outdir DIR        Directory to write results to. (Default ./)

    --cpus INT          Number of processors made available to a single process.
                            If greater than "--max_cpus" it will be set equal to
                            "--max_cpus" (Default: 1)

BBDuk Parameters:
    --adapters FASTA    Illumina adapters to remove (Default: BBmap adapters)

    --adapter_k INT     Kmer length used for finding adapters. Adapters
                            shorter than k will not be found. (Default: 23)

    --phix FASTA        phiX174 reference genome to remove (Default: NC_001422)

    --phix_k INT        Kmer length used for finding phiX174. Contaminants
                            shorter than k will not be found. (Default: 31)

    --ktrim STR         Trim reads to remove bases matching reference kmers.
                            Values:
                                f (do not trim)
                                r (trim to the right, Default)
                                l (trim to the left)

    --mink INT          Look for shorter kmers at read tips down to this
                            length, when k-trimming or masking. 0 means
                            disabled. Enabling this will disable maskmiddle.
                            (Default: 11)

    --hdist INT         Maximum Hamming distance for ref kmers (subs only).
                            Memory use is proportional to (3*K)^hdist.
                            (Default: 1)

    --tpe BOOL          When kmer right-trimming, trim both reads to the
                            minimum length of either.
                            Values:
                                f (do not equally trim)
                                t (equally trim to the right, Default)

    --tbo BOOL          Trim adapters based on where paired reads overlap.
                            Values:
                                f (do not trim by overlap)
                                t (trim by overlap, Default)

    --qtrim STR         Trim read ends to remove bases with quality below
                            trimq. Performed AFTER looking for kmers. Values:
                                rl (trim both ends, Default)
                                f (neither end)
                                r (right end only)
                                l (left end only)
                                w (sliding window)

    --trimq FLOAT       Regions with average quality BELOW this will be
                            trimmed if qtrim is set to something other than f.
                            (Default: 6)

    --maq INT           Reads with average quality (after trimming) below
                            this will be discarded. (Default: 20)

    --minlength INT     Reads shorter than this after trimming will be
                            discarded. Pairs will be discarded if both are
                            shorter. (Default: 35)

    --ftm INT           If positive, right-trim length to be equal to zero,
                            modulo this number. (Default: 0)

    --tossjunk          Discard reads with invalid characters as bases.
                            Values:
                                f (keep all reads)
                                t (toss reads with ambiguous bases, Default)

    --qout STR          Output quality offset.
                            Values:
                                33 (PHRED33 offset quality scores, Default)
                                64 (PHRED64 offset quality scores)
                                auto (keeps the current input offset)

    --xmx STR           This will be passed to Java to set memory usage.
                            Examples:
                                8g will specify 8 gigs of RAM (Default)
                                20g will specify 20 gigs of RAM
                                200m will specify 200 megs of RAM

Lighter Parameters:
    --maxcor INT        Max number of corrections within a 20bp window
                            (Default: 1)

Reformat (BBmap) Parameters:
    --sampleseed INT    Set to a positive number to use that prng seed for
                            sampling (Default 42).

    --keep_cache        Keeps 'work' and '.nextflow' logs, default is to
                            delete on successful completion.
    --version           Print workflow version information
    --help              Show this message and exit
    
illumina-cleanup --version
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [amazing_gilbert] - revision: 0416ba407c
illumina-cleanup 1.0.0

illumina-cleanup --example_fastqs
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [naughty_borg] - revision: 0416ba407c
Printing example input for "--fastqs"

sample  r1      r2
test001 /path/to/fastqs/test_R1.fastq.gz        /path/to/fastqs/test_R2.fastq.gz
test002 /path/to/fastqs/test.fastq.gz

illumina-cleanup --check_fastqs --fastqs example-data/bad-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [lonely_hypatia] - revision: 96c6a1a7ae
4:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test003_R1.fastq.gz exists, and try again
5:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test003_R1.fastq.gz exists, and try again
5:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test002_R2.fastq.gz exists, and try again
Sample name "test002" is not unique, please revise sample names
Verify sample names are unique and/or FASTQ paths are correct
See "--example_fastqs" for an example
Exiting

illumina-cleanup --check_fastqs --fastqs example-data/good-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [astonishing_colden] - revision: 96c6a1a7ae
Printing what would have been processed. Each line consists of an array of
three elements: [SAMPLE_NAME, IS_SINGLE_END, [FASTQ_1, FASTQ_2]]

Found:
[test001, false, [/home/rpetit/illumina-cleanup/test/fastqs/test_R1.fastq.gz, /home/rpetit/illumina-cleanup/test/fastqs/test_R2.fastq.gz]]
[test002, true, [/home/rpetit/illumina-cleanup/test/fastqs/test.fastq.gz]]
```
