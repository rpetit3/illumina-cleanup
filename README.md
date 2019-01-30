A Nextflow powered pipeline for pre-processing Illumina sequence reads.

# illumina-cleanup
*illumina-cleanup* is a simple Nextflow pipeline to remove Illumina adapters, PhiX contaminants, quality filter, base-error correction, and reduced sequence coverage Illumina samples. This is accomplished by making use of [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) (BBDuk, Reformat), [Lighter](https://github.com/mourisl/Lighter) and [fastq-scan](https://github.com/rpetit3/fastq-scan). 

# Installation
### Bioconda
In the works!

### From Souce
The first requirement is to have each of the dependencies (see [*Dependencies*](https://github.com/rpetit3/illumina-cleanup/tree/version-update#dependencies)) installed and available in your PATH.

With all dependencies installed:
```
git clone git@github.com:rpetit3/illumina-cleanup.git
cd illumina-cleanup
bin/illumina-cleanup --version
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [amazing_gilbert] - revision: 0416ba407c
illumina-cleanup 1.0.0
```

# Usage
Since *illumina-cleanup* is essentially a wrapper, many of the runtime parameters are for individual tools. Although, there a few important parameters specific to *illumina-cleanup* which are discussed below.

### `--help` Output
```
./illumina-cleanup
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [distracted_bartik] - revision: 0416ba407c

illumina-cleanup v1.0.0
Usage: illumina-cleanup --fastqs input_fastqs.txt --max_cpus 8

Dependencies: bbmap fastqc fastq-scan lighter nextflow pigz

Required Parameters:
    --fastqs STR        An input file containing the sample name and absolute
                            paths to FASTQs to process.

Optional Parameters:
    --coverage INT      Reduce samples to a given coverage. (Default: 100)

    --genome_size INT   Expected genome size (bp) for all samples. (Default: 0).

    --outdir DIR        Directory to write results to. (Default ./)

    --max_cpus INT       The maximum number of processors this workflow should
                             have access to at any given moment. (Default: 1)

    --cpus INT          Number of processors made available to a single process.
                            If greater than "--max_cpus" it will be set equal to
                            "--max_cpus" (Default: 1)

    --example_fastqs    Print an example of expected input for FASTQs file.

    --check_fastqs      Verify "--fastqs" produces the expected inputs.

    --keep_cache        Keeps 'work' and '.nextflow' logs, default is to
                            delete on successful completion.

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

    --version           Print workflow version information
    --help              Show this message and exit
```

### `--fastqs` a *File of Filenames* for Input FASTQs
*illumina-cleanup* was developed to take in a file with information about the inputs, a *file of filenames* (FOFN). Essentially what this file does is specify the location of FASTQs to be processed, whether they are single-end (SE) or paired-end (PE) reads, and what to name (sample name) the output files. 

While this is an additional step for you, the user, it avoids potenial bugs associated with pattern matching FASTQs to extract a name and the SE vs PE status. Most importantly, by taking this approach, Nextflow can easily queue the analysis of a single FASTQ or hundreds of FASTQs in a single command. There is also the added benefit of knowing which FASTQs were analysed and their location at a later time.

#### Structure of the *FASTQs* FOFN
You can use the `--example_fastqs` to get an example of the expected structure for the input FASTQs FOFN.

```
illumina-cleanup --example_fastqs
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [naughty_borg] - revision: 0416ba407c
Printing example input for "--fastqs"

sample  r1      r2
test001 /path/to/fastqs/test_R1.fastq.gz        /path/to/fastqs/test_R2.fastq.gz
test002 /path/to/fastqs/test.fastq.gz
```

The expected structure is a **tab-delimited** table with three columns:

1. `sample`: A unique prefix, or unique name, to be used for naming output files
2. `r1`: If paired-end, the first pair of reads, else the single-end reads
3. `r2`: If paired-end, the second pair of reads

These three columns are used as the header for the file. In other words, all input FOFNs require their first line to be:
```
sample  r1      r2
```

All lines after the header line, contain unique sample names and location(s) to associated FASTQ file(s). For the purpose of the FOFN, absolute paths are preferred as they are more informative in the future. Relative paths may be ok, but have not been thoroughly tested. 

In the example above, two samples would be cleaned up. Sample `test001` has two FASTQs and would be processed as pair-end reads. While sample `test002` only has a single FASTQ and would be processed as single-end reads.

#### Validate Your FASTQ FOFN
When a FOFN is given, the first thing *illumina-cleanup* does is verify the path all FASTQ files is valid. If all paths are valid, each sample will then be processed, otherwise a list of samples with errors will be output to STDERR. 

If you would like to only validate your FOFN (and not run the full pipeline), you can use the `--check_fastqs` parameter.

##### FOFN without Errors
```
illumina-cleanup --check_fastqs --fastqs example-data/good-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [astonishing_colden] - revision: 96c6a1a7ae
Printing what would have been processed. Each line consists of an array of
three elements: [SAMPLE_NAME, IS_SINGLE_END, [FASTQ_1, FASTQ_2]]

Found:
[test001, false, [/home/rpetit/illumina-cleanup/test/fastqs/test_R1.fastq.gz, /home/rpetit/illumina-cleanup/test/fastqs/test_R2.fastq.gz]]
[test002, true, [/home/rpetit/illumina-cleanup/test/fastqs/test.fastq.gz]]
```
Each sample has passed validation and is put into a three element array:

1. sample - the name for this sample
2. is_single_end - the reads are single-end (true) or paired-end (false)
3. fastq_array - the fastqs associated with the sample

This array is then automatically queued up for proccessing by Nextflow.

##### FOFN with errors
```
illumina-cleanup --check_fastqs --fastqs example-data/bad-fastqs.txt
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [kickass_mestorf] - revision: 222a5ad8b1
LINE 4:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test003_R1.fastq.gz exists, and try again
LINE 5:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test003_R1.fastq.gz exists, and try again
LINE 5:ERROR: Please verify /home/rpetit/illumina-cleanup/test/fastqs/test002_R2.fastq.gz exists, and try again
Sample name "test002" is not unique, please revise sample names
The header line (line 1) does not follow expected structure.
Verify sample names are unique and/or FASTQ paths are correct
See "--example_fastqs" for an example
Exiting
```

In the above example, there are mulitple errors. Lines 4 and 5 (`LINE 4:ERROR` or `LINE 5:ERROR`) suggest that based on the given paths the FASTQs do not exist. The sample name `test002` has been used multiple times, and must be corrected. There is also an issue with the header line that must be looked into.


### `--max_cpus` vs `--cpus`
By default when Nextflow executes, it uses all available cpus are used to queue up processes. As you might imagine, if you are on a single server with multiple users, this approach (use all cpus) might annoy other users. To circumvent this feature, two parmeters have been included `--max_cpus` and `--cpus`.


```
    --max_cpus INT       The maximum number of processors this workflow should
                             have access to at any given moment. The default
                             for Nextflow is to use all available processors.

    --cpus INT          Number of processors made available to a single process.
                            If greater than "--max_cpus" it will be set equal to
                            "--max_cpus" (Default: 1)
 ```
 
What `--max_cpus` does is specify to Nextflow the maximum number of cpus it is allowed to occupy at any given time. `--cpus` on the other hand, specifies how many cpus any given step (adapter trimming, qc, etc...) can occupy. By default **`--max_cpus` is set to 1** and if `--cpus` is set to a value greater than `--max_cpus` it will be set equal to `--max_cpus`. This appoach errs on the side of caution (e.g. not taking up the whole server!) and lets Nextflow figure out how best to queue up tasks.

### `--genome_size` for Error Correction and Coverage Reduction
Base error corrections are handled by *Lighter*. Lighter requires an expected genome size to be given at runtime. If a genome size is not given by the user, the error correction step will be skipped.

Also, a genome size is required to randomly (*reformat.sh* from BBTools) subsample sequences to a given coverage. Again this step will be skipped if a genome size is not given. 

### `--keep_cache` for Intermediate Files
By default, *illumina-cleanup* will remove all intermediate files **after successfully completing all steps**. It does this by removing the `work` and `.nextflow` directories created by Nextflow. Again, this only occurs after all samples have been successfully processed. If there is an error during execution, the cache will remain intact and the process can still be resumed (`-resume`). 

If you would like to retain the cache files, you can use the `--keep_cache` parameters to do so. Please keep in mind there is storage overhead in doing so. The cache will contain multiple uncompressed FASTQ files for each sample that was processed.

### Remaining Parameters
The remaining parameters are mostly associated with specific programs. These parameters are grouped by program in the usage. All available parameters from a given program may not be available in *illumina-cleanup* (Example: BBDUk's extensive set of parameters). For suggested use cases, these parameters may be made available in the future.

# Dependencies
In order to use *illumina-cleanup* you will need to have the following programs installed. The program versions that were tested are given in parentheses.

* __Nextflow (v19.01.0)__  
Used to manage the workflow.  
_Di Tommaso, P., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E., Notredame, C., 2017. [Nextflow enables reproducible computational workflows.](https://www.nature.com/articles/nbt.3820.pdf?origin=ppub) Nat. Biotechnol. 35, 316–319._

* __BBTools (v38.34)__  
Used *BBDuk* and *Reformat* for removing adapters, PhiX contaminants, quality (PHRED-based) filtering, and randomly subsampling.  
_[BBTools Home Page](https://jgi.doe.gov/data-and-tools/bbtools/)_

* __FastQC (v0.11.8)__  
Used to generate visuals of sequence quality and as input for [MultiQC](https://multiqc.info/)  
_Andrews, S. FastQC: a quality control tool for high throughput sequence data.
(http://www.bioinformatics.babraham.ac.uk/projects/fastqc)._

* __fastq-scan (v0.3)__  
Used to get a quick summary of sequence quality in JSON format.  
_[fastq-scan Home Page](https://github.com/rpetit3/fastq-scan)_

* __Lighter (v1.1.2)__  
Used for sequence error correction.   
_Song, L., Florea, L. and Langmead, B., [Lighter: Fast and Memory-efficient Sequencing Error Correction without Counting](http://genomebiology.com/2014/15/11/509/). Genome Biol. 2014 Nov 15;15(11):509._

* __pigz (v2.3.4)__  
Used to speed up G-zip steps.  
_Adler, Mark. "pigz: A parallel implementation of gzip for modern multi-processor, multi-core machines." Jet Propulsion Laboratory (2015)._

# Version
```  
illumina-cleanup --version
N E X T F L O W  ~  version 19.01.0
Launching `/home/rpetit/illumina-cleanup/bin/illumina-cleanup` [amazing_gilbert] - revision: 0416ba407c
illumina-cleanup 1.0.0
```

# Disclaimer
I have only tested this against bacterial sequences. Please verify the default parameters given are valid for your use case. 
