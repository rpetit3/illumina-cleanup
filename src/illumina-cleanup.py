#! /usr/bin/env python
"""
Filter low quality reads, and reduce the read count to a given coverage.
Input is read from STDIN, cleaned up, and printed to STDOUT.
usage: fastq_cleanup [--coverage FLOAT] [--genome_size INT]
                 [--total_read_count INT] [--read_length_cutoff INT]
                 [--paired_reads] [--min_mean_quality INT]
                 [--min_read_length INT] [-i STRING] [-h] [--version]
Example Usage: zcat SOME.fastq.gz | fastq_cleanup --total_read_count 836846
"""
import random
import numpy as np


class IlluminaQC(object):
    """Class for the clean up of FASTQ files."""

    def __init__(self, subsample, is_paired, read_length_cutoff,
                 min_mean_quality, min_read_length):
        """Initialize class variables."""
        self.VERSION = "0.5"
        self.fastq = []
        self.subsample = subsample
        self.is_paired = is_paired
        self.read_length_cutoff = read_length_cutoff
        self.min_mean_quality = min_mean_quality
        self.min_read_length = min_read_length
        self.__temp_read = []
        self.__phred64 = 0
        self.__phred33 = 0
        self.__phredunk = 0
        self.__contains_ns = 0
        self.__total_ns = 0
        self.__read_length = 0
        self.__bad_quality = 0
        self.__n_lengths = 0
        self.__missed_ns = 0

    def read_large_fastq(self, file, fraction):
        """Reduce a large (2 GB+) file to a subset before cleanup."""
        random.seed(123456)
        while 1:
            head = file.readline().rstrip()
            if not head:
                break

            seq = file.readline().rstrip()
            plus = file.readline().rstrip()
            qual = file.readline().rstrip()

            if random.random() <= fraction:
                self.fastq.append(head)
                self.fastq.append(seq)
                self.fastq.append(plus)
                self.fastq.append(qual)

                if self.is_paired:
                    self.fastq.append(file.readline().rstrip())
                    self.fastq.append(file.readline().rstrip())
                    self.fastq.append(file.readline().rstrip())
                    self.fastq.append(file.readline().rstrip())

    def __mean_quality(self, qual):
        """Create a count of the quality score."""
        qual_stats = np.array([ord(j) for j in qual])
        return np.mean(qual_stats) - 33

    def __get_nonambiguous_sequence(self, seq, qual):
        """Return the longest non-ambiguous sequence."""
        seqs = []
        new_seq = []
        quals = []
        new_qual = []
        for index, base in enumerate(seq):
            if base == 'N' and len(new_seq):
                seqs.append(''.join(new_seq))
                quals.append(''.join(new_qual))
                new_seq = []
                new_qual = []
            else:
                new_seq.append(base)
                new_qual.append(qual[index])

        if len(new_seq):
            seqs.append(''.join(new_seq))
            quals.append(''.join(new_qual))

        # Get the largest chunk
        new_seq = ''
        new_qual = ''
        for index, chunk in enumerate(seqs):
            if len(chunk) > len(new_seq):
                new_seq = chunk
                new_qual = quals[index]
            elif len(chunk) == len(new_seq):
                random.seed(123456)
                if random.random() <= 0.5:
                    new_seq = chunk
                    new_qual = quals[index]

        return [new_seq, new_qual]

    def __quality_trim(self, seq, qual):
        """Remove bases below minimum quality."""
        new_seq = []
        new_qual = []
        for index, base in enumerate(seq):
            q = ord(qual[index]) - 33
            if q < self.min_mean_quality and index >= self.min_read_length:
                break
            else:
                new_seq.append(base)
                new_qual.append(qual[index])

        return [''.join(new_seq), ''.join(new_qual)]

    def __test_read(self, index, append=''):
        """Test if the read passes quality filters."""
        head = self.fastq[index].split()[0]
        if append:
            if not head.endswith(append):
                head = "{0}{1}".format(self.fastq[index].split()[0], append)

        seq = self.fastq[index + 1]
        length = 0
        qual = self.fastq[index + 3]

        # If ambiguous nucleotides, get largest subsequence without
        # ambiguous nucleotides.
        if seq.count('N') > 0:
            seq, qual = self.__get_nonambiguous_sequence(seq, qual)
            self.__n_lengths += len(seq)
            self.__total_ns += seq.count('N')
            self.__contains_ns += 1

        if (seq.count('N') == 0 and len(seq) >= self.min_read_length):
            if (self.read_length_cutoff and length > self.read_length_cutoff):
                seq = seq[:self.read_length_cutoff]
                qual = qual[:self.read_length_cutoff]
            else:
                if (self.read_length_cutoff):
                    self.__read_length += 1

            if (self.__mean_quality(qual) >= self.min_mean_quality):
                self.__temp_read.append('{0}\n{1}\n+\n{2}'.format(
                    head,
                    seq,
                    qual
                ))
                length = len(seq)
            else:
                self.__bad_quality += 1
        else:
            self.__missed_ns += 1

        return length

    def generate_order(self, total_read_count):
        """Generate a random (seeded) order of reads to be selected."""
        if self.is_paired and total_read_count % 2 == 0:
            total_read_count = total_read_count / 2
        self.__read_order = list(range(int(total_read_count)))
        if self.subsample:
            random.seed(123456)
            random.shuffle(self.__read_order)

    def clean_up_fastq(self):
        """Test each read, if good print it else move next read."""
        basepair_count = 0
        for i in range(len(self.__read_order)):
            self.__temp_read[:] = []
            if self.is_paired:
                index = self.__read_order[i] * 8 - 8
                read1 = self.__test_read(index, append="/1")
                index = self.__read_order[i] * 8 - 4
                read2 = self.__test_read(index, append="/2")
                if read1 > 0 and read2 > 0:
                    print('\n'.join(self.__temp_read))
                    basepair_count += (read1 + read2)
            else:
                index = self.__read_order[i] * 4 - 4
                read = self.__test_read(index)
                if read > 0:
                    print(self.__temp_read[0])
                    basepair_count += read

            if self.subsample and basepair_count >= self.subsample:
                break


def define_min_quality(stats, cutoff, paired):
    """Determine minimum quality score to filter reads on."""
    # Make it a little more relaxed for paired reads
    x = 2 if paired else 1
    if stats['coverage'] >= (cutoff * 3 * x):
        return int(stats['qual_75th'])
    elif stats['coverage'] >= (cutoff * 2 * x):
        return int(stats['qual_median'])
    elif stats['coverage'] >= (cutoff * 1.5 * x):
        return int(stats['qual_25th'])
    elif stats['qual_mean'] - 2 * stats['qual_std'] >= 20:
        return int(stats['qual_mean'] - 2 * stats['qual_std'])
    else:
        return 20


def define_min_read(stats, cutoff, paired):
    """Determine minimum read length to filter reads on."""
    # Make it a little more relaxed for paired reads
    x = 2 if paired else 1
    if cutoff * 3 * x:
        return int(stats['read_median'])
    elif stats['coverage'] >= (cutoff * 2 * x):
        return int(stats['read_25th'])
    elif stats['coverage'] >= (cutoff * 1.5 * x):
        return int(stats['read_mean'] - (2 * stats['read_std']))
    elif stats['read_mean'] - (2 * stats['read_std']) >= 70:
        return 70
    else:
        return 35


if __name__ == '__main__':
    import sys
    import json
    import argparse as ap

    parser = ap.ArgumentParser(
        prog='illumina-cleanup.py',
        conflict_handler='resolve',
        description=('Filter low quality reads and reduce the read count to a '
                     'given coverage. Input read from STDIN, cleaned up reads '
                     'are print to STDOUT.'))
    group1 = parser.add_argument_group('Options', '')
    group1.add_argument('--stats', metavar="STR", type=str,
                        help='Read statisitcs from fastq_stats',
                        default=False)
    group1.add_argument('--coverage', metavar="INT", type=int,
                        help='Subsample coverage. (Default: Take all)',
                        default=100)
    group1.add_argument('--genome_size', metavar="INT", type=int,
                        help='Estimated genome size. (Default: 2814816)',
                        default=2814816)
    group1.add_argument('--total_read_count', metavar="INT", type=int,
                        help='Total count of input reads.')
    group1.add_argument('--read_length_cutoff', metavar="INT", type=int,
                        help='Trim reads to a certain length.', default=False)
    group1.add_argument('--paired', action='store_true', default=False,
                        help='Input is interleaved paired end reads.', )
    group1.add_argument('--no_length_filter', action='store_true',
                        default=False,
                        help='Do not filter reads based on read lengths.', )
    group1.add_argument('--min_mean_quality', metavar="INT", type=int,
                        help='Minimum mean read quality cutoff. (Default: 20)',
                        default=20)
    group1.add_argument('--min_read_length', metavar="INT", type=int,
                        help='Minimum read length cutoff. (Default: 35bp)',
                        default=35)

    group2 = parser.add_argument_group('Extra', '')
    group2.add_argument('-i', '--jobid', help='Job ID of sequence',
                        metavar="STRING")

    group3 = parser.add_argument_group('Optional', '')
    group3.add_argument('-h', '--help', action='help',
                        help='Show this help message and exit')
    group3.add_argument('--version', action='version', version='%(prog)s v0.1',
                        help='Show program\'s version number and exit')

    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit(1)

    args = parser.parse_args()

    # Read JSON output of fastq_stats
    stats = None
    if args.stats:
        with open(args.stats, 'r') as f:
            json_data = json.load(f)

        stats = json_data["qc_stats"]
        args.min_mean_quality = define_min_quality(
            stats, args.coverage, args.paired
        )
        args.total_read_count = stats['read_total']
        args.min_read_length = stats['read_min']

        if not args.no_length_filter:
            suggested_min = define_min_read(stats, args.coverage, args.paired)
            if suggested_min <= stats['read_max']:
                args.min_read_length = suggested_min

    subsample = args.coverage * args.genome_size if args.coverage else False

    # Process FASTQ
    fq = IlluminaQC(subsample, args.paired,
                    args.read_length_cutoff, args.min_mean_quality,
                    args.min_read_length)

    # If large fastq reduce to random subset of 2.5x coverage cutoff
    if stats['coverage'] > 750:
        fraction = (args.genome_size * 750) / stats['total_bp']
        fq.read_large_fastq(sys.stdin, fraction)
        args.total_read_count = len(fq.fastq) / 4
    else:
        fq.fastq = [line.rstrip() for line in sys.stdin.readlines()]

    fq.generate_order(args.total_read_count)
    fq.clean_up_fastq()
