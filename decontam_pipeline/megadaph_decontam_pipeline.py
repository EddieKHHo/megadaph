#!/usr/bin/env python
"""Megadaph Decontamination Pipeline

Built using Ruffus pipelining infrastructure (See http://www.ruffus.org.uk/)

Pipeline steps are as follows:
    Index Fasta Files -> Align Fastq files ->

Usage:
    megadaph_decontam_pipeline.py [--flowchart=<file.svg>] [-v | --verbose]
-p <nthreads> -a ASSEMBLY_DIR -r1 FORWARD_READ_DIR -r2 REVERSE_READ_DIR
    megadaph_decontam_pipeline.py (-h | --help)
    megadaph_decontam_pipeline.py (-j | --just-print)

Options:
    -h --help  Show this screen
    -j --just-print  Show actions that pipeline would take if run
    -f --flowchart=<file.svg>  Create a flowchart image describing the pipeline
    -v --verbose  Increased logging verbosity
    -p --threads  Number of threads to use
    -a --assemblies=ASSEMBLY_DIR  Location of directory holding .fasta assembly
files.
    -r1 --fwd_reads=FORWARD_READ_DIR  Location of directory holding forward
reads
    -r2 --rev_reads=REVERSE_READ_DIR  Location of directory holding reverse
reads. Can be the same as r1 - direction will be distinguished using .1/.2 in
filenames.


Required Input Files:
    1. A set of paired Fastq files
    2. A corresponding set of assemblies produced from the Fastq files in
    Fasta format

Notes:
    Matched Fastq and Fasta files must share the same prefix and be named using
    dot filename convention (E.g FA_SC.1.trim.fastq). Forward reads must have
    '.1' immediately after shared prefix, reverse reads must have '.2'.
"""

from docopt import docopt

if __name__ == '__main__':
    opt = docopt(__doc__)
    print(opt)
