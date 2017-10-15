#!/usr/bin/env python
"""Megadaph Decontamination Pipeline

Built using Ruffus pipelining infrastructure (See http://www.ruffus.org.uk/)

Pipeline steps are as follows:
    Index Fasta Files -> Align Fastq files ->

Usage:
    megadaph_decontam_pipeline.py [--flowchart=<file.svg>] [--verbose] \
    --threads=<n> --assemblies=<directory> --fwd_reads=<directory> \
    --rev_reads=<directory>
    megadaph_decontam_pipeline.py (-h | --help)
    megadaph_decontam_pipeline.py (-j | --just-print)

Options:
    -h --help  Show this screen
    -j --just-print  Show actions that pipeline would take if run
    -c --flowchart=<file.svg>  Create a flowchart image describing the pipeline
    -v --verbose  Increased logging verbosity
    -p --threads=<nthreads>  Number of threads to use
    -a --assemblies=ASSEMBLY_DIR  Location of directory holding .fasta assembly
files.
    -x --fwd_reads=FORWARD_READ_DIR  Location of directory holding forward
reads
    -y --rev_reads=REVERSE_READ_DIR  Location of directory holding reverse
reads. Can be the same as x - direction will be distinguished using 1./2. in
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

from pathlib import Path
from typing import Dict
from typing import List
from typing import Tuple

from docopt import docopt
import ruffus

import fmbiopy.biofile as biofile
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus


def _get_gzip_status(files: List[Path]) -> bool:
    """Get the gzip status of a group of files

    Files are expected to have the same Gzip status.

    Returns
    -------
    True if all `files` are gzipped, False if all `files` are unzipped.

    Raises
    ------
    ValueError
        If only some files are gzipped
    """
    gzipped = [f.suffix == '.gz' for f in files]
    unzipped = [not x for x in gzipped]
    if all(gzipped):
        return True
    elif all(unzipped):
        return False
    raise ValueError("Only some of the input files are gzipped")


def _validate_args(
        assemblies: List[Path],
        fwd_reads: List[Path],
        rev_reads: List[Path],
        )-> bool:
    """Check that the given input files are valid

    Paired fastq files must share the same prefix and use .1/.2 to distinguish
    read direction.

    Parameters
    ----------
    fwd_reads, rev_reads
        The user inputted fastq file paths
    assemblies
        The user inputted assembly files

    Returns
    -------
    True if successful.
    """

    inputs = [fwd_reads, rev_reads, assemblies]
    gzipped = [_get_gzip_status(input) for input in inputs]
    biofile.MatchedPrefixGroup([
        biofile.BiofileGroup(fwd_reads, filetype='fastq', gzipped=gzipped[0]),
        biofile.BiofileGroup(rev_reads, filetype='fastq', gzipped=gzipped[1]),
        biofile.BiofileGroup(assemblies, filetype='fasta', gzipped=gzipped[2])
        ]).validate()
    return True


def _extract_input_files(opt: Dict[str, str]) -> List[Tuple[Path, Path, Path]]:
    """Process command line arguments

    Input file prefixes are checked to ensure they can be matched correctly and
    have the correct file extensions.

    Returns
    -------
    Dict[str, List[str]]
        A dictionary of the given input files
    """

    # Gather the inputs
    inputs: List[str] = [
            opt['--assemblies'], opt['--fwd_reads'], opt['--rev_reads']]
    # Convert to Paths
    paths = [Path(direc).resolve() for direc in inputs]
    # The expected filetypes
    extensions = [
            biofile.Fasta.accepted_extensions,
            biofile.Fastq.accepted_extensions,
            biofile.Fastq.accepted_extensions]
    # The expected substrings
    substrings = [None, '1.', '2.']

    # Get the paths to all input files
    file_paths: List[List[Path]] = []
    for path, typ, substring in zip(paths, extensions, substrings):
        file_paths.append(fmpaths.find(path, typ, substring))

    assemblies = file_paths[0]
    fwd_reads = file_paths[1]
    rev_reads = file_paths[2]
    _validate_args(assemblies, fwd_reads, rev_reads)

    return list(zip(assemblies, fwd_reads, rev_reads))


def _main() -> bool:
    """Main function

    Returns
    -------
    True if runs without error
    """
    opt = docopt(__doc__)
    input_files = _extract_input_files(opt)

    # Get suffix of the fasta files, so that we can determine the output name
    fasta_suffix = input_files[0][0].suffix

    # Pipeline start
    pipe = ruffus.Pipeline(name="decontam")
    pipe.transform(
            task_func=fmruffus.PairedBowtie2Align,
            name='alignment',
            filter=ruffus.suffix(fasta_suffix),
            input=input_files,
            output='.bam',
            output_dir='PairedBowtie2align')
    pipe.run()
    return True


if __name__ == '__main__':
    _main()
